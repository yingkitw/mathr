//! Minimal HTTP server for the math notebook web UI.
//!
//! Serves a single-page app at `GET /` and provides a JSON API:
//! - `POST /api/eval` — evaluate a cell expression
//! - `GET /api/notebook` — return current notebook JSON
//! - `POST /api/notebook` — replace notebook JSON

use crate::error::{MathError, Result};
use crate::eval::Context;
use crate::notebook::Notebook;
use crate::repl;
use std::io::{BufRead, BufReader, Read, Write};
use std::net::{TcpListener, TcpStream};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

/// Shared notebook state.
pub struct NotebookServer {
    notebook: Arc<Mutex<Notebook>>,
    file_path: Option<PathBuf>,
}

impl NotebookServer {
    pub fn new(notebook: Notebook, file_path: Option<PathBuf>) -> Self {
        Self {
            notebook: Arc::new(Mutex::new(notebook)),
            file_path,
        }
    }

    /// Start serving on the given port. Blocks until Ctrl-C.
    pub fn serve(&self, port: u16) -> Result<()> {
        let listener = TcpListener::bind(format!("127.0.0.1:{}", port)).map_err(|e| {
            MathError::InvalidArgument(format!("cannot bind port {}: {}", port, e))
        })?;
        eprintln!("mathr notebook server: http://127.0.0.1:{}", port);
        eprintln!("Press Ctrl-C to stop.");

        let notebook = self.notebook.clone();
        let file_path = self.file_path.clone();

        for stream in listener.incoming() {
            let stream = match stream {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("connection error: {}", e);
                    continue;
                }
            };
            let nb = notebook.clone();
            let fp = file_path.clone();
            std::thread::spawn(move || {
                if let Err(e) = handle_connection(stream, nb, fp) {
                    eprintln!("request error: {}", e);
                }
            });
        }
        Ok(())
    }
}

fn handle_connection(
    stream: TcpStream,
    notebook: Arc<Mutex<Notebook>>,
    file_path: Option<PathBuf>,
) -> Result<()> {
    let mut reader = BufReader::new(stream.try_clone().map_err(|e| {
        MathError::InvalidArgument(format!("io error: {}", e))
    })?);
    let mut writer = stream;

    // Read request line
    let mut request_line = String::new();
    reader
        .read_line(&mut request_line)
        .map_err(|e| MathError::InvalidArgument(format!("io error: {}", e)))?;
    let request_line = request_line.trim().to_string();

    // Read headers
    let mut content_length = 0usize;
    loop {
        let mut header = String::new();
        let n = reader
            .read_line(&mut header)
            .map_err(|e| MathError::InvalidArgument(format!("io error: {}", e)))?;
        if n == 0 || header.trim().is_empty() {
            break;
        }
        let lower = header.to_lowercase();
        if lower.starts_with("content-length:") {
            content_length = lower["content-length:".len()..]
                .trim()
                .parse()
                .unwrap_or(0);
        }
    }

    // Read body if present
    let mut body = String::new();
    if content_length > 0 {
        let mut buf = vec![0u8; content_length];
        reader
            .read_exact(&mut buf)
            .map_err(|e| MathError::InvalidArgument(format!("io error: {}", e)))?;
        body = String::from_utf8_lossy(&buf).to_string();
    }

    // Route
    let parts: Vec<&str> = request_line.split_whitespace().collect();
    if parts.len() < 2 {
        send_response(&mut writer, 400, "text/plain", "Bad Request")?;
        return Ok(());
    }
    let method = parts[0];
    let path = parts[1];

    match (method, path) {
        ("GET", "/") => {
            send_response(&mut writer, 200, "text/html; charset=utf-8", WEB_UI_HTML)?;
        }
        ("GET", "/api/notebook") => {
            let nb = notebook.lock().unwrap();
            send_response(&mut writer, 200, "application/json", &nb.to_json())?;
        }
        ("POST", "/api/eval") => {
            let expr = extract_json_field(&body, "input").unwrap_or_default();
            let ctx = Context::standard();
            let steps_result = repl::dispatch_steps(&expr, ctx);
            let (output, steps_json) = match steps_result {
                Ok(steps) => {
                    let last = steps.last().map(|s| s.as_str()).unwrap_or("").to_string();
                    let steps_arr: Vec<String> = steps
                        .iter()
                        .map(|s| json_escape(s))
                        .collect();
                    (last, format!("[{}]", steps_arr.join(",")))
                }
                Err(e) => {
                    let msg = format!("Error: {}", e);
                    (msg.clone(), format!("[{}]", json_escape(&msg)))
                }
            };
            let json = format!(
                "{{\"input\": {}, \"output\": {}, \"steps\": {}}}",
                json_escape(&expr),
                json_escape(&output),
                steps_json
            );
            send_response(&mut writer, 200, "application/json", &json)?;
        }
        ("POST", "/api/notebook") => {
            let mut nb = notebook.lock().unwrap();
            *nb = crate::notebook::parse_notebook_json(&body)?;
            if let Some(ref fp) = file_path {
                let _ = nb.save(fp);
            }
            send_response(&mut writer, 200, "application/json", "{\"status\":\"saved\"}")?;
        }
        ("POST", "/api/save") => {
            let nb = notebook.lock().unwrap();
            if let Some(ref fp) = file_path {
                match nb.save(fp) {
                    Ok(_) => send_response(
                        &mut writer,
                        200,
                        "application/json",
                        "{\"status\":\"saved\"}",
                    )?,
                    Err(e) => send_response(
                        &mut writer,
                        500,
                        "application/json",
                        &format!("{{\"error\": {}}}", json_escape(&e.to_string())),
                    )?,
                }
            } else {
                send_response(
                    &mut writer,
                    200,
                    "application/json",
                    &nb.to_json(),
                )?;
            }
        }
        _ => {
            send_response(&mut writer, 404, "text/plain", "Not Found")?;
        }
    }

    Ok(())
}

fn send_response(
    writer: &mut TcpStream,
    status: u16,
    content_type: &str,
    body: &str,
) -> Result<()> {
    let status_text = match status {
        200 => "OK",
        400 => "Bad Request",
        404 => "Not Found",
        500 => "Internal Server Error",
        _ => "OK",
    };
    let response = format!(
        "HTTP/1.1 {} {}\r\nContent-Type: {}\r\nContent-Length: {}\r\nAccess-Control-Allow-Origin: *\r\nAccess-Control-Allow-Methods: GET, POST, OPTIONS\r\nAccess-Control-Allow-Headers: Content-Type\r\nConnection: close\r\n\r\n{}",
        status,
        status_text,
        content_type,
        body.len(),
        body
    );
    writer
        .write_all(response.as_bytes())
        .map_err(|e| MathError::InvalidArgument(format!("io error: {}", e)))?;
    Ok(())
}

fn extract_json_field(json: &str, field: &str) -> Option<String> {
    let pattern = format!("\"{}\"", field);
    let pos = json.find(&pattern)?;
    let rest = &json[pos + pattern.len()..];
    let colon = rest.find(':')?;
    let rest = &rest[colon + 1..];
    let rest = rest.trim_start();
    if !rest.starts_with('"') {
        return None;
    }
    let rest = &rest[1..];
    let mut result = String::new();
    let mut escaped = false;
    for ch in rest.chars() {
        if escaped {
            match ch {
                'n' => result.push('\n'),
                't' => result.push('\t'),
                '"' => result.push('"'),
                '\\' => result.push('\\'),
                c => result.push(c),
            }
            escaped = false;
        } else if ch == '\\' {
            escaped = true;
        } else if ch == '"' {
            break;
        } else {
            result.push(ch);
        }
    }
    Some(result)
}

fn json_escape(s: &str) -> String {
    let mut result = String::from("\"");
    for ch in s.chars() {
        match ch {
            '"' => result.push_str("\\\""),
            '\\' => result.push_str("\\\\"),
            '\n' => result.push_str("\\n"),
            '\t' => result.push_str("\\t"),
            '\r' => result.push_str("\\r"),
            c => result.push(c),
        }
    }
    result.push('"');
    result
}

const WEB_UI_HTML: &str = include_str!("webui.html");
