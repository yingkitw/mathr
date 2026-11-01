//! A Model Context Protocol (MCP) server exposing mathr's engine as tools.
//!
//! The server speaks **JSON-RPC 2.0 over stdio**: newline-delimited JSON
//! messages on stdin/stdout, per the MCP stdio transport
//! (<https://modelcontextprotocol.io/specification/2025-06-18/basic/transports>).
//! Messages must not contain embedded newlines.
//!
//! Launch it with:
//!
//! ```text
//! mathr mcp
//! ```
//!
//! Point an MCP client at it with a server config such as:
//!
//! ```json
//! { "mcpServers": { "mathr": { "command": "mathr", "args": ["mcp"] } } }
//! ```
//!
//! # Design
//!
//! Every tool is a thin wrapper that builds a mathr command string and runs
//! it through [`crate::repl::dispatch_str`]. The REPL is the single source
//! of truth for behaviour, so the MCP surface cannot drift from the CLI.
//!
//! Cross-port pattern from the `unsolve` TypeScript port, which exposes the
//! same engine as 14 MCP tools over Streamable HTTP
//! (`../unsolve/api/mcp.ts`).
//!
//! # Protocol surface
//!
//! - `initialize` — capability negotiation; replies with the server's
//!   protocol version, tool capability, and server info.
//! - `notifications/initialized` — accepted and ignored.
//! - `tools/list` — the tool catalogue with JSON Schemas.
//! - `tools/call` — invoke a tool.
//! - `ping` — liveness check.
//!
//! Tool *execution* failures are reported as a successful JSON-RPC response
//! carrying `isError: true` (the MCP convention), while protocol-level
//! problems (unknown tool, bad arguments) use JSON-RPC error codes.

use crate::error::{MathError, Result};
use std::io::{BufRead, Write};

/// The MCP protocol version this server implements.
pub const PROTOCOL_VERSION: &str = "2025-06-18";

/// The server name reported in `initialize`.
pub const SERVER_NAME: &str = "mathr";

/// Server version reported in `initialize`.
pub const SERVER_VERSION: &str = env!("CARGO_PKG_VERSION");

mod json;
mod tools;

pub use json::{parse_json, Json};
pub use tools::{call_tool, tools, Tool};

// =========================================================================
// JSON-RPC dispatch
// =========================================================================

/// A JSON-RPC error code.
#[derive(Debug, Clone)]
pub struct RpcError {
    pub code: i64,
    pub message: String,
}

/// The outcome of handling one JSON-RPC message.
///
/// Spelled out rather than using `Result` because mathr's crate-level
/// `Result<T>` alias has a single type parameter, and this layer needs the
/// two-parameter form.
pub type RpcResult = std::result::Result<Json, RpcError>;

impl RpcError {
    fn new(code: i64, message: impl Into<String>) -> RpcError {
        RpcError {
            code,
            message: message.into(),
        }
    }
}

impl std::fmt::Display for RpcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

/// JSON-RPC error codes used by the server.
pub const METHOD_NOT_FOUND: i64 = -32601;
pub const INVALID_PARAMS: i64 = -32602;
pub const INTERNAL_ERROR: i64 = -32603;

/// Handle one incoming JSON-RPC message and return the response, if any.
///
/// Returns `None` for notifications (messages without an `id`), which must
/// not be answered.
pub fn handle_message(msg: &Json) -> Option<RpcResult> {
    let method = match msg.get("method").and_then(|m| m.as_str()) {
        Some(m) => m,
        None => {
            // Not a request we can answer (e.g. a client response).
            return None;
        }
    };
    let id = msg.get("id").cloned();
    let params = msg.get("params").cloned().unwrap_or(Json::Null);

    let outcome: RpcResult = match method {
        "initialize" => Ok(initialize_result()),
        "notifications/initialized" | "initialized" => return None,
        "ping" => Ok(Json::obj(vec![])),
        "tools/list" => Ok(tools_list_result()),
        "tools/call" => call_tool_result(&params),
        other => Err(RpcError::new(
            METHOD_NOT_FOUND,
            format!("method not found: {}", other),
        )),
    };

    let id = id?;
    let payload = match outcome {
        Ok(result) => Json::obj(vec![
            ("jsonrpc".into(), Json::str("2.0")),
            ("id".into(), id),
            ("result".into(), result),
        ]),
        Err(e) => Json::obj(vec![
            ("jsonrpc".into(), Json::str("2.0")),
            ("id".into(), id),
            (
                "error".into(),
                Json::obj(vec![
                    ("code".into(), Json::Num(e.code as f64)),
                    ("message".into(), Json::str(e.message)),
                ]),
            ),
        ]),
    };
    Some(Ok(payload))
}

/// The `initialize` result payload.
fn initialize_result() -> Json {
    Json::obj(vec![
        ("protocolVersion".into(), Json::str(PROTOCOL_VERSION)),
        (
            "capabilities".into(),
            Json::obj(vec![("tools".into(), Json::obj(vec![]))]),
        ),
        (
            "serverInfo".into(),
            Json::obj(vec![
                ("name".into(), Json::str(SERVER_NAME)),
                ("version".into(), Json::str(SERVER_VERSION)),
            ]),
        ),
        (
            "instructions".into(),
            Json::str(
                "mathr is a deterministic math engine: symbolic calculus, equation solving, \
                 exact linear algebra, number theory and statistics. No LLM is involved in \
                 the solve path. Prefer the 'evaluate' tool for plain expressions and \
                 'solve' for equations.",
            ),
        ),
    ])
}

/// The `tools/list` result payload.
fn tools_list_result() -> Json {
    let list: Vec<Json> = tools()
        .iter()
        .map(|t| {
            Json::obj(vec![
                ("name".into(), Json::str(t.name)),
                ("description".into(), Json::str(t.description)),
                ("inputSchema".into(), (t.input_schema)()),
            ])
        })
        .collect();
    Json::obj(vec![("tools".into(), Json::Array(list))])
}

/// Execute a `tools/call` request.
///
/// A tool that fails produces a *successful* JSON-RPC response whose result
/// carries `isError: true` — the MCP convention for execution failures.
/// Protocol-level problems (unknown tool, missing arguments) use JSON-RPC
/// error codes instead.
fn call_tool_result(params: &Json) -> RpcResult {
    let name = params
        .get("name")
        .and_then(|v| v.as_str())
        .ok_or_else(|| RpcError::new(INVALID_PARAMS, "missing tool name"))?;
    let known = tools().iter().any(|t| t.name == name);
    if !known {
        return Err(RpcError::new(
            INVALID_PARAMS,
            format!("unknown tool: {}", name),
        ));
    }
    let args = params.get("arguments").cloned().unwrap_or(Json::Null);
    match call_tool(name, &args) {
        Ok(text) => Ok(Json::obj(vec![
            (
                "content".into(),
                Json::Array(vec![Json::obj(vec![
                    ("type".into(), Json::str("text")),
                    ("text".into(), Json::str(text)),
                ])]),
            ),
            ("isError".into(), Json::Bool(false)),
        ])),
        Err(e) => Ok(Json::obj(vec![
            (
                "content".into(),
                Json::Array(vec![Json::obj(vec![
                    ("type".into(), Json::str("text")),
                    ("text".into(), Json::str(format!("error: {}", e))),
                ])]),
            ),
            ("isError".into(), Json::Bool(true)),
        ])),
    }
}

// =========================================================================
// stdio server loop
// =========================================================================

/// Run the MCP server over stdin/stdout until the input stream closes.
///
/// Each line of stdin is one JSON-RPC message; each response is written to
/// stdout as a single line. Nothing else may be written to stdout — the
/// MCP spec forbids it. Diagnostics go to stderr.
pub fn serve_stdio() -> Result<()> {
    let stdin = std::io::stdin();
    let mut stdout = std::io::stdout();
    let mut line = String::new();
    loop {
        line.clear();
        let n = stdin.lock().read_line(&mut line)?;
        if n == 0 {
            // Input closed: the client is shutting us down.
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let response = match parse_json(trimmed) {
            Ok(msg) => handle_message(&msg),
            Err(e) => Some(Err(RpcError::new(
                INVALID_PARAMS,
                format!("invalid JSON: {}", e),
            ))),
        };
        if let Some(outcome) = response {
            let payload = match outcome {
                Ok(v) => v,
                Err(e) => Json::obj(vec![
                    ("jsonrpc".into(), Json::str("2.0")),
                    ("id".into(), Json::Null),
                    (
                        "error".into(),
                        Json::obj(vec![
                            ("code".into(), Json::Num(e.code as f64)),
                            ("message".into(), Json::str(e.message)),
                        ]),
                    ),
                ]),
            };
            writeln!(stdout, "{}", payload.to_json_string())?;
            stdout.flush()?;
        }
    }
    Ok(())
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    /// The arguments object for a tool call.
    fn arg(v: Json) -> Json {
        v
    }

    /// Invoke a tool and return its JSON-RPC result.
    fn call(name: &str, args: Json) -> RpcResult {
        let params = Json::obj(vec![
            ("name".into(), Json::str(name)),
            ("arguments".into(), args),
        ]);
        call_tool_result(&params)
    }

    fn text_of(result: &Json) -> String {
        result
            .get("content")
            .and_then(|c| c.as_array())
            .and_then(|a| a.first())
            .and_then(|c| c.get("text"))
            .and_then(|t| t.as_str())
            .unwrap_or("")
            .to_string()
    }

    // ---- JSON layer -------------------------------------------------------

    #[test]
    fn json_parses_object() {
        let v = parse_json(r#"{"a": 1, "b": "two", "c": [1, 2, 3]}"#).unwrap();
        assert_eq!(v.get("a").and_then(|x| x.as_f64()), Some(1.0));
        assert_eq!(v.get("b").and_then(|x| x.as_str()), Some("two"));
        assert_eq!(v.get("c").and_then(|x| x.as_array()).unwrap().len(), 3);
    }

    #[test]
    fn json_parses_nested_and_escapes() {
        let v = parse_json(r#"{"a": {"b": [1, {"c": "d\"e"}]}}"#).unwrap();
        let inner = v.get("a").unwrap().get("b").unwrap().as_array().unwrap();
        let obj = &inner[1];
        assert_eq!(obj.get("c").and_then(|x| x.as_str()), Some("d\"e"));
    }

    #[test]
    fn json_parses_unicode_escape() {
        let v = parse_json(r#"{"s": "π"}"#).unwrap();
        assert_eq!(v.get("s").and_then(|x| x.as_str()), Some("π"));
    }

    #[test]
    fn json_parses_surrogate_pair() {
        let v = parse_json(r#"{"s": "😀"}"#).unwrap();
        assert_eq!(v.get("s").and_then(|x| x.as_str()), Some("😀"));
    }

    #[test]
    fn json_parses_negative_and_exponent_numbers() {
        let v = parse_json(r#"{"a": -1.5, "b": 2e3}"#).unwrap();
        assert_eq!(v.get("a").and_then(|x| x.as_f64()), Some(-1.5));
        assert_eq!(v.get("b").and_then(|x| x.as_f64()), Some(2000.0));
    }

    #[test]
    fn json_rejects_trailing_content() {
        assert!(parse_json(r#"{"a":1} extra"#).is_err());
    }

    #[test]
    fn json_rejects_unterminated_string() {
        assert!(parse_json(r#"{"a": "b"#).is_err());
    }

    #[test]
    fn json_serializes_escapes() {
        let v = Json::str("a\"b\\c\nd");
        assert_eq!(v.to_json_string(), r#""a\"b\\c\nd""#);
    }

    #[test]
    fn json_serializes_integral_numbers_without_fraction() {
        assert_eq!(Json::Num(3.0).to_json_string(), "3");
        assert_eq!(Json::Num(3.5).to_json_string(), "3.5");
    }

    #[test]
    fn json_serializes_non_finite_as_null() {
        assert_eq!(Json::Num(f64::NAN).to_json_string(), "null");
        assert_eq!(Json::Num(f64::INFINITY).to_json_string(), "null");
    }

    #[test]
    fn json_roundtrip() {
        let src = r#"{"a":[1,2,{"b":"c"}],"d":true,"e":null}"#;
        let v = parse_json(src).unwrap();
        let again = parse_json(&v.to_json_string()).unwrap();
        assert_eq!(v, again);
    }

    // ---- protocol ---------------------------------------------------------

    #[test]
    fn initialize_reports_protocol_and_tools() {
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":1,"method":"initialize"}"#).unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        let result = res.get("result").unwrap();
        assert_eq!(
            result.get("protocolVersion").and_then(|v| v.as_str()),
            Some(PROTOCOL_VERSION)
        );
        assert!(result.get("capabilities").unwrap().get("tools").is_some());
        assert_eq!(
            result.get("serverInfo").unwrap().get("name").and_then(|v| v.as_str()),
            Some("mathr")
        );
    }

    #[test]
    fn initialize_echoes_request_id() {
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":42,"method":"initialize"}"#).unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        assert_eq!(res.get("id").and_then(|v| v.as_f64()), Some(42.0));
    }

    #[test]
    fn notification_gets_no_response() {
        let msg =
            parse_json(r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#).unwrap();
        assert!(handle_message(&msg).is_none());
    }

    #[test]
    fn unknown_method_is_method_not_found() {
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":1,"method":"nope"}"#).unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        let err = res.get("error").unwrap();
        assert_eq!(err.get("code").and_then(|v| v.as_f64()), Some(METHOD_NOT_FOUND as f64));
    }

    #[test]
    fn unknown_tool_is_invalid_params() {
        let res = call("not-a-tool", Json::obj(vec![])).unwrap_err();
        assert_eq!(res.code, INVALID_PARAMS);
    }

    #[test]
    fn tools_list_includes_expected_tools() {
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#).unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        let list = res
            .get("result")
            .unwrap()
            .get("tools")
            .unwrap()
            .as_array()
            .unwrap()
            .to_vec();
        let names: Vec<&str> = list
            .iter()
            .filter_map(|t| t.get("name").and_then(|n| n.as_str()))
            .collect();
        for expected in [
            "evaluate",
            "differentiate",
            "integrate",
            "solve",
            "nsolve",
            "simplify",
            "taylor",
            "limit",
            "poly-roots",
            "number-theory",
            "stats",
            "matrix",
            "minimize",
            "apart",
        ] {
            assert!(names.contains(&expected), "missing tool {}", expected);
        }
    }

    #[test]
    fn every_tool_has_a_schema() {
        for t in tools() {
            let schema = (t.input_schema)();
            assert_eq!(
                schema.get("type").and_then(|v| v.as_str()),
                Some("object"),
                "tool {} lacks an object schema",
                t.name
            );
            assert!(schema.get("properties").is_some(), "tool {} lacks properties", t.name);
            assert!(!t.description.is_empty(), "tool {} lacks a description", t.name);
        }
    }

    #[test]
    fn tool_names_are_unique() {
        let names: Vec<&str> = tools().iter().map(|t| t.name).collect();
        let mut sorted = names.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(names.len(), sorted.len());
    }

    // ---- tool execution ---------------------------------------------------

    #[test]
    fn evaluate_tool_returns_number() {
        let res = call("evaluate", arg(Json::obj(vec![("expression".into(), Json::str("2+3"))]))).unwrap();
        assert_eq!(text_of(&res), "5");
        assert_eq!(res.get("isError").and_then(|v| v.as_bool()), Some(false));
    }

    #[test]
    fn evaluate_tool_is_exact_for_rationals() {
        let res = call(
            "evaluate",
            arg(Json::obj(vec![("expression".into(), Json::str("1/2 + 3/4"))])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "5/4");
    }

    #[test]
    fn evaluate_tool_accepts_variables() {
        let res = call(
            "evaluate",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("x^2")),
                (
                    "variables".into(),
                    Json::obj(vec![("x".into(), Json::Num(3.0))]),
                ),
            ])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "9");
    }

    #[test]
    fn differentiate_tool_works() {
        let res = call(
            "differentiate",
            arg(Json::obj(vec![("expression".into(), Json::str("x^3"))])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "3*x^2");
    }

    #[test]
    fn differentiate_tool_honours_variable() {
        let res = call(
            "differentiate",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("y*x^2")),
                ("variable".into(), Json::str("y")),
            ])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "x^2");
    }

    #[test]
    fn integrate_tool_works() {
        let res = call(
            "integrate",
            arg(Json::obj(vec![("expression".into(), Json::str("x^2"))])),
        )
        .unwrap();
        assert!(text_of(&res).contains("x^3"), "got {}", text_of(&res));
    }

    #[test]
    fn integrate_tool_with_bounds_is_numeric() {
        // With both bounds, the tool computes the definite integral.
        let res = call(
            "integrate",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("sin(x)")),
                ("from".into(), Json::Num(0.0)),
                ("to".into(), Json::Num(std::f64::consts::PI)),
            ])),
        )
        .unwrap();
        assert!(text_of(&res).contains("= 2"), "got {}", text_of(&res));
    }

    #[test]
    fn integrate_tool_accepts_expression_bounds() {
        let res = call(
            "integrate",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("sin(x)")),
                ("from".into(), Json::str("0")),
                ("to".into(), Json::str("pi/2")),
            ])),
        )
        .unwrap();
        assert!(text_of(&res).contains("= 1"), "got {}", text_of(&res));
    }

    #[test]
    fn minimize_tool_accepts_expression_bounds() {
        let res = call(
            "minimize",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("(x-3)^2")),
                ("from".into(), Json::str("0")),
                ("to".into(), Json::str("10")),
            ])),
        )
        .unwrap();
        assert!(text_of(&res).contains("x* = 3"), "got {}", text_of(&res));
    }

    #[test]
    fn solve_tool_works() {
        let res = call(
            "solve",
            arg(Json::obj(vec![("equation".into(), Json::str("x^2 - 4"))])),
        )
        .unwrap();
        assert!(text_of(&res).contains("2"), "got {}", text_of(&res));
    }

    #[test]
    fn nsolve_tool_works() {
        let res = call(
            "nsolve",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("sin(x)")),
                ("from".into(), Json::Num(-7.0)),
                ("to".into(), Json::Num(7.0)),
            ])),
        )
        .unwrap();
        let text = text_of(&res);
        // Five roots: -2pi, -pi, 0, pi, 2pi.
        assert!(text.contains("6.2831853072"), "got {}", text);
        assert!(text.contains("3.1415926536"), "got {}", text);
        assert!(text.contains("roots:"), "got {}", text);
        // Count the roots: five comma-separated numbers.
        let list = text.strip_prefix("roots: ").expect("roots list");
        assert_eq!(list.split(", ").count(), 5, "got {}", text);
    }

    #[test]
    fn nsolve_tool_rejects_poles() {
        let res = call(
            "nsolve",
            arg(Json::obj(vec![("expression".into(), Json::str("1/(x-1)"))])),
        )
        .unwrap();
        assert!(
            text_of(&res).contains("no real roots"),
            "got {}",
            text_of(&res)
        );
    }

    #[test]
    fn simplify_tool_works() {
        let res = call(
            "simplify",
            arg(Json::obj(vec![("expression".into(), Json::str("2*x + 3*x"))])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "5*x");
    }

    #[test]
    fn expand_tool_works() {
        let res = call(
            "expand",
            arg(Json::obj(vec![("expression".into(), Json::str("(x + y)^2"))])),
        )
        .unwrap();
        let text = text_of(&res);
        assert!(text.contains("x^2"), "got {}", text);
        assert!(text.contains("y^2"), "got {}", text);
    }

    #[test]
    fn msolve_tool_works() {
        let res = call(
            "msolve",
            arg(Json::obj(vec![(
                "equations".into(),
                Json::str("2x + y = 5; x + 3y = 10"),
            )])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "x = 1, y = 3");
    }

    #[test]
    fn msolve_tool_keeps_fractions_exact() {
        let res = call(
            "msolve",
            arg(Json::obj(vec![(
                "equations".into(),
                Json::str("(1/2)x + y = 1; x + (1/3)y = 1"),
            )])),
        )
        .unwrap();
        assert!(text_of(&res).contains("4/5"), "got {}", text_of(&res));
    }

    #[test]
    fn expand_in_tool_works() {
        let res = call(
            "expand-in",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("(x + sin(y))^2")),
                ("variable".into(), Json::str("x")),
            ])),
        )
        .unwrap();
        assert!(text_of(&res).contains("sin(y)"), "got {}", text_of(&res));
        assert!(text_of(&res).contains("x^2"), "got {}", text_of(&res));
    }

    #[test]
    fn apart_tool_works() {        let res = call(
            "apart",
            arg(Json::obj(vec![("expression".into(), Json::str("1/(x^2-1)"))])),
        )
        .unwrap();
        assert!(!text_of(&res).is_empty());
    }

    #[test]
    fn taylor_tool_works() {
        let res = call(
            "taylor",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("exp(x)")),
                ("terms".into(), Json::Num(3.0)),
            ])),
        )
        .unwrap();
        let text = text_of(&res);
        assert!(text.contains("x"), "got {}", text);
    }

    #[test]
    fn limit_tool_works() {
        let res = call(
            "limit",
            arg(Json::obj(vec![("expression".into(), Json::str("sin(x)/x"))])),
        )
        .unwrap();
        assert!(text_of(&res).ends_with("= 1"), "got {}", text_of(&res));
    }

    #[test]
    fn poly_roots_tool_works() {
        let res = call(
            "poly-roots",
            arg(Json::obj(vec![(
                "coefficients".into(),
                Json::Array(vec![Json::Num(1.0), Json::Num(-5.0), Json::Num(6.0)]),
            )])),
        )
        .unwrap();
        let text = text_of(&res);
        assert!(text.contains('2') && text.contains('3'), "got {}", text);
    }

    #[test]
    fn number_theory_tool_works() {
        let res = call(
            "number-theory",
            arg(Json::obj(vec![
                ("operation".into(), Json::str("gcd")),
                ("a".into(), Json::Num(48.0)),
                ("b".into(), Json::Num(36.0)),
            ])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "12");
    }

    #[test]
    fn stats_tool_works() {
        let res = call(
            "stats",
            arg(Json::obj(vec![(
                "data".into(),
                Json::Array(vec![
                    Json::Num(1.0),
                    Json::Num(2.0),
                    Json::Num(3.0),
                    Json::Num(4.0),
                ]),
            )])),
        )
        .unwrap();
        assert!(text_of(&res).contains("mean"), "got {}", text_of(&res));
    }

    #[test]
    fn matrix_det_tool_works() {
        let res = call(
            "matrix",
            arg(Json::obj(vec![
                ("operation".into(), Json::str("det")),
                (
                    "matrix".into(),
                    Json::Array(vec![
                        Json::Array(vec![Json::Num(4.0), Json::Num(7.0)]),
                        Json::Array(vec![Json::Num(2.0), Json::Num(6.0)]),
                    ]),
                ),
            ])),
        )
        .unwrap();
        assert_eq!(text_of(&res), "det = 10");
    }

    #[test]
    fn minimize_tool_works() {
        let res = call(
            "minimize",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("(x-3)^2")),
                ("from".into(), Json::Num(0.0)),
                ("to".into(), Json::Num(10.0)),
            ])),
        )
        .unwrap();
        assert!(!text_of(&res).is_empty());
    }

    #[test]
    fn missing_argument_is_tool_error() {
        let res = call("evaluate", arg(Json::obj(vec![]))).unwrap();
        assert_eq!(res.get("isError").and_then(|v| v.as_bool()), Some(true));
        assert!(text_of(&res).contains("error"), "got {}", text_of(&res));
    }

    #[test]
    fn bad_expression_is_tool_error_not_protocol_error() {
        let res = call(
            "evaluate",
            arg(Json::obj(vec![("expression".into(), Json::str("@@@"))])),
        )
        .unwrap();
        assert_eq!(res.get("isError").and_then(|v| v.as_bool()), Some(true));
    }

    #[test]
    fn error_response_carries_id() {
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":7,"method":"tools/call","params":{"name":"nope"}}"#)
            .unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        assert_eq!(res.get("id").and_then(|v| v.as_f64()), Some(7.0));
        assert!(res.get("error").is_some());
    }

    #[test]
    fn tool_call_response_is_success_with_content() {
        let msg = parse_json(
            r#"{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"simplify","arguments":{"expression":"x + 0"}}}"#,
        )
        .unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        let result = res.get("result").unwrap();
        assert_eq!(result.get("isError").and_then(|v| v.as_bool()), Some(false));
        assert_eq!(text_of(result), "x");
    }

    #[test]
    fn response_is_single_line() {
        // The stdio transport requires one message per line.
        let msg = parse_json(r#"{"jsonrpc":"2.0","id":1,"method":"tools/list"}"#).unwrap();
        let res = handle_message(&msg).unwrap().unwrap();
        assert!(!res.to_json_string().contains('\n'));
    }

    #[test]
    fn responses_do_not_contain_raw_newlines() {
        // Even a multi-line tool result must be escaped inside the JSON text.
        let res = call(
            "taylor",
            arg(Json::obj(vec![
                ("expression".into(), Json::str("exp(x)")),
                ("terms".into(), Json::Num(3.0)),
            ])),
        )
        .unwrap();
        assert!(!res.to_json_string().contains('\n'));
    }
}
