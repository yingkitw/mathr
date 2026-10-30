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

// =========================================================================
// Minimal JSON layer
// =========================================================================

/// A JSON value. Hand-rolled to keep mathr dependency-free (same convention
/// as [`crate::serialize`], which also parses JSON without `serde`).
#[derive(Debug, Clone, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Array(Vec<Json>),
    /// Object members, in insertion order.
    Object(Vec<(String, Json)>),
}

impl Json {
    /// Build a JSON object from an iterator of members.
    pub fn obj<I: IntoIterator<Item = (String, Json)>>(members: I) -> Json {
        Json::Object(members.into_iter().collect())
    }

    /// Build a JSON string value.
    pub fn str(s: impl Into<String>) -> Json {
        Json::Str(s.into())
    }

    /// Look up an object member by key.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    /// Borrow the value as a string.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    /// Borrow the value as a number.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(n) => Some(*n),
            _ => None,
        }
    }

    /// Borrow the value as a boolean.
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Json::Bool(b) => Some(*b),
            _ => None,
        }
    }

    /// Borrow the value as an array.
    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Array(items) => Some(items),
            _ => None,
        }
    }

    /// Serialize to a compact JSON string.
    ///
    /// Named `to_json_string` rather than `to_string` so it does not shadow
    /// the `Display` implementation.
    pub fn to_json_string(&self) -> String {
        let mut out = String::new();
        self.write(&mut out);
        out
    }

    fn write(&self, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(true) => out.push_str("true"),
            Json::Bool(false) => out.push_str("false"),
            Json::Num(n) => {
                if n.is_finite() {
                    // Integral values print without a fractional part so the
                    // output matches what a JSON reader expects for integers.
                    if *n == n.trunc() && n.abs() < 1e15 {
                        out.push_str(&format!("{}", *n as i64));
                    } else {
                        out.push_str(&format!("{}", n));
                    }
                } else {
                    // JSON has no NaN/Infinity; null is the conventional stand-in.
                    out.push_str("null");
                }
            }
            Json::Str(s) => write_json_string(s, out),
            Json::Array(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    item.write(out);
                }
                out.push(']');
            }
            Json::Object(members) => {
                out.push('{');
                for (i, (k, v)) in members.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    write_json_string(k, out);
                    out.push(':');
                    v.write(out);
                }
                out.push('}');
            }
        }
    }
}

impl std::fmt::Display for Json {
    /// Renders the value as compact JSON.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = String::new();
        self.write(&mut out);
        f.write_str(&out)
    }
}

/// Write a JSON string literal, escaping the characters JSON requires.
fn write_json_string(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Parse a JSON document.
pub fn parse_json(src: &str) -> Result<Json> {
    let mut p = JsonParser {
        bytes: src.as_bytes(),
        pos: 0,
    };
    let v = p.value()?;
    p.skip_ws();
    if p.pos != p.bytes.len() {
        return Err(MathError::Parse(format!(
            "unexpected trailing JSON at byte {}",
            p.pos
        )));
    }
    Ok(v)
}

struct JsonParser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> JsonParser<'a> {
    fn skip_ws(&mut self) {
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn value(&mut self) -> Result<Json> {
        self.skip_ws();
        match self.bytes.get(self.pos) {
            None => Err(MathError::Parse("unexpected end of JSON".into())),
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => self.string().map(Json::Str),
            Some(b't') => self.literal("true").map(|_| Json::Bool(true)),
            Some(b'f') => self.literal("false").map(|_| Json::Bool(false)),
            Some(b'n') => self.literal("null").map(|_| Json::Null),
            Some(c) if *c == b'-' || c.is_ascii_digit() => self.number(),
            Some(c) => Err(MathError::Parse(format!(
                "unexpected JSON char '{}' at byte {}",
                *c as char, self.pos
            ))),
        }
    }

    fn literal(&mut self, word: &str) -> Result<()> {
        if self.bytes[self.pos..].starts_with(word.as_bytes()) {
            self.pos += word.len();
            Ok(())
        } else {
            Err(MathError::Parse(format!(
                "invalid JSON literal at byte {}",
                self.pos
            )))
        }
    }

    fn number(&mut self) -> Result<Json> {
        let start = self.pos;
        if self.bytes.get(self.pos) == Some(&b'-') {
            self.pos += 1;
        }
        while self.pos < self.bytes.len()
            && (self.bytes[self.pos].is_ascii_digit() || self.bytes[self.pos] == b'.')
        {
            self.pos += 1;
        }
        // Exponent part.
        if self.pos < self.bytes.len() && (self.bytes[self.pos] | 0x20) == b'e' {
            self.pos += 1;
            if self.pos < self.bytes.len()
                && (self.bytes[self.pos] == b'+' || self.bytes[self.pos] == b'-')
            {
                self.pos += 1;
            }
            while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }
        let text = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| MathError::Parse("invalid UTF-8 in JSON number".into()))?;
        text.parse::<f64>()
            .map(Json::Num)
            .map_err(|_| MathError::Parse(format!("invalid JSON number '{}'", text)))
    }

    fn string(&mut self) -> Result<String> {
        // Skip the opening quote.
        self.pos += 1;
        let mut out = String::new();
        loop {
            let Some(&c) = self.bytes.get(self.pos) else {
                return Err(MathError::Parse("unterminated JSON string".into()));
            };
            self.pos += 1;
            match c {
                b'"' => return Ok(out),
                b'\\' => {
                    let Some(&esc) = self.bytes.get(self.pos) else {
                        return Err(MathError::Parse("unterminated JSON escape".into()));
                    };
                    self.pos += 1;
                    match esc {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => {
                            let cp = self.hex4()?;
                            // Combine a surrogate pair when present.
                            if (0xD800..0xDC00).contains(&cp)
                                && self.bytes.get(self.pos) == Some(&b'\\')
                                && self.bytes.get(self.pos + 1) == Some(&b'u')
                            {
                                self.pos += 2;
                                let low = self.hex4()?;
                                if (0xDC00..0xE000).contains(&low) {
                                    let combined =
                                        0x10000 + ((cp - 0xD800) << 10) + (low - 0xDC00);
                                    out.push(
                                        char::from_u32(combined).unwrap_or('\u{FFFD}'),
                                    );
                                } else {
                                    out.push('\u{FFFD}');
                                    out.push(char::from_u32(low).unwrap_or('\u{FFFD}'));
                                }
                            } else {
                                out.push(char::from_u32(cp).unwrap_or('\u{FFFD}'));
                            }
                        }
                        other => {
                            return Err(MathError::Parse(format!(
                                "invalid JSON escape '\\{}'",
                                other as char
                            )))
                        }
                    }
                }
                _ => {
                    // Collect the full UTF-8 sequence.
                    let len = utf8_len(c);
                    let end = (self.pos - 1 + len).min(self.bytes.len());
                    let chunk = std::str::from_utf8(&self.bytes[self.pos - 1..end])
                        .map_err(|_| MathError::Parse("invalid UTF-8 in JSON string".into()))?;
                    out.push_str(chunk);
                    self.pos = end;
                }
            }
        }
    }

    fn hex4(&mut self) -> Result<u32> {
        if self.pos + 4 > self.bytes.len() {
            return Err(MathError::Parse("truncated JSON \\u escape".into()));
        }
        let text = std::str::from_utf8(&self.bytes[self.pos..self.pos + 4])
            .map_err(|_| MathError::Parse("invalid JSON \\u escape".into()))?;
        let cp = u32::from_str_radix(text, 16)
            .map_err(|_| MathError::Parse(format!("invalid JSON \\u escape '{}'", text)))?;
        self.pos += 4;
        Ok(cp)
    }

    fn array(&mut self) -> Result<Json> {
        self.pos += 1; // [
        let mut items = Vec::new();
        self.skip_ws();
        if self.bytes.get(self.pos) == Some(&b']') {
            self.pos += 1;
            return Ok(Json::Array(items));
        }
        loop {
            items.push(self.value()?);
            self.skip_ws();
            match self.bytes.get(self.pos) {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Json::Array(items));
                }
                _ => {
                    return Err(MathError::Parse(format!(
                        "expected ',' or ']' in JSON array at byte {}",
                        self.pos
                    )))
                }
            }
        }
    }

    fn object(&mut self) -> Result<Json> {
        self.pos += 1; // {
        let mut members = Vec::new();
        self.skip_ws();
        if self.bytes.get(self.pos) == Some(&b'}') {
            self.pos += 1;
            return Ok(Json::Object(members));
        }
        loop {
            self.skip_ws();
            if self.bytes.get(self.pos) != Some(&b'"') {
                return Err(MathError::Parse(format!(
                    "expected string key in JSON object at byte {}",
                    self.pos
                )));
            }
            let key = self.string()?;
            self.skip_ws();
            if self.bytes.get(self.pos) != Some(&b':') {
                return Err(MathError::Parse(format!(
                    "expected ':' in JSON object at byte {}",
                    self.pos
                )));
            }
            self.pos += 1;
            let val = self.value()?;
            members.push((key, val));
            self.skip_ws();
            match self.bytes.get(self.pos) {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Json::Object(members));
                }
                _ => {
                    return Err(MathError::Parse(format!(
                        "expected ',' or '}}' in JSON object at byte {}",
                        self.pos
                    )))
                }
            }
        }
    }
}

/// Number of bytes in a UTF-8 sequence given its first byte.
fn utf8_len(first: u8) -> usize {
    match first {
        0x00..=0x7F => 1,
        0xC0..=0xDF => 2,
        0xE0..=0xEF => 3,
        _ => 4,
    }
}

// =========================================================================
// Tool catalogue
// =========================================================================

/// A tool definition: name, description and JSON Schema for its arguments.
pub struct Tool {
    pub name: &'static str,
    pub description: &'static str,
    pub input_schema: fn() -> Json,
}

const EXPR_DESC: &str = "Math expression. Plain syntax (x^2 + 3*x) or LaTeX (\\frac{1}{2}, \\sin(x)).";
const VAR_DESC: &str =
    "Variable to operate on (default: auto-detected, usually x).";

/// The tools exposed by the server.
///
/// The set mirrors the `unsolve` TypeScript port's MCP surface, adapted to
/// mathr's command vocabulary.
pub fn tools() -> Vec<Tool> {
    vec![
        Tool {
            name: "evaluate",
            description: "Evaluate a math expression. Exact rational arithmetic for fractions (1/2 + 3/4 -> 5/4); numeric results otherwise (sin(pi/4) -> 0.7071...). Supports LaTeX input, constants pi/e/tau, and functions sin/cos/tan/exp/ln/log/sqrt/gamma and more.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variables".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("object")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Optional variable substitutions, e.g. {\"x\": 2}",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "differentiate",
            description: "Symbolic differentiation with step-by-step working (sum, product, quotient, chain, power rules). Returns the derivative in plain syntax, e.g. diff(sin(x^2)) -> cos(x^2)*2*x.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "integrate",
            description: "Integration: a symbolic antiderivative (+ C) by default, or a numeric definite integral when both 'from' and 'to' are given. Bounds may be numbers or expressions such as \"pi\" or \"pi/2\".",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Lower bound of a definite integral (needs 'to' as well)",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Upper bound of a definite integral (needs 'from' as well)",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "solve",
            description: "Solve an equation for a variable. Exact roots for linear/quadratic polynomials (rational arithmetic), rational roots + numeric roots for higher degrees, and robust numeric root-finding for non-polynomial equations.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "equation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Equation like 'x^2 - 5*x + 6 = 0' or an expression (assumed = 0)",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("equation")])),
                ])
            },
        },
        Tool {
            name: "nsolve",
            description: "Find ALL real roots of an equation by scanning an interval for sign changes and bisecting each one. Unlike 'solve', which returns a single root near an initial guess, this returns every root: nsolve(sin(x)) over [-10, 10] gives 0, +/-pi, +/-2pi, +/-3pi. Poles are rejected (a sign flip across a pole is not a root).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Lower bound of the scan (default -30)"),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Upper bound of the scan (default 30)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "simplify",
            description: "Simplify an expression with constant folding and algebraic identities (collecting like terms, 0 + x -> x, x/1 -> x, exact fraction arithmetic).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "taylor",
            description: "Symbolic Taylor series expansion of a function around a point with N terms, e.g. exp(x) around 0 -> 1 + x + 0.5*x^2 + ...",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "at".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Expansion point (default 0)"),
                                    ),
                                ]),
                            ),
                            (
                                "terms".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    (
                                        "description".into(),
                                        Json::str("Number of terms (default 5)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "limit",
            description: "Limit of an expression as a variable approaches a point. Handles direct substitution, 0/0 and infinity/infinity via recursive L'Hopital's rule, numeric probing, and limits at +/- infinity. Example: sin(x)/x at 0 -> 1.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "at".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Point to approach (default 0)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "poly-roots",
            description: "All roots of a polynomial with numeric coefficients - real AND complex - via the Durand-Kerner iteration. Coefficients in descending order: [1, -5, 6] solves x^2 - 5x + 6 = 0 -> roots 2, 3.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "coefficients".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("number"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Polynomial coefficients, leading (highest degree) first",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("coefficients")]),
                    ),
                ])
            },
        },
        Tool {
            name: "number-theory",
            description: "Integer number theory: gcd, lcm, factor (prime factorization), is-prime (deterministic Miller-Rabin), fib (fast doubling, exact BigInt), binom (exact binomial coefficient).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "operation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "enum".into(),
                                        Json::Array(vec![
                                            Json::str("gcd"),
                                            Json::str("lcm"),
                                            Json::str("factor"),
                                            Json::str("is-prime"),
                                            Json::str("fib"),
                                            Json::str("binom"),
                                        ]),
                                    ),
                                ]),
                            ),
                            (
                                "a".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    ("description".into(), Json::str("First integer")),
                                ]),
                            ),
                            (
                                "b".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    (
                                        "description".into(),
                                        Json::str("Second integer (gcd/lcm/binom)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("operation"), Json::str("a")]),
                    ),
                ])
            },
        },
        Tool {
            name: "stats",
            description: "Descriptive statistics over a numeric dataset: n, sum, mean, median, min, max, range, sample variance, and sample standard deviation.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "data".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("number"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str("Numeric dataset"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("data")])),
                ])
            },
        },
        Tool {
            name: "matrix",
            description: "Matrix operations on a square/rectangular numeric matrix: det, inv, transpose, linsolve (solve A*x = b given the augmented n x (n+1) matrix [A | b]), trace, eigen (eigenvalues + eigenvectors of a real symmetric matrix).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "operation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "enum".into(),
                                        Json::Array(vec![
                                            Json::str("det"),
                                            Json::str("inv"),
                                            Json::str("transpose"),
                                            Json::str("linsolve"),
                                            Json::str("trace"),
                                            Json::str("eigen"),
                                        ]),
                                    ),
                                ]),
                            ),
                            (
                                "matrix".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("array"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Matrix as rows of numbers, e.g. [[4, 7], [2, 6]]; for linsolve, the augmented matrix [[2, 1, 5], [1, 3, 10]]",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("operation"), Json::str("matrix")]),
                    ),
                ])
            },
        },
        Tool {
            name: "minimize",
            description: "Minimize a single-variable function over a bounded interval via golden-section search (assumes the function is unimodal on the interval). Returns the argmin x and the minimum value f(x).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Lower bound of the search interval (a number, or an expression such as \"pi\")",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Upper bound of the search interval (a number, or an expression such as \"2*pi\")",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![
                            Json::str("expression"),
                            Json::str("from"),
                            Json::str("to"),
                        ]),
                    ),
                ])
            },
        },
        Tool {
            name: "apart",
            description: "Partial fraction decomposition of a rational function whose denominator factors into linear factors (e.g. 1/(x^2-1) -> 1/2/(x-1) - 1/2/(x+1)). Handles improper fractions via polynomial long division.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "expand",
            description: "Expand a polynomial expression into a collected sum of monomials (multivariate), e.g. (x + y)^2 -> x^2 + 2*x*y + y^2.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
    ]
}

// =========================================================================
// Tool execution
// =========================================================================

/// Run a tool by name with the given arguments, returning the text result.
///
/// Errors are returned as `Err` and reported by the caller as an MCP tool
/// execution error (`isError: true`).
pub fn call_tool(name: &str, args: &Json) -> Result<String> {
    match name {
        "evaluate" => {
            let expr = require_str(args, "expression")?;
            let mut ctx = crate::eval::Context::standard();
            if let Some(Json::Object(members)) = args.get("variables") {
                for (k, v) in members {
                    let Some(n) = v.as_f64() else {
                        return Err(MathError::Eval(format!(
                            "variable '{}' must be a number",
                            k
                        )));
                    };
                    ctx.set(k, n);
                }
            }
            let e = crate::parser::Parser::parse(&expr)?;
            // Prefer an exact rational result when the expression is rational,
            // matching the REPL's plain-evaluation path.
            if let Some(r) = crate::rational::eval_rational(&e) {
                return Ok(r.to_string());
            }
            let v = crate::eval::eval(&e, &ctx)?;
            Ok(crate::repl::format_value(v))
        }
        "differentiate" => {
            let expr = require_str(args, "expression")?;
            let cmd = format!("diff {}", expr);
            run_with_var(cmd, args)
        }
        "integrate" => {
            let expr = require_str(args, "expression")?;
            // With both bounds, compute the numeric definite integral.
            if args.get("from").is_some() && args.get("to").is_some() {
                let from = bound_text(args, "from")?;
                let to = bound_text(args, "to")?;
                let cmd = format!("int {} from {} to {}", expr, from, to);
                return run_with_var(cmd, args);
            }
            // Otherwise the symbolic antiderivative.
            let cmd = format!("integrate {}", expr);
            run_with_var(cmd, args)
        }
        "solve" => {
            let expr = require_str(args, "equation")?;
            let cmd = format!("solve {}", expr);
            run_with_var(cmd, args)
        }
        "nsolve" => {
            // All real roots by scanning for sign changes.
            let expr = require_str(args, "expression")?;
            let mut cmd = format!("nsolve {}", expr);
            if let Some(from) = args.get("from").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(from)));
            }
            if let Some(to) = args.get("to").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(to)));
            }
            run_with_var(cmd, args)
        }
        "simplify" => {
            let expr = require_str(args, "expression")?;
            run_command(&format!("simplify {}", expr))
        }
        "expand" => {
            let expr = require_str(args, "expression")?;
            run_command(&format!("expand {}", expr))
        }
        "apart" => {
            let expr = require_str(args, "expression")?;
            run_with_var(format!("apart {}", expr), args)
        }
        "taylor" => {
            let expr = require_str(args, "expression")?;
            let at = args.get("at").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let terms = args
                .get("terms")
                .and_then(|v| v.as_f64())
                .unwrap_or(5.0)
                .clamp(1.0, 64.0) as usize;
            let mut cmd = format!("taylor {}", expr);
            if at != 0.0 {
                cmd.push_str(&format!(" {}", crate::repl::format_value(at)));
            }
            cmd.push_str(&format!(" {}", terms));
            run_with_var(cmd, args)
        }
        "limit" => {
            let expr = require_str(args, "expression")?;
            let at = args.get("at").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let mut cmd = format!("limit {}", expr);
            cmd.push_str(&format!(" {}", crate::repl::format_value(at)));
            run_with_var(cmd, args)
        }
        "poly-roots" => {
            let coeffs = args
                .get("coefficients")
                .and_then(|v| v.as_array())
                .ok_or_else(|| {
                    MathError::Eval("poly-roots needs a 'coefficients' array".into())
                })?;
            let mut parts = Vec::with_capacity(coeffs.len());
            for c in coeffs {
                let n = c.as_f64().ok_or_else(|| {
                    MathError::Eval("poly-roots coefficients must be numbers".into())
                })?;
                parts.push(crate::repl::format_value(n));
            }
            run_command(&format!("poly-roots {}", parts.join(" ")))
        }
        "number-theory" => {
            let op = require_str(args, "operation")?;
            let a = require_num(args, "a")?;
            let mut cmd = format!("{} {}", op, crate::repl::format_value(a));
            if let Some(b) = args.get("b").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(b)));
            }
            run_command(&cmd)
        }
        "stats" => {
            let data = args
                .get("data")
                .and_then(|v| v.as_array())
                .ok_or_else(|| MathError::Eval("stats needs a 'data' array".into()))?;
            let mut parts = Vec::with_capacity(data.len());
            for d in data {
                let n = d
                    .as_f64()
                    .ok_or_else(|| MathError::Eval("stats data must be numbers".into()))?;
                parts.push(crate::repl::format_value(n));
            }
            run_command(&format!("stats {}", parts.join(" ")))
        }
        "matrix" => {
            let op = require_str(args, "operation")?;
            let rows = args
                .get("matrix")
                .and_then(|v| v.as_array())
                .ok_or_else(|| MathError::Eval("matrix needs a 'matrix' array".into()))?;
            // mathr's matrix syntax separates rows with `|` and entries with
            // commas or whitespace.
            let mut lines = Vec::with_capacity(rows.len());
            for row in rows {
                let cells = row.as_array().ok_or_else(|| {
                    MathError::Eval("matrix must be an array of rows".into())
                })?;
                let mut parts = Vec::with_capacity(cells.len());
                for c in cells {
                    let n = c.as_f64().ok_or_else(|| {
                        MathError::Eval("matrix entries must be numbers".into())
                    })?;
                    parts.push(crate::repl::format_value(n));
                }
                lines.push(parts.join(", "));
            }
            run_command(&format!("{} {}", op, lines.join(" | ")))
        }
        "minimize" => {
            let expr = require_str(args, "expression")?;
            // The bounds may be numbers or expressions ("pi", "pi/2").
            let from = bound_text(args, "from")?;
            let to = bound_text(args, "to")?;
            let cmd = format!("minimize {} from {} to {}", expr, from, to);
            run_with_var(cmd, args)
        }
        other => Err(MathError::Eval(format!("unknown tool: {}", other))),
    }
}

/// Read a bound argument as source text.
///
/// A number is formatted the REPL way; a string is passed through so
/// expressions like `pi` or `pi/2` work.
fn bound_text(args: &Json, key: &str) -> Result<String> {
    match args.get(key) {
        Some(Json::Num(n)) => Ok(crate::repl::format_value(*n)),
        Some(Json::Str(s)) => Ok(s.trim().to_string()),
        _ => Err(MathError::Eval(format!("missing required argument '{}'", key))),
    }
}

/// Fetch a required string argument.
fn require_str(args: &Json, key: &str) -> Result<String> {
    args.get(key)
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
        .ok_or_else(|| MathError::Eval(format!("missing required argument '{}'", key)))
}

/// Fetch a required numeric argument.
fn require_num(args: &Json, key: &str) -> Result<f64> {
    args.get(key)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| MathError::Eval(format!("missing required argument '{}'", key)))
}

/// Run a command, optionally appending `for <variable>`.
fn run_with_var(mut cmd: String, args: &Json) -> Result<String> {
    if let Some(var) = args.get("variable").and_then(|v| v.as_str()) {
        if !var.is_empty() {
            cmd.push_str(&format!(" for {}", var));
        }
    }
    run_command(&cmd)
}

/// Run a mathr command string and return its output.
fn run_command(cmd: &str) -> Result<String> {
    let ctx = crate::eval::Context::standard();
    match crate::repl::dispatch_str(cmd, ctx)? {
        Some(out) => Ok(out),
        None => Ok(String::new()),
    }
}

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
    fn apart_tool_works() {
        let res = call(
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
