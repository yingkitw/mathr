//! Expression serialization for [`Expr`] in three interchangeable formats:
//!
//! - **S-expressions** (Lisp-like prefix notation): [`to_sexpr`] / [`from_sexpr`]
//! - **JSON** (nested objects): [`to_json`] / [`from_json`]
//! - **RPN** (postfix, stack-based): [`to_rpn`] / [`from_rpn`]
//!
//! All three round-trip: `from_*(to_*(e)) == e`. No external dependencies —
//! the JSON reader is a small hand-rolled recursive-descent parser, matching
//! the rest of the crate's "no serde" convention.
//!
//! ## Example
//! ```
//! use mathr::expr::Expr;
//! use mathr::serialize::{to_sexpr, from_sexpr};
//! let e = Expr::add(Expr::mul(Expr::num(2.0), Expr::var("x")), Expr::num(1.0));
//! let s = to_sexpr(&e);
//! assert_eq!(s, "(add (mul (num 2) (var x)) (num 1))");
//! let e2 = from_sexpr(&s).unwrap();
//! assert!(e.equals(&e2));
//! ```

use crate::error::{MathError, Result};
use crate::expr::Expr;

// =========================================================================
// Shared number formatting
// =========================================================================

/// Format an `f64` the way the rest of mathr does: integers as `i64`,
/// non-finite as `NaN` / `inf` / `-inf`. Used by all three formats so the
/// textual forms stay consistent with `Expr::Display`.
fn fmt_num(n: f64) -> String {
    if n.is_nan() {
        "NaN".to_string()
    } else if n.is_infinite() {
        if n > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        }
    } else if n.fract() == 0.0 && n.abs() < 1e15 {
        format!("{}", n as i64)
    } else {
        format!("{}", n)
    }
}

/// Parse a number string produced by [`fmt_num`] (or any valid f64 literal).
fn parse_num(s: &str) -> Result<f64> {
    match s {
        "NaN" => Ok(f64::NAN),
        "inf" => Ok(f64::INFINITY),
        "-inf" => Ok(f64::NEG_INFINITY),
        _ => s
            .parse::<f64>()
            .map_err(|_| MathError::Parse(format!("invalid number: {}", s))),
    }
}

// =========================================================================
// S-expressions
// =========================================================================

/// Convert an [`Expr`] to a Lisp-like S-expression.
///
/// Form: `(num 2)`, `(var x)`, `(neg e)`, `(add a b)`, `(sub a b)`,
/// `(mul a b)`, `(div a b)`, `(pow a b)`, `(func sin x)`.
pub fn to_sexpr(e: &Expr) -> String {
    let mut s = String::new();
    write_sexpr(e, &mut s);
    s
}

fn write_sexpr(e: &Expr, out: &mut String) {
    match e {
        Expr::Num(n) => out.push_str(&format!("(num {})", fmt_num(*n))),
        Expr::Var(v) => out.push_str(&format!("(var {})", v)),
        Expr::Neg(a) => {
            out.push_str("(neg ");
            write_sexpr(a, out);
            out.push(')');
        }
        Expr::Add(a, b) => {
            out.push_str("(add ");
            write_sexpr(a, out);
            out.push(' ');
            write_sexpr(b, out);
            out.push(')');
        }
        Expr::Sub(a, b) => {
            out.push_str("(sub ");
            write_sexpr(a, out);
            out.push(' ');
            write_sexpr(b, out);
            out.push(')');
        }
        Expr::Mul(a, b) => {
            out.push_str("(mul ");
            write_sexpr(a, out);
            out.push(' ');
            write_sexpr(b, out);
            out.push(')');
        }
        Expr::Div(a, b) => {
            out.push_str("(div ");
            write_sexpr(a, out);
            out.push(' ');
            write_sexpr(b, out);
            out.push(')');
        }
        Expr::Pow(a, b) => {
            out.push_str("(pow ");
            write_sexpr(a, out);
            out.push(' ');
            write_sexpr(b, out);
            out.push(')');
        }
        Expr::Func(name, args) => {
            out.push_str("(func ");
            out.push_str(name);
            for a in args {
                out.push(' ');
                write_sexpr(a, out);
            }
            out.push(')');
        }
    }
}

/// Parse an S-expression produced by [`to_sexpr`] back into an [`Expr`].
pub fn from_sexpr(input: &str) -> Result<Expr> {
    let tokens = sexpr_tokenize(input)?;
    let mut p = SexprParser { tokens, pos: 0 };
    let e = p.parse_expr()?;
    if p.pos < p.tokens.len() {
        return Err(MathError::Parse(format!(
            "unexpected trailing S-expr token at {}",
            p.pos
        )));
    }
    Ok(e)
}

#[derive(Debug, Clone)]
enum SexprTok {
    Open,        // (
    Close,       // )
    Atom(String), // bare token
}

fn sexpr_tokenize(input: &str) -> Result<Vec<SexprTok>> {
    let mut toks = Vec::new();
    let bytes = input.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let c = bytes[i];
        if c.is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if c == b'(' {
            toks.push(SexprTok::Open);
            i += 1;
            continue;
        }
        if c == b')' {
            toks.push(SexprTok::Close);
            i += 1;
            continue;
        }
        // Atom: read until whitespace or paren
        let start = i;
        while i < bytes.len() && !bytes[i].is_ascii_whitespace() && bytes[i] != b'(' && bytes[i] != b')' {
            i += 1;
        }
        let atom = &input[start..i];
        if atom.is_empty() {
            return Err(MathError::Parse("empty S-expr atom".into()));
        }
        toks.push(SexprTok::Atom(atom.to_string()));
    }
    Ok(toks)
}

struct SexprParser {
    tokens: Vec<SexprTok>,
    pos: usize,
}

impl SexprParser {
    fn parse_expr(&mut self) -> Result<Expr> {
        match self.tokens.get(self.pos) {
            Some(SexprTok::Open) => {
                self.pos += 1;
                let head = match self.tokens.get(self.pos) {
                    Some(SexprTok::Atom(h)) => h.clone(),
                    _ => return Err(MathError::Parse("expected head atom after '('".into())),
                };
                self.pos += 1;
                let mut args: Vec<Expr> = Vec::new();
                while !matches!(self.tokens.get(self.pos), Some(SexprTok::Close) | None) {
                    args.push(self.parse_expr()?);
                }
                match self.tokens.get(self.pos) {
                    Some(SexprTok::Close) => self.pos += 1,
                    _ => return Err(MathError::Parse("missing ')' in S-expr".into())),
                }
                build_node(&head, args)
            }
            Some(SexprTok::Atom(a)) => {
                // Bare atom: a number or a variable (lenient — accept `2` or `x`)
                let a = a.clone();
                self.pos += 1;
                if let Ok(n) = a.parse::<f64>() {
                    Ok(Expr::num(n))
                } else {
                    Ok(Expr::var(a))
                }
            }
            Some(SexprTok::Close) => Err(MathError::Parse("unexpected ')'".into())),
            None => Err(MathError::Parse("unexpected end of S-expr".into())),
        }
    }
}

fn build_node(head: &str, args: Vec<Expr>) -> Result<Expr> {
    match head {
        "num" => {
            if args.len() != 1 {
                return Err(MathError::Parse("(num ...) needs one atom".into()));
            }
            match args.into_iter().next().unwrap() {
                Expr::Var(s) => parse_num(&s).map(Expr::num),
                Expr::Num(n) => Ok(Expr::num(n)),
                _ => Err(MathError::Parse("(num ...) needs a numeric atom".into())),
            }
        }
        "var" => {
            if args.len() != 1 {
                return Err(MathError::Parse("(var ...) needs one atom".into()));
            }
            match args.into_iter().next().unwrap() {
                Expr::Var(s) => Ok(Expr::var(s)),
                Expr::Num(n) if n.fract() == 0.0 => Ok(Expr::var(format!("{}", n as i64))),
                _ => Err(MathError::Parse("(var ...) needs a name atom".into())),
            }
        }
        "neg" => unary(args, Expr::neg),
        "add" => binary(args, Expr::add),
        "sub" => binary(args, Expr::sub),
        "mul" => binary(args, Expr::mul),
        "div" => binary(args, Expr::div),
        "pow" => binary(args, Expr::pow),
        "func" => {
            if args.is_empty() {
                return Err(MathError::Parse("(func ...) needs name + args".into()));
            }
            let mut it = args.into_iter();
            let name = match it.next().unwrap() {
                Expr::Var(s) => s,
                Expr::Num(n) if n.fract() == 0.0 => format!("{}", n as i64),
                other => return Err(MathError::Parse(format!("(func ...) bad name: {}", other))),
            };
            let rest: Vec<Expr> = it.collect();
            Ok(Expr::func(name, rest))
        }
        other => Err(MathError::Parse(format!("unknown S-expr head: {}", other))),
    }
}

fn unary(args: Vec<Expr>, f: fn(Expr) -> Expr) -> Result<Expr> {
    if args.len() != 1 {
        return Err(MathError::Parse("unary node needs one arg".into()));
    }
    Ok(f(args.into_iter().next().unwrap()))
}

fn binary(args: Vec<Expr>, f: fn(Expr, Expr) -> Expr) -> Result<Expr> {
    if args.len() != 2 {
        return Err(MathError::Parse("binary node needs two args".into()));
    }
    let mut it = args.into_iter();
    let a = it.next().unwrap();
    let b = it.next().unwrap();
    Ok(f(a, b))
}

// =========================================================================
// JSON
// =========================================================================

/// Convert an [`Expr`] to a JSON string of nested objects.
///
/// Schema: `{"t":"num","v":2}`, `{"t":"var","v":"x"}`, `{"t":"neg","e":{...}}`,
/// `{"t":"add","a":{...},"b":{...}}` (likewise `sub`/`mul`/`div`/`pow`),
/// `{"t":"func","n":"sin","a":[{...},...]}`. Non-finite numbers serialize
/// their `v` as a string (`"NaN"`, `"inf"`, `"-inf"`).
pub fn to_json(e: &Expr) -> String {
    let mut s = String::new();
    write_json(e, &mut s);
    s
}

fn write_json(e: &Expr, out: &mut String) {
    match e {
        Expr::Num(n) => {
            if n.is_finite() {
                out.push_str(&format!("{{\"t\":\"num\",\"v\":{}}}", fmt_num(*n)));
            } else {
                out.push_str(&format!("{{\"t\":\"num\",\"v\":\"{}\"}}", fmt_num(*n)));
            }
        }
        Expr::Var(v) => {
            out.push_str("{\"t\":\"var\",\"v\":");
            json_string(v, out);
            out.push('}');
        }
        Expr::Neg(a) => {
            out.push_str("{\"t\":\"neg\",\"e\":");
            write_json(a, out);
            out.push('}');
        }
        Expr::Add(a, b) => {
            out.push_str("{\"t\":\"add\",\"a\":");
            write_json(a, out);
            out.push_str(",\"b\":");
            write_json(b, out);
            out.push('}');
        }
        Expr::Sub(a, b) => {
            out.push_str("{\"t\":\"sub\",\"a\":");
            write_json(a, out);
            out.push_str(",\"b\":");
            write_json(b, out);
            out.push('}');
        }
        Expr::Mul(a, b) => {
            out.push_str("{\"t\":\"mul\",\"a\":");
            write_json(a, out);
            out.push_str(",\"b\":");
            write_json(b, out);
            out.push('}');
        }
        Expr::Div(a, b) => {
            out.push_str("{\"t\":\"div\",\"a\":");
            write_json(a, out);
            out.push_str(",\"b\":");
            write_json(b, out);
            out.push('}');
        }
        Expr::Pow(a, b) => {
            out.push_str("{\"t\":\"pow\",\"a\":");
            write_json(a, out);
            out.push_str(",\"b\":");
            write_json(b, out);
            out.push('}');
        }
        Expr::Func(name, args) => {
            out.push_str("{\"t\":\"func\",\"n\":");
            json_string(name, out);
            out.push_str(",\"a\":[");
            for (i, a) in args.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                write_json(a, out);
            }
            out.push_str("]}");
        }
    }
}

fn json_string(s: &str, out: &mut String) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

/// Parse a JSON string produced by [`to_json`] back into an [`Expr`].
pub fn from_json(input: &str) -> Result<Expr> {
    let mut p = JsonParser {
        bytes: input.as_bytes(),
        pos: 0,
    };
    p.skip_ws();
    let v = p.parse_value()?;
    p.skip_ws();
    if p.pos < p.bytes.len() {
        return Err(MathError::Parse(format!(
            "unexpected trailing JSON at byte {}",
            p.pos
        )));
    }
    json_value_to_expr(&v)
}

#[derive(Debug)]
enum JsonValue {
    Null,
    #[allow(dead_code)]
    Bool(bool),
    Num(f64),
    Str(String),
    Array(Vec<JsonValue>),
    Object(Vec<(String, JsonValue)>),
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

    fn parse_value(&mut self) -> Result<JsonValue> {
        self.skip_ws();
        match self.bytes.get(self.pos) {
            None => Err(MathError::Parse("unexpected end of JSON".into())),
            Some(b'{') => self.parse_object(),
            Some(b'[') => self.parse_array(),
            Some(b'"') => self.parse_string().map(JsonValue::Str),
            Some(b't') | Some(b'f') => self.parse_bool(),
            Some(b'n') => self.parse_null(),
            Some(c) if *c == b'-' || *c == b'+' || c.is_ascii_digit() => self.parse_number(),
            Some(_) => Err(MathError::Parse(format!(
                "unexpected JSON char at byte {}",
                self.pos
            ))),
        }
    }

    fn parse_object(&mut self) -> Result<JsonValue> {
        self.pos += 1; // {
        let mut entries = Vec::new();
        self.skip_ws();
        if self.bytes.get(self.pos) == Some(&b'}') {
            self.pos += 1;
            return Ok(JsonValue::Object(entries));
        }
        loop {
            self.skip_ws();
            let key = self.parse_string()?;
            self.skip_ws();
            if self.bytes.get(self.pos) != Some(&b':') {
                return Err(MathError::Parse("expected ':' in JSON object".into()));
            }
            self.pos += 1;
            let val = self.parse_value()?;
            entries.push((key, val));
            self.skip_ws();
            match self.bytes.get(self.pos) {
                Some(b',') => {
                    self.pos += 1;
                }
                Some(b'}') => {
                    self.pos += 1;
                    break;
                }
                _ => return Err(MathError::Parse("expected ',' or '}' in JSON object".into())),
            }
        }
        Ok(JsonValue::Object(entries))
    }

    fn parse_array(&mut self) -> Result<JsonValue> {
        self.pos += 1; // [
        let mut items = Vec::new();
        self.skip_ws();
        if self.bytes.get(self.pos) == Some(&b']') {
            self.pos += 1;
            return Ok(JsonValue::Array(items));
        }
        loop {
            let v = self.parse_value()?;
            items.push(v);
            self.skip_ws();
            match self.bytes.get(self.pos) {
                Some(b',') => {
                    self.pos += 1;
                }
                Some(b']') => {
                    self.pos += 1;
                    break;
                }
                _ => return Err(MathError::Parse("expected ',' or ']' in JSON array".into())),
            }
        }
        Ok(JsonValue::Array(items))
    }

    fn parse_string(&mut self) -> Result<String> {
        if self.bytes.get(self.pos) != Some(&b'"') {
            return Err(MathError::Parse("expected '\"' to start JSON string".into()));
        }
        self.pos += 1;
        let mut s = String::new();
        while self.pos < self.bytes.len() {
            let c = self.bytes[self.pos];
            self.pos += 1;
            match c {
                b'"' => return Ok(s),
                b'\\' => {
                    let esc = self
                        .bytes
                        .get(self.pos)
                        .ok_or_else(|| MathError::Parse("unterminated JSON escape".into()))?;
                    self.pos += 1;
                    match *esc {
                        b'"' => s.push('"'),
                        b'\\' => s.push('\\'),
                        b'/' => s.push('/'),
                        b'n' => s.push('\n'),
                        b't' => s.push('\t'),
                        b'r' => s.push('\r'),
                        b'b' => s.push('\u{0008}'),
                        b'f' => s.push('\u{000C}'),
                        b'u' => {
                            if self.pos + 4 > self.bytes.len() {
                                return Err(MathError::Parse("bad \\u escape in JSON".into()));
                            }
                            let hex = std::str::from_utf8(&self.bytes[self.pos..self.pos + 4])
                                .map_err(|_| MathError::Parse("bad \\u escape in JSON".into()))?;
                            let code = u32::from_str_radix(hex, 16)
                                .map_err(|_| MathError::Parse("bad \\u escape in JSON".into()))?;
                            self.pos += 4;
                            if let Some(ch) = char::from_u32(code) {
                                s.push(ch);
                            }
                        }
                        _ => return Err(MathError::Parse("bad JSON escape".into())),
                    }
                }
                _ => {
                    // Read a full UTF-8 char starting at this byte
                    if c < 0x80 {
                        s.push(c as char);
                    } else {
                        // Back up and decode the multi-byte sequence
                        self.pos -= 1;
                        let rest = &self.bytes[self.pos..];
                        let mut buf = [0u8; 4];
                        let len = std::cmp::min(4, rest.len());
                        buf[..len].copy_from_slice(&rest[..len]);
                        match std::str::from_utf8(&buf[..len]) {
                            Ok(decoded) if !decoded.is_empty() => {
                                let ch = decoded.chars().next().unwrap();
                                let ch_len = ch.len_utf8();
                                s.push(ch);
                                self.pos += ch_len;
                            }
                            _ => {
                                // Truncated; consume one byte to make progress
                                s.push(c as char);
                            }
                        }
                    }
                }
            }
        }
        Err(MathError::Parse("unterminated JSON string".into()))
    }

    fn parse_bool(&mut self) -> Result<JsonValue> {
        if self.bytes[self.pos..].starts_with(b"true") {
            self.pos += 4;
            Ok(JsonValue::Bool(true))
        } else if self.bytes[self.pos..].starts_with(b"false") {
            self.pos += 5;
            Ok(JsonValue::Bool(false))
        } else {
            Err(MathError::Parse("bad JSON bool".into()))
        }
    }

    fn parse_null(&mut self) -> Result<JsonValue> {
        if self.bytes[self.pos..].starts_with(b"null") {
            self.pos += 4;
            Ok(JsonValue::Null)
        } else {
            Err(MathError::Parse("bad JSON null".into()))
        }
    }

    fn parse_number(&mut self) -> Result<JsonValue> {
        let start = self.pos;
        if self.bytes.get(self.pos) == Some(&b'-') || self.bytes.get(self.pos) == Some(&b'+') {
            self.pos += 1;
        }
        while self.pos < self.bytes.len()
            && (self.bytes[self.pos].is_ascii_digit()
                || matches!(self.bytes[self.pos], b'.' | b'e' | b'E' | b'+' | b'-'))
        {
            self.pos += 1;
        }
        let txt = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| MathError::Parse("bad JSON number".into()))?;
        txt.parse::<f64>()
            .map(JsonValue::Num)
            .map_err(|_| MathError::Parse(format!("invalid JSON number: {}", txt)))
    }
}

fn json_value_to_expr(v: &JsonValue) -> Result<Expr> {
    match v {
        JsonValue::Object(entries) => {
            let t = json_field_str(entries, "t")?;
            match t.as_str() {
                "num" => {
                    let v = json_field_value(entries, "v")?;
                    match v {
                        JsonValue::Num(n) => Ok(Expr::num(*n)),
                        JsonValue::Str(s) => parse_num(s).map(Expr::num),
                        _ => Err(MathError::Parse("num.v must be number or string".into())),
                    }
                }
                "var" => {
                    let v = json_field_str(entries, "v")?;
                    Ok(Expr::var(v))
                }
                "neg" => {
                    let e = json_field_value(entries, "e")?;
                    Ok(Expr::neg(json_value_to_expr(e)?))
                }
                "add" => json_binary(entries, Expr::add),
                "sub" => json_binary(entries, Expr::sub),
                "mul" => json_binary(entries, Expr::mul),
                "div" => json_binary(entries, Expr::div),
                "pow" => json_binary(entries, Expr::pow),
                "func" => {
                    let name = json_field_str(entries, "n")?;
                    let a = json_field_value(entries, "a")?;
                    match a {
                        JsonValue::Array(items) => {
                            let args: Result<Vec<Expr>> =
                                items.iter().map(json_value_to_expr).collect();
                            Ok(Expr::func(name, args?))
                        }
                        _ => Err(MathError::Parse("func.a must be an array".into())),
                    }
                }
                other => Err(MathError::Parse(format!("unknown JSON node type: {}", other))),
            }
        }
        _ => Err(MathError::Parse("JSON expression must be an object".into())),
    }
}

fn json_binary(entries: &[(String, JsonValue)], f: fn(Expr, Expr) -> Expr) -> Result<Expr> {
    let a = json_field_value(entries, "a")?;
    let b = json_field_value(entries, "b")?;
    Ok(f(json_value_to_expr(a)?, json_value_to_expr(b)?))
}

fn json_field_str(entries: &[(String, JsonValue)], key: &str) -> Result<String> {
    for (k, v) in entries {
        if k == key {
            if let JsonValue::Str(s) = v {
                return Ok(s.clone());
            }
            return Err(MathError::Parse(format!("JSON field '{}' must be a string", key)));
        }
    }
    Err(MathError::Parse(format!("missing JSON field '{}'", key)))
}

fn json_field_value<'a>(
    entries: &'a [(String, JsonValue)],
    key: &str,
) -> Result<&'a JsonValue> {
    for (k, v) in entries {
        if k == key {
            return Ok(v);
        }
    }
    Err(MathError::Parse(format!("missing JSON field '{}'", key)))
}

// =========================================================================
// RPN (postfix)
// =========================================================================

/// Convert an [`Expr`] to a space-separated postfix (RPN) string.
///
/// Operands are pushed first, then the operator. Binary ops: `a b +`.
/// Unary neg: `a neg`. Functions use a `<name>:<arity>` call token, e.g.
/// `x sin:1` for `sin(x)`, `x y pow:2` for `pow(x, y)`.
pub fn to_rpn(e: &Expr) -> String {
    let mut toks: Vec<String> = Vec::new();
    write_rpn(e, &mut toks);
    toks.join(" ")
}

fn write_rpn(e: &Expr, out: &mut Vec<String>) {
    match e {
        Expr::Num(n) => out.push(fmt_num(*n)),
        Expr::Var(v) => out.push(v.clone()),
        Expr::Neg(a) => {
            write_rpn(a, out);
            out.push("neg".to_string());
        }
        Expr::Add(a, b) => rpn_binary(a, b, "+", out),
        Expr::Sub(a, b) => rpn_binary(a, b, "-", out),
        Expr::Mul(a, b) => rpn_binary(a, b, "*", out),
        Expr::Div(a, b) => rpn_binary(a, b, "/", out),
        Expr::Pow(a, b) => rpn_binary(a, b, "^", out),
        Expr::Func(name, args) => {
            for a in args {
                write_rpn(a, out);
            }
            out.push(format!("{}:{}", name, args.len()));
        }
    }
}

fn rpn_binary(a: &Expr, b: &Expr, op: &str, out: &mut Vec<String>) {
    write_rpn(a, out);
    write_rpn(b, out);
    out.push(op.to_string());
}

/// Parse an RPN string produced by [`to_rpn`] back into an [`Expr`].
pub fn from_rpn(input: &str) -> Result<Expr> {
    let mut stack: Vec<Expr> = Vec::new();
    for tok in input.split_whitespace() {
        if tok.is_empty() {
            continue;
        }
        match tok {
            "+" => rpn_pop2(&mut stack, Expr::add)?,
            "-" => rpn_pop2(&mut stack, Expr::sub)?,
            "*" => rpn_pop2(&mut stack, Expr::mul)?,
            "/" => rpn_pop2(&mut stack, Expr::div)?,
            "^" => rpn_pop2(&mut stack, Expr::pow)?,
            "neg" => {
                let a = stack
                    .pop()
                    .ok_or_else(|| MathError::Parse("RPN 'neg' on empty stack".into()))?;
                stack.push(Expr::neg(a));
            }
            other => {
                if let Some((name, arity_str)) = other.rsplit_once(':') {
                    if let Ok(arity) = arity_str.parse::<usize>() {
                        if stack.len() < arity {
                            return Err(MathError::Parse(format!(
                                "RPN call '{}' needs {} args",
                                other, arity
                            )));
                        }
                        let mut args: Vec<Expr> = (0..arity).map(|_| stack.pop().unwrap()).collect();
                        args.reverse(); // pops come out in reverse order; restore left-to-right
                        stack.push(Expr::func(name.to_string(), args));
                        continue;
                    }
                }
                // Operand: number or variable
                if let Ok(n) = tok.parse::<f64>() {
                    stack.push(Expr::num(n));
                } else if tok == "NaN" || tok == "inf" || tok == "-inf" {
                    stack.push(Expr::num(parse_num(tok)?));
                } else {
                    stack.push(Expr::var(tok.to_string()));
                }
            }
        }
    }
    if stack.len() == 1 {
        Ok(stack.pop().unwrap())
    } else if stack.is_empty() {
        Err(MathError::Parse("RPN produced empty expression".into()))
    } else {
        Err(MathError::Parse(format!(
            "RPN left {} items on stack",
            stack.len()
        )))
    }
}

fn rpn_pop2(stack: &mut Vec<Expr>, f: fn(Expr, Expr) -> Expr) -> Result<()> {
    if stack.len() < 2 {
        return Err(MathError::Parse("RPN binary op on short stack".into()));
    }
    let b = stack.pop().unwrap();
    let a = stack.pop().unwrap();
    stack.push(f(a, b));
    Ok(())
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use crate::expr::Expr;

    fn sample() -> Expr {
        // 2*x + sin(x) - 1
        Expr::sub(
            Expr::add(
                Expr::mul(Expr::num(2.0), Expr::var("x")),
                Expr::func("sin", vec![Expr::var("x")]),
            ),
            Expr::num(1.0),
        )
    }

    fn deep() -> Expr {
        // (x + 1) ^ (y * 2) / sqrt(z)
        Expr::div(
            Expr::pow(
                Expr::add(Expr::var("x"), Expr::num(1.0)),
                Expr::mul(Expr::var("y"), Expr::num(2.0)),
            ),
            Expr::func("sqrt", vec![Expr::var("z")]),
        )
    }

    fn multi_arg() -> Expr {
        // pow(x, 2) + max(a, b)
        Expr::add(
            Expr::func("pow", vec![Expr::var("x"), Expr::num(2.0)]),
            Expr::func("max", vec![Expr::var("a"), Expr::var("b")]),
        )
    }

    // ---- S-expr ----

    #[test]
    fn sexpr_basic() {
        let e = Expr::add(Expr::mul(Expr::num(2.0), Expr::var("x")), Expr::num(1.0));
        assert_eq!(to_sexpr(&e), "(add (mul (num 2) (var x)) (num 1))");
    }

    #[test]
    fn sexpr_roundtrip_sample() {
        let e = sample();
        let s = to_sexpr(&e);
        let e2 = from_sexpr(&s).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_roundtrip_deep() {
        let e = deep();
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_roundtrip_multi_arg() {
        let e = multi_arg();
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_neg_and_func() {
        let e = Expr::neg(Expr::func("sin", vec![Expr::var("x")]));
        assert_eq!(to_sexpr(&e), "(neg (func sin (var x)))");
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_non_finite() {
        let cases = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
        for n in cases {
            let e = Expr::num(n);
            let s = to_sexpr(&e);
            let e2 = from_sexpr(&s).unwrap();
            match (e, e2) {
                (Expr::Num(a), Expr::Num(b)) => {
                    assert_eq!(a.is_nan(), b.is_nan());
                    if !a.is_nan() {
                        assert_eq!(a, b);
                    }
                }
                _ => panic!("expected Num"),
            }
        }
    }

    #[test]
    fn sexpr_bad_input() {
        assert!(from_sexpr("(add (num 1))").is_err()); // wrong arity
        assert!(from_sexpr("(bogus (num 1))").is_err()); // unknown head
        assert!(from_sexpr("(num 1").is_err()); // unclosed
    }

    // ---- JSON ----

    #[test]
    fn json_basic() {
        let e = Expr::add(Expr::mul(Expr::num(2.0), Expr::var("x")), Expr::num(1.0));
        let j = to_json(&e);
        assert!(j.contains("\"t\":\"add\""));
        assert!(j.contains("\"t\":\"mul\""));
        assert!(j.contains("\"v\":\"x\""));
        assert!(j.contains("\"v\":2"));
    }

    #[test]
    fn json_roundtrip_sample() {
        let e = sample();
        let e2 = from_json(&to_json(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_roundtrip_deep() {
        let e = deep();
        let e2 = from_json(&to_json(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_roundtrip_multi_arg() {
        let e = multi_arg();
        let e2 = from_json(&to_json(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_non_finite() {
        for n in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let e = Expr::num(n);
            let j = to_json(&e);
            assert!(j.contains("\"v\":\""));
            let e2 = from_json(&j).unwrap();
            match (e, e2) {
                (Expr::Num(a), Expr::Num(b)) => {
                    assert_eq!(a.is_nan(), b.is_nan());
                    if !a.is_nan() {
                        assert_eq!(a, b);
                    }
                }
                _ => panic!("expected Num"),
            }
        }
    }

    #[test]
    fn json_string_escaping() {
        // Variable name with a quote and newline — must survive round-trip.
        let e = Expr::var("a\"b\nc");
        let j = to_json(&e);
        assert!(j.contains("\\\""));
        assert!(j.contains("\\n"));
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_bad_input() {
        assert!(from_json("{\"t\":\"add\",\"a\":{\"t\":\"num\",\"v\":1}}").is_err()); // missing b
        assert!(from_json("{\"t\":\"bogus\"}").is_err());
        assert!(from_json("not json").is_err());
    }

    // ---- RPN ----

    #[test]
    fn rpn_basic() {
        // 2*x + 1 → "2 x * 1 +"
        let e = Expr::add(Expr::mul(Expr::num(2.0), Expr::var("x")), Expr::num(1.0));
        assert_eq!(to_rpn(&e), "2 x * 1 +");
    }

    #[test]
    fn rpn_roundtrip_sample() {
        let e = sample();
        let e2 = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_roundtrip_deep() {
        let e = deep();
        let e2 = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_roundtrip_multi_arg() {
        let e = multi_arg();
        let r = to_rpn(&e);
        assert!(r.contains("pow:2"));
        assert!(r.contains("max:2"));
        let e2 = from_rpn(&r).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_neg_and_func() {
        let e = Expr::neg(Expr::func("sin", vec![Expr::var("x")]));
        assert_eq!(to_rpn(&e), "x sin:1 neg");
        let e2 = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_non_finite() {
        for n in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let e = Expr::num(n);
            let e2 = from_rpn(&to_rpn(&e)).unwrap();
            match (e, e2) {
                (Expr::Num(a), Expr::Num(b)) => {
                    assert_eq!(a.is_nan(), b.is_nan());
                    if !a.is_nan() {
                        assert_eq!(a, b);
                    }
                }
                _ => panic!("expected Num"),
            }
        }
    }

    #[test]
    fn rpn_bad_input() {
        assert!(from_rpn("1 +").is_err()); // short stack
        assert!(from_rpn("1 2 + 3").is_err()); // leftover
        assert!(from_rpn("a b pow:3").is_err()); // arity mismatch
    }

    // ---- Cross-format consistency ----

    #[test]
    fn all_formats_agree() {
        let e = sample();
        let from_s = from_sexpr(&to_sexpr(&e)).unwrap();
        let from_j = from_json(&to_json(&e)).unwrap();
        let from_r = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&from_s));
        assert!(e.equals(&from_j));
        assert!(e.equals(&from_r));
        // And all three parsed forms agree with each other
        assert!(from_s.equals(&from_j));
        assert!(from_j.equals(&from_r));
    }

    // ---- Additional edge cases ----

    fn nested_funcs() -> Expr {
        // sin(cos(x^2)) + exp(-x)
        Expr::add(
            Expr::func("sin", vec![Expr::func("cos", vec![Expr::pow(Expr::var("x"), Expr::num(2.0))])]),
            Expr::func("exp", vec![Expr::neg(Expr::var("x"))]),
        )
    }

    fn division_expr() -> Expr {
        // (a + b) / (c - d)
        Expr::div(
            Expr::add(Expr::var("a"), Expr::var("b")),
            Expr::sub(Expr::var("c"), Expr::var("d")),
        )
    }

    #[test]
    fn sexpr_nested_functions() {
        let e = nested_funcs();
        let s = to_sexpr(&e);
        assert!(s.contains("(func sin"));
        assert!(s.contains("(func cos"));
        assert!(s.contains("(func exp"));
        let e2 = from_sexpr(&s).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_division() {
        let e = division_expr();
        let s = to_sexpr(&e);
        assert!(s.contains("(div "));
        let e2 = from_sexpr(&s).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_negative_number() {
        let e = Expr::num(-5.0);
        assert_eq!(to_sexpr(&e), "(num -5)");
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_decimal_number() {
        let e = Expr::num(3.14);
        assert_eq!(to_sexpr(&e), "(num 3.14)");
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_bare_atom_number() {
        // Lenient path: bare "42" should parse as Num(42)
        let e = from_sexpr("42").unwrap();
        assert!(e.equals(&Expr::num(42.0)));
    }

    #[test]
    fn sexpr_bare_atom_variable() {
        // Lenient path: bare "x" should parse as Var("x")
        let e = from_sexpr("x").unwrap();
        assert!(e.equals(&Expr::var("x")));
    }

    #[test]
    fn sexpr_whitespace_tolerance() {
        let e = sample();
        // Extra spaces between tokens should still parse
        let s = to_sexpr(&e).replace(' ', "  ");
        let e2 = from_sexpr(&s).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn sexpr_func_no_args() {
        let e = Expr::func("foo", vec![]);
        assert_eq!(to_sexpr(&e), "(func foo)");
        let e2 = from_sexpr(&to_sexpr(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_nested_functions() {
        let e = nested_funcs();
        let j = to_json(&e);
        assert!(j.contains("\"t\":\"func\""));
        assert!(j.contains("\"n\":\"sin\""));
        assert!(j.contains("\"n\":\"cos\""));
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_division() {
        let e = division_expr();
        let j = to_json(&e);
        assert!(j.contains("\"t\":\"div\""));
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_func_no_args() {
        let e = Expr::func("foo", vec![]);
        let j = to_json(&e);
        assert!(j.contains("\"a\":[]"));
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_single_number() {
        let e = Expr::num(42.0);
        let j = to_json(&e);
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn json_single_variable() {
        let e = Expr::var("x");
        let j = to_json(&e);
        let e2 = from_json(&j).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_nested_functions() {
        let e = nested_funcs();
        let r = to_rpn(&e);
        // x 2 ^ cos:1 sin:1  x neg exp:1  +
        assert!(r.contains("cos:1"));
        assert!(r.contains("sin:1"));
        assert!(r.contains("exp:1"));
        let e2 = from_rpn(&r).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_division() {
        let e = division_expr();
        let r = to_rpn(&e);
        // a b + c d - /
        assert!(r.contains("/"));
        let e2 = from_rpn(&r).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_multiple_spaces() {
        let e = sample();
        let r = to_rpn(&e).replace(' ', "  ");
        let e2 = from_rpn(&r).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_func_zero_arity() {
        let e = Expr::func("foo", vec![]);
        let r = to_rpn(&e);
        assert_eq!(r, "foo:0");
        let e2 = from_rpn(&r).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_single_number() {
        let e = Expr::num(42.0);
        assert_eq!(to_rpn(&e), "42");
        let e2 = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_single_variable() {
        let e = Expr::var("x");
        assert_eq!(to_rpn(&e), "x");
        let e2 = from_rpn(&to_rpn(&e)).unwrap();
        assert!(e.equals(&e2));
    }

    #[test]
    fn rpn_empty_input() {
        assert!(from_rpn("").is_err());
        assert!(from_rpn("   ").is_err());
    }

    #[test]
    fn sexpr_empty_input() {
        assert!(from_sexpr("").is_err());
        assert!(from_sexpr("   ").is_err());
    }

    #[test]
    fn json_empty_input() {
        assert!(from_json("").is_err());
        assert!(from_json("   ").is_err());
    }

    // ---- Round-trip preserves evaluation result ----

    #[test]
    fn roundtrip_preserves_eval() {
        use crate::eval::{eval, Context};
        use crate::parser::Parser;
        let original = Parser::parse("2*x + sin(x)").unwrap();
        let ctx = Context::standard().with("x", 1.5);
        let val_orig = eval(&original, &ctx).unwrap();
        // Round-trip through all formats and verify eval matches
        for (text, parse) in [
            (to_sexpr(&original), from_sexpr as fn(&str) -> Result<Expr>),
            (to_json(&original), from_json as fn(&str) -> Result<Expr>),
            (to_rpn(&original), from_rpn as fn(&str) -> Result<Expr>),
        ] {
            let parsed = parse(&text).unwrap();
            let val_parsed = eval(&parsed, &ctx).unwrap();
            assert_abs_diff_eq!(val_orig, val_parsed, epsilon = 1e-12);
        }
    }

    // ---- Deeply nested stress test ----

    #[test]
    fn deeply_nested_stress() {
        // Build ((((x+1)+1)+1)+1) — 5 levels of nesting
        let mut e = Expr::var("x");
        for _ in 0..5 {
            e = Expr::add(e, Expr::num(1.0));
        }
        for (text, parse) in [
            (to_sexpr(&e), from_sexpr as fn(&str) -> Result<Expr>),
            (to_json(&e), from_json as fn(&str) -> Result<Expr>),
            (to_rpn(&e), from_rpn as fn(&str) -> Result<Expr>),
        ] {
            let parsed = parse(&text).unwrap();
            assert!(e.equals(&parsed), "deep nest failed: {}", text);
        }
    }

    #[test]
    fn all_formats_roundtrip_nested_funcs() {
        let e = nested_funcs();
        for (text, parse) in [
            (to_sexpr(&e), from_sexpr as fn(&str) -> Result<Expr>),
            (to_json(&e), from_json as fn(&str) -> Result<Expr>),
            (to_rpn(&e), from_rpn as fn(&str) -> Result<Expr>),
        ] {
            let parsed = parse(&text).unwrap();
            assert!(e.equals(&parsed), "nested funcs round-trip failed: {}", text);
        }
    }
}
