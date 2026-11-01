//! Minimal JSON layer for the MCP server.
//!
//! Hand-rolled (no `serde`) to keep mathr's dependency footprint flat,
//! matching the `crate::serialize` convention. Handles escapes, surrogate
//! pairs, exponents, and maps non-finite numbers to `null`.

use crate::error::{MathError, Result};

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
