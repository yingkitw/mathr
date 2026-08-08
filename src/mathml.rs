//! MathML (Presentation) export and import for [`Expr`].
//!
//! - [`to_mathml`] converts an `Expr` AST to a W3C Presentation MathML string.
//! - [`from_mathml`] parses a subset of Presentation MathML back into an `Expr`.
//!
//! ## Export example
//! ```
//! use mathr::expr::Expr;
//! use mathr::mathml::to_mathml;
//! let e = Expr::add(Expr::mul(Expr::num(2.0), Expr::var("x")), Expr::num(1.0));
//! let ml = to_mathml(&e);
//! assert!(ml.contains("<mi>x</mi>"));
//! ```
//!
//! ## Supported MathML elements (import)
//! `<mn>`, `<mi>`, `<mo>`, `<mrow>`, `<mfrac>`, `<msup>`, `<msub>`,
//! `<msqrt>`, `<mroot>`, `<mtext>`, `<mstyle>`, `<mfenced>`

use crate::error::{MathError, Result};
use crate::expr::Expr;

// =========================================================================
// Export: Expr → Presentation MathML
// =========================================================================

/// Convert an [`Expr`] to a W3C Presentation MathML string (without the
/// wrapping `<math>` element — use [`to_mathml_doc`] for a full document).
pub fn to_mathml(e: &Expr) -> String {
    let mut out = String::new();
    write_expr(e, &mut out);
    out
}

/// Convert an [`Expr`] to a complete `<math>` document.
pub fn to_mathml_doc(e: &Expr) -> String {
    format!("<math xmlns=\"http://www.w3.org/1998/Math/MathML\">{}</math>", to_mathml(e))
}

fn write_expr(e: &Expr, out: &mut String) {
    match e {
        Expr::Num(n) => {
            if n.is_nan() {
                out.push_str("<mi>NaN</mi>");
            } else if n.is_infinite() {
                if *n > 0.0 {
                    out.push_str("<mi>&#x221E;</mi>"); // ∞
                } else {
                    out.push_str("<mo>-</mo><mi>&#x221E;</mi>");
                }
            } else if n.fract() == 0.0 && n.abs() < 1e15 {
                out.push_str(&format!("<mn>{}</mn>", *n as i64));
            } else {
                out.push_str(&format!("<mn>{}</mn>", n));
            }
        }
        Expr::Var(v) => match v.as_str() {
            "pi" => out.push_str("<mi>&#x03C0;</mi>"), // π
            "e" => out.push_str("<mi>e</mi>"),
            "tau" => out.push_str("<mi>&#x03C4;</mi>"), // τ
            "inf" => out.push_str("<mi>&#x221E;</mi>"), // ∞
            _ => out.push_str(&format!("<mi>{}</mi>", v)),
        },
        Expr::Neg(inner) => {
            out.push_str("<mrow><mo>-</mo>");
            write_expr(inner, out);
            out.push_str("</mrow>");
        }
        Expr::Add(a, b) => {
            out.push_str("<mrow>");
            write_expr(a, out);
            out.push_str("<mo>+</mo>");
            write_expr(b, out);
            out.push_str("</mrow>");
        }
        Expr::Sub(a, b) => {
            out.push_str("<mrow>");
            write_expr(a, out);
            out.push_str("<mo>-</mo>");
            write_expr(b, out);
            out.push_str("</mrow>");
        }
        Expr::Mul(a, b) => {
            out.push_str("<mrow>");
            write_expr(a, out);
            out.push_str("<mo>&#x22C5;</mo>"); // ⋅
            write_expr(b, out);
            out.push_str("</mrow>");
        }
        Expr::Div(a, b) => {
            out.push_str("<mfrac>");
            write_expr(a, out);
            write_expr(b, out);
            out.push_str("</mfrac>");
        }
        Expr::Pow(a, b) => {
            out.push_str("<msup>");
            write_expr(a, out);
            write_expr(b, out);
            out.push_str("</msup>");
        }
        Expr::Func(name, args) => {
            // Special cases for standard MathML function elements
            match (name.as_str(), args.len()) {
                ("sqrt", 1) => {
                    out.push_str("<msqrt>");
                    write_expr(&args[0], out);
                    out.push_str("</msqrt>");
                }
                ("factorial", 1) => {
                    out.push_str("<mrow>");
                    write_expr(&args[0], out);
                    out.push_str("<mo>!</mo>");
                    out.push_str("</mrow>");
                }
                ("abs", 1) => {
                    out.push_str("<mrow><mo>|</mo>");
                    write_expr(&args[0], out);
                    out.push_str("<mo>|</mo></mrow>");
                }
                ("C", 2) => {
                    // Binomial: <mfenced open="(" close=")"><mfrac>n</mfrac>k</mfenced>
                    out.push_str("<mfenced open=\"(\" close=\")\"><mfrac>");
                    write_expr(&args[0], out);
                    write_expr(&args[1], out);
                    out.push_str("</mfrac></mfenced>");
                }
                _ => {
                    out.push_str("<mrow>");
                    out.push_str(&format!("<mi>{}</mi>", name));
                    out.push_str("<mfenced open=\"(\" close=\")\">");
                    for (i, a) in args.iter().enumerate() {
                        if i > 0 {
                            out.push_str("<mo>,</mo>");
                        }
                        write_expr(a, out);
                    }
                    out.push_str("</mfenced>");
                    out.push_str("</mrow>");
                }
            }
        }
    }
}

// =========================================================================
// Import: Presentation MathML → Expr
// =========================================================================

/// Parse a Presentation MathML string into an [`Expr`].
///
/// Supports a practical subset: `<mn>`, `<mi>`, `<mo>`, `<mrow>`,
/// `<mfrac>`, `<msup>`, `<msub>`, `<msqrt>`, `<mroot>`,
/// `<mtext>`, `<mstyle>`, `<mfenced>`.
pub fn from_mathml(input: &str) -> Result<Expr> {
    let tokens = tokenize(input)?;
    let mut parser = MathMLParser { tokens, pos: 0 };
    let result = parser.parse_element()?;
    if parser.pos < parser.tokens.len() {
        return Err(MathError::Parse(format!(
            "unexpected trailing MathML at token {}",
            parser.pos
        )));
    }
    Ok(result)
}

#[derive(Debug, Clone)]
enum Token {
    Open(String),  // <tag ...>
    Close(String), // </tag>
    Text(String),  // text content
    SelfClose(String), // <tag .../>
}

fn tokenize(input: &str) -> Result<Vec<Token>> {
    let mut tokens = Vec::new();
    let bytes = input.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i].is_ascii_whitespace() {
            i += 1;
            continue;
        }
        if bytes[i] == b'<' {
            // Comment?
            if input[i..].starts_with("<!--") {
                if let Some(end) = input[i + 4..].find("-->") {
                    i += 4 + end + 3;
                    continue;
                }
                return Err(MathError::Parse("unterminated comment in MathML".into()));
            }
            // Closing tag?
            if input[i..].starts_with("</") {
                let end = input[i + 2..]
                    .find('>')
                    .ok_or_else(|| MathError::Parse("unterminated closing tag".into()))?;
                let tag = input[i + 2..i + 2 + end].trim().to_string();
                tokens.push(Token::Close(tag));
                i += 2 + end + 1;
                continue;
            }
            // Opening or self-closing tag
            let end = input[i + 1..]
                .find('>')
                .ok_or_else(|| MathError::Parse("unterminated opening tag".into()))?;
            let raw = &input[i + 1..i + 1 + end];
            let self_closing = raw.ends_with('/');
            let raw = raw.trim_end_matches('/').trim();
            // Extract tag name (first word, ignoring attributes)
            let tag_name = raw
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            if tag_name.is_empty() {
                return Err(MathError::Parse("empty tag name in MathML".into()));
            }
            if self_closing {
                tokens.push(Token::SelfClose(tag_name));
            } else {
                tokens.push(Token::Open(tag_name));
            }
            i += 1 + end + 1;
        } else {
            // Text content
            let start = i;
            while i < bytes.len() && bytes[i] != b'<' {
                i += 1;
            }
            let text = input[start..i].trim();
            if !text.is_empty() {
                tokens.push(Token::Text(decode_entities(text)));
            }
        }
    }
    Ok(tokens)
}

fn decode_entities(s: &str) -> String {
    s.replace("&#x221E;", "inf")
        .replace("&#x03C0;", "pi")
        .replace("&#x03C4;", "tau")
        .replace("&#x22C5;", "*")
        .replace("&infin;", "inf")
        .replace("&pi;", "pi")
        .replace("&tau;", "tau")
        .replace("&times;", "*")
        .replace("&minus;", "-")
        .replace("&plusmn;", "+-")
        .replace("&amp;", "&")
        .replace("&lt;", "<")
        .replace("&gt;", ">")
}

struct MathMLParser {
    tokens: Vec<Token>,
    pos: usize,
}

impl MathMLParser {
    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn parse_element(&mut self) -> Result<Expr> {
        match self.peek() {
            Some(Token::Open(tag)) | Some(Token::SelfClose(tag)) => {
                let tag = tag.clone();
                let is_self_close = matches!(self.peek(), Some(Token::SelfClose(_)));
                self.pos += 1;
                if is_self_close {
                    return self.parse_self_close(&tag);
                }
                let result = self.parse_tag_body(&tag)?;
                // Expect closing tag
                if let Some(Token::Close(close_tag)) = self.peek() {
                    if close_tag == &tag {
                        self.pos += 1;
                        return Ok(result);
                    }
                }
                Err(MathError::Parse(format!(
                    "expected </{}> at token {}",
                    tag, self.pos
                )))
            }
            Some(Token::Text(_)) => {
                // Bare text — try to parse as expression
                if let Token::Text(t) = &self.tokens[self.pos] {
                    self.pos += 1;
                    crate::parser::Parser::parse(t)
                } else {
                    unreachable!()
                }
            }
            _ => Err(MathError::Parse(format!(
                "unexpected token at {}",
                self.pos
            ))),
        }
    }

    fn parse_self_close(&mut self, tag: &str) -> Result<Expr> {
        match tag {
            "mn" | "mi" | "mo" | "mtext" => {
                Err(MathError::Parse(format!("<{}> cannot be self-closing", tag)))
            }
            _ => Ok(Expr::num(0.0)), // empty element → 0
        }
    }

    fn parse_tag_body(&mut self, tag: &str) -> Result<Expr> {
        match tag {
            "math" => self.parse_element(),
            "mn" => self.parse_text_as_number(),
            "mi" => self.parse_text_as_var(),
            "mo" => self.parse_text_as_op(),
            "mtext" => self.parse_text_as_var(),
            "mrow" => self.parse_mrow(),
            "mfrac" => self.parse_mfrac(),
            "msup" => self.parse_msup(),
            "msub" => self.parse_msub(),
            "msqrt" => self.parse_msqrt(),
            "mroot" => self.parse_mroot(),
            "mstyle" => self.parse_element(),
            "mfenced" => self.parse_mfenced(),
            _ => {
                // Unknown tag — try to parse its children as a sequence
                self.parse_mrow()
            }
        }
    }

    fn parse_text_as_number(&mut self) -> Result<Expr> {
        let text = match self.peek() {
            Some(Token::Text(t)) => t.clone(),
            _ => return Err(MathError::Parse("expected text in <mn>".into())),
        };
        self.pos += 1;
        let n: f64 = text
            .parse()
            .map_err(|_| MathError::Parse(format!("invalid number: {}", text)))?;
        Ok(Expr::num(n))
    }

    fn parse_text_as_var(&mut self) -> Result<Expr> {
        let text = match self.peek() {
            Some(Token::Text(t)) => t.clone(),
            _ => return Ok(Expr::num(0.0)), // Empty <mi></mi>
        };
        self.pos += 1;
        match text.as_str() {
            "pi" | "π" => Ok(Expr::var("pi")),
            "tau" | "τ" => Ok(Expr::var("tau")),
            "e" => Ok(Expr::var("e")),
            "inf" | "∞" => Ok(Expr::num(f64::INFINITY)),
            _ => Ok(Expr::var(text)),
        }
    }

    fn parse_text_as_op(&mut self) -> Result<Expr> {
        let text = match self.peek() {
            Some(Token::Text(t)) => t.clone(),
            _ => return Err(MathError::Parse("expected text in <mo>".into())),
        };
        self.pos += 1;
        Ok(Expr::var(text))
    }

    fn parse_mrow(&mut self) -> Result<Expr> {
        let mut children = Vec::new();
        while !matches!(self.peek(), Some(Token::Close(_))) {
            if self.peek().is_none() {
                break;
            }
            children.push(self.parse_element()?);
        }
        if children.is_empty() {
            return Ok(Expr::num(0.0));
        }
        if children.len() == 1 {
            return Ok(children.into_iter().next().unwrap());
        }
        // Build expression from sequence of operands and operators
        build_sequence(children)
    }

    fn parse_mfrac(&mut self) -> Result<Expr> {
        let num = self.parse_element()?;
        let den = self.parse_element()?;
        Ok(Expr::div(num, den))
    }

    fn parse_msup(&mut self) -> Result<Expr> {
        let base = self.parse_element()?;
        let exp = self.parse_element()?;
        Ok(Expr::pow(base, exp))
    }

    fn parse_msub(&mut self) -> Result<Expr> {
        // Subscript: treat as multiplication by a variable with subscript name
        // e.g. x_1 → var "x_1"
        let base = self.parse_element()?;
        let sub = self.parse_element()?;
        if let (Expr::Var(name), Expr::Num(n)) = (&base, &sub) {
            if n.fract() == 0.0 {
                return Ok(Expr::var(format!("{}_{}", name, *n as i64)));
            }
        }
        // Fallback: treat subscript as implicit multiplication
        Ok(Expr::mul(base, sub))
    }

    fn parse_msqrt(&mut self) -> Result<Expr> {
        let inner = self.parse_element()?;
        Ok(Expr::func("sqrt", vec![inner]))
    }

    fn parse_mroot(&mut self) -> Result<Expr> {
        // <mroot>radicand index</mroot> → pow(radicand, 1/index)
        let radicand = self.parse_element()?;
        let index = self.parse_element()?;
        if let Expr::Num(n) = &index {
            if n.fract() == 0.0 && *n > 0.0 {
                let inv = 1.0 / *n;
                return Ok(Expr::pow(radicand, Expr::num(inv)));
            }
        }
        Ok(Expr::pow(radicand, Expr::div(Expr::num(1.0), index)))
    }

    fn parse_mfenced(&mut self) -> Result<Expr> {
        let mut children = Vec::new();
        while !matches!(self.peek(), Some(Token::Close(_))) {
            if self.peek().is_none() {
                break;
            }
            children.push(self.parse_element()?);
        }
        if children.is_empty() {
            return Ok(Expr::num(0.0));
        }
        if children.len() == 1 {
            return Ok(children.into_iter().next().unwrap());
        }
        // Check if it looks like a function call: <mi>f</mi> <mfenced>...</mfenced>
        // or a binomial: <mfrac> inside
        // For now, treat as a parenthesized sequence
        build_sequence(children)
    }
}

/// Build an expression from a flat sequence of operands and operators.
fn build_sequence(children: Vec<Expr>) -> Result<Expr> {
    // First pass: detect function calls — a Var (function name) immediately
    // followed by another operand (the argument list from <mfenced>) becomes Func.
    // We detect this by looking for Var(name) followed by a non-operator Expr
    // where name is a known function name.
    let known_funcs = [
        "sin", "cos", "tan", "asin", "acos", "atan",
        "sinh", "cosh", "tanh", "exp", "ln", "log", "log2", "log10",
        "sqrt", "cbrt", "abs", "floor", "ceil", "round", "sign",
        "gamma", "erf", "erfc", "sinc", "fract",
        "gcd", "lcm", "C", "factorial", "mod", "min", "max", "pow",
        "bessel_j0", "bessel_j1", "bessel_j",
    ];

    // Separate into operands and operators
    let mut operands: Vec<Expr> = Vec::new();
    let mut ops: Vec<String> = Vec::new();
    let mut i = 0;
    while i < children.len() {
        let child = &children[i];
        match child {
            Expr::Var(s) if matches!(s.as_str(), "+" | "-" | "*" | "/" | "^" | "⋅" | "·") => {
                ops.push(s.clone());
                i += 1;
            }
            Expr::Var(name) if known_funcs.contains(&name.as_str()) => {
                // Check if next child is an operand (not an operator) → function call
                if i + 1 < children.len() {
                    let next = &children[i + 1];
                    let is_op = matches!(next, Expr::Var(s) if matches!(s.as_str(), "+" | "-" | "*" | "/" | "^" | "⋅" | "·"));
                    if !is_op {
                        // Function call: name(arg)
                        // The next operand may be a single arg or a sequence (from mfenced with commas)
                        // For simplicity, treat the next operand as the single argument
                        let arg = children[i + 1].clone();
                        operands.push(Expr::func(name.clone(), vec![arg]));
                        i += 2;
                        continue;
                    }
                }
                operands.push(child.clone());
                i += 1;
            }
            _ => {
                operands.push(child.clone());
                i += 1;
            }
        }
    }
    if operands.is_empty() {
        return Err(MathError::Parse("no operands in MathML sequence".into()));
    }
    if ops.is_empty() {
        // All operands with no operators → implicit multiplication (if >1)
        if operands.len() == 1 {
            return Ok(operands.into_iter().next().unwrap());
        }
        return Ok(operands.into_iter().reduce(Expr::mul).unwrap());
    }
    if operands.len() != ops.len() + 1 {
        return Err(MathError::Parse(format!(
            "operand/op mismatch: {} operands, {} ops",
            operands.len(),
            ops.len()
        )));
    }
    // Apply operators with precedence: ^ > * / > + -
    let mut operands = operands;
    let mut ops = ops;
    // First pass: ^
    let mut i = 0;
    while i < ops.len() {
        if ops[i] == "^" {
            let right = operands.remove(i + 1);
            let left = operands.remove(i);
            operands.insert(i, Expr::pow(left, right));
            ops.remove(i);
        } else {
            i += 1;
        }
    }
    // Second pass: * / ⋅ ·
    let mut i = 0;
    while i < ops.len() {
        if ops[i] == "*" || ops[i] == "⋅" || ops[i] == "·" {
            let right = operands.remove(i + 1);
            let left = operands.remove(i);
            operands.insert(i, Expr::mul(left, right));
            ops.remove(i);
        } else if ops[i] == "/" {
            let right = operands.remove(i + 1);
            let left = operands.remove(i);
            operands.insert(i, Expr::div(left, right));
            ops.remove(i);
        } else {
            i += 1;
        }
    }
    // Third pass: + -
    let mut i = 0;
    while i < ops.len() {
        if ops[i] == "+" {
            let right = operands.remove(i + 1);
            let left = operands.remove(i);
            operands.insert(i, Expr::add(left, right));
            ops.remove(i);
        } else if ops[i] == "-" {
            let right = operands.remove(i + 1);
            let left = operands.remove(i);
            operands.insert(i, Expr::sub(left, right));
            ops.remove(i);
        } else {
            i += 1;
        }
    }
    if operands.len() == 1 {
        Ok(operands.into_iter().next().unwrap())
    } else {
        Ok(operands.into_iter().reduce(Expr::mul).unwrap())
    }
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn export_number() {
        let ml = to_mathml(&Expr::num(42.0));
        assert_eq!(ml, "<mn>42</mn>");
    }

    #[test]
    fn export_var() {
        let ml = to_mathml(&Expr::var("x"));
        assert_eq!(ml, "<mi>x</mi>");
    }

    #[test]
    fn export_pi() {
        let ml = to_mathml(&Expr::var("pi"));
        assert!(ml.contains("&#x03C0;"));
    }

    #[test]
    fn export_add() {
        let e = Expr::add(Expr::var("x"), Expr::num(1.0));
        let ml = to_mathml(&e);
        assert!(ml.contains("<mi>x</mi>"));
        assert!(ml.contains("<mn>1</mn>"));
        assert!(ml.contains("<mo>+</mo>"));
    }

    #[test]
    fn export_mul() {
        let e = Expr::mul(Expr::num(2.0), Expr::var("x"));
        let ml = to_mathml(&e);
        assert!(ml.contains("<mo>&#x22C5;</mo>")); // ⋅
    }

    #[test]
    fn export_frac() {
        let e = Expr::div(Expr::num(1.0), Expr::num(2.0));
        let ml = to_mathml(&e);
        assert!(ml.contains("<mfrac>"));
        assert!(ml.contains("</mfrac>"));
    }

    #[test]
    fn export_pow() {
        let e = Expr::pow(Expr::var("x"), Expr::num(2.0));
        let ml = to_mathml(&e);
        assert!(ml.contains("<msup>"));
    }

    #[test]
    fn export_sqrt() {
        let e = Expr::func("sqrt", vec![Expr::var("x")]);
        let ml = to_mathml(&e);
        assert!(ml.contains("<msqrt>"));
    }

    #[test]
    fn export_factorial() {
        let e = Expr::func("factorial", vec![Expr::num(5.0)]);
        let ml = to_mathml(&e);
        assert!(ml.contains("<mo>!</mo>"));
    }

    #[test]
    fn export_abs() {
        let e = Expr::func("abs", vec![Expr::var("x")]);
        let ml = to_mathml(&e);
        assert!(ml.contains("<mo>|</mo>"));
    }

    #[test]
    fn export_binomial() {
        let e = Expr::func("C", vec![Expr::num(5.0), Expr::num(2.0)]);
        let ml = to_mathml(&e);
        assert!(ml.contains("<mfrac>"));
        assert!(ml.contains("<mfenced"));
    }

    #[test]
    fn export_function() {
        let e = Expr::func("sin", vec![Expr::var("x")]);
        let ml = to_mathml(&e);
        assert!(ml.contains("<mi>sin</mi>"));
        assert!(ml.contains("<mfenced"));
    }

    #[test]
    fn export_doc() {
        let e = Expr::var("x");
        let ml = to_mathml_doc(&e);
        assert!(ml.starts_with("<math"));
        assert!(ml.contains("xmlns"));
        assert!(ml.ends_with("</math>"));
    }

    #[test]
    fn import_number() {
        let e = from_mathml("<mn>42</mn>").unwrap();
        assert_eq!(e, Expr::num(42.0));
    }

    #[test]
    fn import_var() {
        let e = from_mathml("<mi>x</mi>").unwrap();
        assert_eq!(e, Expr::var("x"));
    }

    #[test]
    fn import_pi() {
        let e = from_mathml("<mi>&#x03C0;</mi>").unwrap();
        assert_eq!(e, Expr::var("pi"));
    }

    #[test]
    fn import_frac() {
        let e = from_mathml("<mfrac><mn>1</mn><mn>2</mn></mfrac>").unwrap();
        assert_eq!(e, Expr::div(Expr::num(1.0), Expr::num(2.0)));
    }

    #[test]
    fn import_sup() {
        let e = from_mathml("<msup><mi>x</mi><mn>2</mn></msup>").unwrap();
        assert_eq!(e, Expr::pow(Expr::var("x"), Expr::num(2.0)));
    }

    #[test]
    fn import_sqrt() {
        let e = from_mathml("<msqrt><mi>x</mi></msqrt>").unwrap();
        assert_eq!(e, Expr::func("sqrt", vec![Expr::var("x")]));
    }

    #[test]
    fn import_mrow_add() {
        let ml = "<mrow><mi>x</mi><mo>+</mo><mn>1</mn></mrow>";
        let e = from_mathml(ml).unwrap();
        assert_eq!(e, Expr::add(Expr::var("x"), Expr::num(1.0)));
    }

    #[test]
    fn import_mrow_mul() {
        let ml = "<mrow><mn>2</mn><mo>&#x22C5;</mo><mi>x</mi></mrow>";
        let e = from_mathml(ml).unwrap();
        assert_eq!(e, Expr::mul(Expr::num(2.0), Expr::var("x")));
    }

    #[test]
    fn import_mrow_precedence() {
        // 2 + 3 * x → 2 + (3*x)
        let ml = "<mrow><mn>2</mn><mo>+</mo><mn>3</mn><mo>&#x22C5;</mo><mi>x</mi></mrow>";
        let e = from_mathml(ml).unwrap();
        assert_eq!(
            e,
            Expr::add(Expr::num(2.0), Expr::mul(Expr::num(3.0), Expr::var("x")))
        );
    }

    #[test]
    fn import_full_doc() {
        let ml = "<math xmlns=\"http://www.w3.org/1998/Math/MathML\"><mfrac><mn>1</mn><mn>2</mn></mfrac></math>";
        let e = from_mathml(ml).unwrap();
        assert_eq!(e, Expr::div(Expr::num(1.0), Expr::num(2.0)));
    }

    #[test]
    fn roundtrip_simple() {
        let original = Expr::add(Expr::var("x"), Expr::num(1.0));
        let ml = to_mathml(&original);
        let parsed = from_mathml(&ml).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn roundtrip_frac() {
        let original = Expr::div(Expr::num(1.0), Expr::num(2.0));
        let ml = to_mathml(&original);
        let parsed = from_mathml(&ml).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn roundtrip_pow() {
        let original = Expr::pow(Expr::var("x"), Expr::num(2.0));
        let ml = to_mathml(&original);
        let parsed = from_mathml(&ml).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn roundtrip_sqrt() {
        let original = Expr::func("sqrt", vec![Expr::var("x")]);
        let ml = to_mathml(&original);
        let parsed = from_mathml(&ml).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn roundtrip_mul() {
        let original = Expr::mul(Expr::num(2.0), Expr::var("x"));
        let ml = to_mathml(&original);
        let parsed = from_mathml(&ml).unwrap();
        assert_eq!(parsed, original);
    }

    #[test]
    fn import_mroot() {
        // <mroot><mn>8</mn><mn>3</mn></mroot> → 8^(1/3)
        let e = from_mathml("<mroot><mn>8</mn><mn>3</mn></mroot>").unwrap();
        assert_eq!(e, Expr::pow(Expr::num(8.0), Expr::num(1.0 / 3.0)));
    }

    #[test]
    fn import_msub() {
        let e = from_mathml("<msub><mi>x</mi><mn>1</mn></msub>").unwrap();
        assert_eq!(e, Expr::var("x_1"));
    }
}
