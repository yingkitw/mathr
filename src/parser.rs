use std::borrow::Cow;

use crate::error::{MathError, Result};
use crate::expr::Expr;

/// Recursive-descent parser for math expressions.
///
/// Supports:
///   - Numbers (integer and decimal, including scientific notation)
///   - Identifiers and constants (`pi`, `e`, `tau`, `inf`)
///   - Binary operators: `+ - * / ^` with standard precedence,
///     `^` right-associative
///   - Unary minus
///   - Function calls: `sin(x)`, `log(x, 10)`
///   - Parenthesised sub-expressions
///   - Implicit multiplication: `2x`, `3(x+1)`, `(x)(y)`
///   - Whitespace is ignored
///   - LaTeX/TeX input: `\frac{a}{b}`, `\sqrt{x}`, `\sin(x)`,
///     `\pi`, `\cdot`, `\left(`, `\right)`, `^{...}`, `_{...}`
pub struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Parser<'a> {
    pub fn new(src: &'a str) -> Self {
        Self { src: src.as_bytes(), pos: 0 }
    }

    pub fn parse(src: &str) -> Result<Expr> {
        let src = strip_math_delimiters(src);
        let mut p = Parser::new(&src);
        let e = p.parse_expr()?;
        p.skip_ws();
        if p.pos < p.src.len() {
            return Err(MathError::Parse(format!(
                "unexpected trailing characters at byte {}",
                p.pos
            )));
        }
        Ok(e)
    }

    // --- helpers ---------------------------------------------------------

    fn peek(&self) -> Option<u8> {
        self.src.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<u8> {
        let c = self.peek()?;
        self.pos += 1;
        Some(c)
    }

    fn skip_ws(&mut self) {
        while let Some(c) = self.peek() {
            if (c as char).is_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    fn expect(&mut self, c: u8) -> Result<()> {
        self.skip_ws();
        match self.peek() {
            Some(x) if x == c => {
                self.pos += 1;
                Ok(())
            }
            Some(x) => Err(MathError::Parse(format!(
                "expected '{}' but found '{}'",
                c as char, x as char
            ))),
            None => Err(MathError::Parse(format!("expected '{}' but found end of input", c as char))),
        }
    }

    fn read_ident(&mut self) -> String {
        let start = self.pos;
        while let Some(c) = self.peek() {
            if (c as char).is_ascii_alphanumeric() || c == b'_' {
                self.pos += 1;
            } else {
                break;
            }
        }
        std::str::from_utf8(&self.src[start..self.pos])
            .unwrap_or("")
            .to_string()
    }

    fn read_number(&mut self) -> Result<f64> {
        let start = self.pos;
        // integer or decimal part
        while let Some(c) = self.peek() {
            if (c as char).is_ascii_digit() {
                self.pos += 1;
            } else {
                break;
            }
        }
        // fractional part
        if self.peek() == Some(b'.') {
            self.pos += 1;
            while let Some(c) = self.peek() {
                if (c as char).is_ascii_digit() {
                    self.pos += 1;
                } else {
                    break;
                }
            }
        }
        // exponent
        if matches!(self.peek(), Some(b'e') | Some(b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+') | Some(b'-')) {
                self.pos += 1;
            }
            while let Some(c) = self.peek() {
                if (c as char).is_ascii_digit() {
                    self.pos += 1;
                } else {
                    break;
                }
            }
        }
        let s = std::str::from_utf8(&self.src[start..self.pos]).unwrap_or("0");
        s.parse::<f64>().map_err(|e| {
            MathError::Parse(format!("could not parse number '{}': {}", s, e))
        })
    }

    // --- grammar ---------------------------------------------------------

    fn parse_expr(&mut self) -> Result<Expr> {
        let mut left = self.parse_term()?;
        loop {
            self.skip_ws();
            // Stop if we see \right (TeX closing delimiter)
            if self.peek() == Some(b'\\') {
                let saved = self.pos;
                self.pos += 1;
                let cmd = self.read_tex_command_name();
                self.pos = saved;
                if cmd == "right" {
                    break;
                }
            }
            match self.peek() {
                Some(b'+') => {
                    self.pos += 1;
                    let right = self.parse_term()?;
                    left = Expr::add(left, right);
                }
                Some(b'-') => {
                    self.pos += 1;
                    let right = self.parse_term()?;
                    left = Expr::sub(left, right);
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_term(&mut self) -> Result<Expr> {
        let mut left = self.parse_factor()?;
        loop {
            self.skip_ws();
            // Handle \cdot as explicit multiplication
            if self.peek() == Some(b'\\') {
                let saved = self.pos;
                self.pos += 1;
                let cmd = self.read_tex_command_name();
                if cmd == "cdot" || cmd == "times" {
                    let right = self.parse_factor()?;
                    left = Expr::mul(left, right);
                    continue;
                } else {
                    // Not a multiplication command — restore position
                    self.pos = saved;
                    // Stop on \right or other non-atom commands
                    if cmd == "right" {
                        break;
                    }
                    // Check if it's an atom starter for implicit multiplication
                    if self.is_atom_starter(Some(b'\\')) {
                        let right = self.parse_factor()?;
                        left = Expr::mul(left, right);
                        continue;
                    }
                    break;
                }
            }
            match self.peek() {
                Some(b'*') | Some(b'/') => {
                    let op = self.advance().unwrap();
                    let right = self.parse_factor()?;
                    left = if op == b'*' { Expr::mul(left, right) } else { Expr::div(left, right) };
                }
                // implicit multiplication: 2x, 2(x+1), x sin(y) (rare but allowed)
                Some(c) if self.is_atom_starter(Some(c)) => {
                    let right = self.parse_factor()?;
                    left = Expr::mul(left, right);
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_factor(&mut self) -> Result<Expr> {
        let base = self.parse_unary()?;
        self.skip_ws();
        if self.peek() == Some(b'^') {
            self.pos += 1;
            self.skip_ws();
            // Handle ^{...} (TeX brace group) or ^atom
            let exp = if self.peek() == Some(b'{') {
                self.parse_brace_group()?
            } else {
                self.parse_factor()? // right-associative
            };
            Ok(Expr::pow(base, exp))
        } else {
            Ok(base)
        }
    }

    fn parse_unary(&mut self) -> Result<Expr> {
        self.skip_ws();
        if self.peek() == Some(b'+') {
            self.pos += 1;
            return self.parse_unary();
        }
        if self.peek() == Some(b'-') {
            self.pos += 1;
            let e = self.parse_unary()?;
            return Ok(Expr::neg(e));
        }
        self.parse_atom()
    }

    fn is_atom_starter(&self, c: Option<u8>) -> bool {
        match c {
            // Parentheses can always be implicitly multiplied: 2(x+1)
            Some(b'(') => true,
            // Backslash starts a TeX command — can be implicitly multiplied
            Some(b'\\') => true,
            // Digits and identifiers can be implicitly multiplied: 2x, x sin(y)
            Some(c) => (c as char).is_ascii_digit() || self.starts_ident_c(c),
            // minus is NOT an atom starter here — it must be handled as a
            // binary subtraction at the expr level. Allowing implicit
            // multiplication on `-` makes `x^3 - 2*x` parse as `x^3 * -2 * x`.
            None => false,
        }
    }

    fn starts_ident_c(&self, c: u8) -> bool {
        matches!(c, b'a'..=b'z' | b'A'..=b'Z' | b'_')
    }

    fn parse_atom(&mut self) -> Result<Expr> {
        self.skip_ws();
        match self.peek() {
            Some(c) if (c as char).is_ascii_digit() => Ok(Expr::num(self.read_number()?)),
            Some(b'(') => {
                self.pos += 1;
                let e = self.parse_expr()?;
                self.expect(b')')?;
                Ok(e)
            }
            Some(b'\\') => self.parse_tex_command(),
            Some(b'{') => {
                self.pos += 1;
                let e = self.parse_expr()?;
                self.expect(b'}')?;
                Ok(e)
            }
            Some(c) if self.starts_ident_c(c) => {
                let name = self.read_ident();
                self.skip_ws();
                if self.peek() == Some(b'(') {
                    self.pos += 1;
                    let args = self.parse_arg_list()?;
                    self.expect(b')')?;
                    Ok(Expr::func(name, args))
                } else {
                    Ok(constant_or_var(&name))
                }
            }
            Some(c) => Err(MathError::Parse(format!(
                "unexpected character '{}' at position {}",
                c as char, self.pos
            ))),
            None => Err(MathError::Parse("unexpected end of input".into())),
        }
    }

    /// Parse a LaTeX command starting with `\`.
    fn parse_tex_command(&mut self) -> Result<Expr> {
        self.pos += 1; // consume '\'
        let cmd = self.read_tex_command_name();
        match cmd.as_str() {
            "frac" => {
                let num = self.parse_brace_group()?;
                let den = self.parse_brace_group()?;
                Ok(Expr::div(num, den))
            }
            "sqrt" => {
                let arg = self.parse_brace_group()?;
                Ok(Expr::func("sqrt", vec![arg]))
            }
            "pi" => Ok(Expr::num(std::f64::consts::PI)),
            "tau" => Ok(Expr::num(std::f64::consts::TAU)),
            "infty" => Ok(Expr::num(f64::INFINITY)),
            "cdot" | "times" => {
                // These are multiplication operators — shouldn't appear here
                // as atom starters, but if they do, treat as error
                Err(MathError::Parse(format!("unexpected \\{} at position {}", cmd, self.pos)))
            }
            "left" => {
                // \left( ... \right) — just parse the parenthesized group
                self.skip_ws();
                if self.peek() == Some(b'(') {
                    self.pos += 1;
                    let e = self.parse_expr()?;
                    self.skip_ws();
                    // Expect 
                    self.expect_tex_right()?;
                    Ok(e)
                } else if self.peek() == Some(b'[') {
                    self.pos += 1;
                    let e = self.parse_expr()?;
                    self.skip_ws();
                    self.expect_tex_right()?;
                    Ok(e)
                } else {
                    Err(MathError::Parse(format!("\\left must be followed by ( or [ at position {}", self.pos)))
                }
            }
            "right" => {
                Err(MathError::Parse(format!("unexpected \\right at position {}", self.pos)))
            }
            "log" => {
                // Could be \log, \log_2, \log_{10}
                self.skip_ws();
                if self.peek() == Some(b'_') {
                    self.pos += 1;
                    let base_str = self.read_subscript();
                    self.skip_ws();
                    let arg = if self.peek() == Some(b'{') {
                        self.parse_brace_group()?
                    } else if self.peek() == Some(b'(') {
                        self.pos += 1;
                        let args = self.parse_arg_list()?;
                        self.expect(b')')?;
                        if args.len() == 1 { args.into_iter().next().unwrap() }
                        else { return Ok(Expr::func("log", args)); }
                    } else {
                        self.parse_factor()?
                    };
                    // Try to parse base as f64, fall back to expression
                    if let Ok(base) = base_str.parse::<f64>() {
                        Ok(Expr::func("log", vec![arg, Expr::num(base)]))
                    } else {
                        let base_expr = Parser::parse(&base_str)?;
                        Ok(Expr::func("log", vec![arg, base_expr]))
                    }
                } else {
                    let arg = if self.peek() == Some(b'{') {
                        self.parse_brace_group()?
                    } else if self.peek() == Some(b'(') {
                        self.pos += 1;
                        let args = self.parse_arg_list()?;
                        self.expect(b')')?;
                        if args.len() == 1 { args.into_iter().next().unwrap() }
                        else { return Ok(Expr::func("log", args)); }
                    } else {
                        self.parse_factor()?
                    };
                    Ok(Expr::func("log", vec![arg]))
                }
            }
            "sin" | "cos" | "tan" | "asin" | "acos" | "atan"
            | "sinh" | "cosh" | "tanh" | "exp" | "ln"
            | "log2" | "log10" | "abs" | "floor" | "ceil" | "round"
            | "sign" | "cbrt" | "fract" | "gamma" | "erf" | "erfc" | "sinc" => {
                self.skip_ws();
                // Function may take argument in braces, parens, or just next atom
                let arg = if self.peek() == Some(b'{') {
                    self.parse_brace_group()?
                } else if self.peek() == Some(b'(') {
                    self.pos += 1;
                    let args = self.parse_arg_list()?;
                    self.expect(b')')?;
                    if args.len() == 1 { args.into_iter().next().unwrap() }
                    else { return Ok(Expr::func(cmd, args)); }
                } else {
                    self.parse_factor()?
                };
                Ok(Expr::func(cmd, vec![arg]))
            }
            "Gamma" => {
                self.skip_ws();
                let arg = if self.peek() == Some(b'{') {
                    self.parse_brace_group()?
                } else {
                    self.parse_factor()?
                };
                Ok(Expr::func("gamma", vec![arg]))
            }
            "operatorname" => {
                // \operatorname{erf}(x) etc.
                self.skip_ws();
                let name = self.parse_brace_group_str()?;
                self.skip_ws();
                if self.peek() == Some(b'(') {
                    self.pos += 1;
                    let args = self.parse_arg_list()?;
                    self.expect(b')')?;
                    Ok(Expr::func(name.to_lowercase(), args))
                } else if self.peek() == Some(b'{') {
                    let arg = self.parse_brace_group()?;
                    Ok(Expr::func(name.to_lowercase(), vec![arg]))
                } else {
                    let arg = self.parse_factor()?;
                    Ok(Expr::func(name.to_lowercase(), vec![arg]))
                }
            }
            "text" => {
                // \text{foo} — treat as variable name
                let name = self.parse_brace_group_str()?;
                Ok(constant_or_var(&name))
            }
            _ => {
                // Unknown command — try treating it as a variable name
                Ok(Expr::var(cmd))
            }
        }
    }

    /// Read a TeX command name after `\` (letters only, stops at non-letter).
    fn read_tex_command_name(&mut self) -> String {
        let start = self.pos;
        while let Some(c) = self.peek() {
            if (c as char).is_ascii_alphabetic() {
                self.pos += 1;
            } else {
                break;
            }
        }
        std::str::from_utf8(&self.src[start..self.pos])
            .unwrap_or("")
            .to_string()
    }

    /// Parse a `{...}` group and return the expression inside.
    fn parse_brace_group(&mut self) -> Result<Expr> {
        self.skip_ws();
        self.expect(b'{')?;
        let e = self.parse_expr()?;
        self.expect(b'}')?;
        Ok(e)
    }

    /// Parse a `{...}` group and return the raw string inside.
    fn parse_brace_group_str(&mut self) -> Result<String> {
        self.skip_ws();
        self.expect(b'{')?;
        let start = self.pos;
        let mut depth = 1;
        while let Some(c) = self.peek() {
            if c == b'{' { depth += 1; }
            else if c == b'}' { depth -= 1; if depth == 0 { break; } }
            self.pos += 1;
        }
        let s = std::str::from_utf8(&self.src[start..self.pos])
            .unwrap_or("")
            .trim()
            .to_string();
        self.expect(b'}')?;
        Ok(s)
    }

    /// Read a subscript (either a single char or a brace group), return as string.
    fn read_subscript(&mut self) -> String {
        if self.peek() == Some(b'{') {
            let start = self.pos + 1;
            self.pos += 1;
            let mut depth = 1;
            while let Some(c) = self.peek() {
                if c == b'{' { depth += 1; }
                else if c == b'}' { depth -= 1; if depth == 0 { break; } }
                self.pos += 1;
            }
            let s = std::str::from_utf8(&self.src[start..self.pos]).unwrap_or("").to_string();
            let _ = self.advance(); // consume }
            s
        } else if let Some(c) = self.advance() {
            (c as char).to_string()
        } else {
            String::new()
        }
    }

    /// Expect `\right)` or `\right]`
    fn expect_tex_right(&mut self) -> Result<()> {
        self.skip_ws();
        if self.peek() == Some(b'\\') {
            self.pos += 1;
            let cmd = self.read_tex_command_name();
            if cmd == "right" {
                self.skip_ws();
                // Consume the closing bracket (either ) or ])
                if matches!(self.peek(), Some(b')') | Some(b']')) {
                    self.pos += 1;
                    return Ok(());
                }
            }
        }
        Err(MathError::Parse(format!("expected \\right) at position {}", self.pos)))
    }

    fn parse_arg_list(&mut self) -> Result<Vec<Expr>> {
        let mut args = Vec::new();
        self.skip_ws();
        if self.peek() == Some(b')') {
            return Ok(args);
        }
        args.push(self.parse_expr()?);
        loop {
            self.skip_ws();
            if self.peek() != Some(b',') {
                break;
            }
            self.pos += 1;
            args.push(self.parse_expr()?);
        }
        Ok(args)
    }
}

/// Strip Markdown/LaTeX math delimiters from a string.
/// Handles: `$$...$$`, `$...$`, `\[...\]`, `\(...\)`.
/// Returns a borrowed slice when no delimiters are present (avoids allocation).
fn strip_math_delimiters(src: &str) -> Cow<'_, str> {
    let trimmed = src.trim();
    // $$...$$ (display math)
    if trimmed.starts_with("$$") && trimmed.ends_with("$$") && trimmed.len() > 4 {
        return Cow::Owned(trimmed[2..trimmed.len() - 2].trim().to_string());
    }
    // \[...\] (display math)
    if trimmed.starts_with("\\[") && trimmed.ends_with("\\]") && trimmed.len() > 4 {
        return Cow::Owned(trimmed[2..trimmed.len() - 2].trim().to_string());
    }
    // \(...\) (inline math)
    if trimmed.starts_with("\\(") && trimmed.ends_with("\\)") && trimmed.len() > 4 {
        return Cow::Owned(trimmed[2..trimmed.len() - 2].trim().to_string());
    }
    // $...$ (inline math) — single $ on each side, not $$
    if trimmed.starts_with('$') && trimmed.ends_with('$') && trimmed.len() > 2 {
        return Cow::Owned(trimmed[1..trimmed.len() - 1].trim().to_string());
    }
    Cow::Borrowed(src)
}

fn constant_or_var(name: &str) -> Expr {
    match name {
        "pi" | "PI" | "Pi" => Expr::num(std::f64::consts::PI),
        "e" => Expr::num(std::f64::consts::E),
        "tau" => Expr::num(std::f64::consts::TAU),
        "inf" | "Inf" | "Infinity" => Expr::num(f64::INFINITY),
        "nan" | "NaN" => Expr::num(f64::NAN),
        _ => Expr::var(name),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> Expr {
        Parser::parse(s).unwrap()
    }

    #[test]
    fn numbers_and_vars() {
        assert_eq!(p("3"), Expr::num(3.0));
        assert_eq!(p("3.14"), Expr::num(3.14));
        assert_eq!(p("x"), Expr::var("x"));
        assert_eq!(p("pi"), Expr::num(std::f64::consts::PI));
    }

    #[test]
    fn precedence() {
        assert_eq!(p("1 + 2*3"), Expr::add(Expr::num(1.0), Expr::mul(Expr::num(2.0), Expr::num(3.0))));
        assert_eq!(p("2^3^2"), Expr::pow(Expr::num(2.0), Expr::pow(Expr::num(3.0), Expr::num(2.0))));
    }

    #[test]
    fn unary_and_parens() {
        assert_eq!(p("-x"), Expr::neg(Expr::var("x")));
        assert_eq!(p("-(1+2)"), Expr::neg(Expr::add(Expr::num(1.0), Expr::num(2.0))));
    }

    #[test]
    fn functions() {
        assert_eq!(p("sin(x)"), Expr::func("sin", vec![Expr::var("x")]));
        assert_eq!(
            p("log(2, 8)"),
            Expr::func("log", vec![Expr::num(2.0), Expr::num(8.0)])
        );
    }

    #[test]
    fn implicit_multiply() {
        assert_eq!(p("2x"), Expr::mul(Expr::num(2.0), Expr::var("x")));
        assert_eq!(
            p("3(x+1)"),
            Expr::mul(Expr::num(3.0), Expr::add(Expr::var("x"), Expr::num(1.0)))
        );
        assert_eq!(
            p("(x)(y)"),
            Expr::mul(Expr::var("x"), Expr::var("y"))
        );
    }

    #[test]
    fn scientific_notation() {
        assert_eq!(p("1.5e3"), Expr::num(1500.0));
        assert_eq!(p("2E-2"), Expr::num(0.02));
    }

    #[test]
    fn complex_expression() {
        let e = Parser::parse("sin(x)^2 + cos(x)^2").unwrap();
        // Just ensure it parses and round-trips.
        let _ = format!("{}", e);
    }

    #[test]
    fn tex_frac() {
        assert_eq!(
            p(r"\frac{1}{2}"),
            Expr::div(Expr::num(1.0), Expr::num(2.0))
        );
        assert_eq!(
            p(r"\frac{x+1}{x-1}"),
            Expr::div(
                Expr::add(Expr::var("x"), Expr::num(1.0)),
                Expr::sub(Expr::var("x"), Expr::num(1.0))
            )
        );
    }

    #[test]
    fn tex_sqrt() {
        assert_eq!(
            p(r"\sqrt{16}"),
            Expr::func("sqrt", vec![Expr::num(16.0)])
        );
    }

    #[test]
    fn tex_pi() {
        assert_eq!(p(r"\pi"), Expr::num(std::f64::consts::PI));
    }

    #[test]
    fn tex_cdot() {
        assert_eq!(
            p(r"2 \cdot 3"),
            Expr::mul(Expr::num(2.0), Expr::num(3.0))
        );
    }

    #[test]
    fn tex_left_right() {
        assert_eq!(
            p(r"\left( 1 + 2 \right) \cdot 3"),
            Expr::mul(
                Expr::add(Expr::num(1.0), Expr::num(2.0)),
                Expr::num(3.0)
            )
        );
    }

    #[test]
    fn tex_pow_braces() {
        assert_eq!(
            p(r"x^{2}"),
            Expr::pow(Expr::var("x"), Expr::num(2.0))
        );
        assert_eq!(
            p(r"x^{2+y}"),
            Expr::pow(Expr::var("x"), Expr::add(Expr::num(2.0), Expr::var("y")))
        );
    }

    #[test]
    fn tex_sin_brace() {
        assert_eq!(
            p(r"\sin{x}"),
            Expr::func("sin", vec![Expr::var("x")])
        );
        assert_eq!(
            p(r"\sin(x)"),
            Expr::func("sin", vec![Expr::var("x")])
        );
    }

    #[test]
    fn tex_log_subscript() {
        assert_eq!(
            p(r"\log_2{8}"),
            Expr::func("log", vec![Expr::num(8.0), Expr::num(2.0)])
        );
    }

    #[test]
    fn tex_gamma() {
        assert_eq!(
            p(r"\Gamma{0.5}"),
            Expr::func("gamma", vec![Expr::num(0.5)])
        );
    }

    #[test]
    fn tex_operatorname() {
        assert_eq!(
            p(r"\operatorname{erf}(1.0)"),
            Expr::func("erf", vec![Expr::num(1.0)])
        );
    }

    #[test]
    fn tex_input_roundtrip() {
        // Parse TeX input, verify it matches plain-text equivalent
        assert_eq!(
            Parser::parse(r"\frac{x+1}{x-1}").unwrap(),
            Parser::parse("(x+1)/(x-1)").unwrap()
        );
        assert_eq!(
            Parser::parse(r"\sin(x^2) + \cos(2 \cdot x)").unwrap(),
            Parser::parse("sin(x^2) + cos(2*x)").unwrap()
        );
        assert_eq!(
            Parser::parse(r"\sqrt{x} + \Gamma(0.5)").unwrap(),
            Parser::parse("sqrt(x) + gamma(0.5)").unwrap()
        );
    }

    #[test]
    fn markdown_inline_math() {
        assert_eq!(
            Parser::parse(r"$\frac{1}{2} + \frac{3}{4}$").unwrap(),
            Parser::parse(r"\frac{1}{2} + \frac{3}{4}").unwrap()
        );
        assert_eq!(
            Parser::parse("$\\sin(\\pi / 4)$").unwrap(),
            Parser::parse("sin(pi/4)").unwrap()
        );
    }

    #[test]
    fn markdown_display_math() {
        assert_eq!(
            Parser::parse("$$\\frac{x^2 - 1}{x - 1}$$").unwrap(),
            Parser::parse("(x^2-1)/(x-1)").unwrap()
        );
    }

    #[test]
    fn latex_bracket_delimiters() {
        assert_eq!(
            Parser::parse(r"\[\sqrt{16}\]").unwrap(),
            Parser::parse("sqrt(16)").unwrap()
        );
        assert_eq!(
            Parser::parse(r"\(\log_2{8}\)").unwrap(),
            Parser::parse("log(8, 2)").unwrap()
        );
    }

    #[test]
    fn tex_nested_frac() {
        assert_eq!(
            Parser::parse(r"\frac{\frac{1}{2}}{\frac{3}{4}}").unwrap(),
            Parser::parse("(1/2)/(3/4)").unwrap()
        );
    }

    #[test]
    fn tex_implicit_multiplication_with_tex_commands() {
        assert_eq!(
            Parser::parse(r"2\pi").unwrap(),
            Expr::mul(Expr::num(2.0), Expr::num(std::f64::consts::PI))
        );
        assert_eq!(
            Parser::parse(r"3\sqrt{16}").unwrap(),
            Expr::mul(Expr::num(3.0), Expr::func("sqrt", vec![Expr::num(16.0)]))
        );
    }

    #[test]
    fn tex_mixed_plain_and_tex() {
        assert_eq!(
            Parser::parse(r"x + \frac{1}{2}").unwrap(),
            Parser::parse("x + 1/2").unwrap()
        );
        assert_eq!(
            Parser::parse(r"\sin(x) * \cos(x)").unwrap(),
            Parser::parse("sin(x) * cos(x)").unwrap()
        );
    }

    #[test]
    fn tex_chained_operations() {
        assert_eq!(
            Parser::parse(r"\frac{1}{2} + \frac{1}{3} + \frac{1}{6}").unwrap(),
            Parser::parse("1/2 + 1/3 + 1/6").unwrap()
        );
    }

    #[test]
    fn tex_pow_with_expression() {
        assert_eq!(
            Parser::parse(r"x^{2y}").unwrap(),
            Expr::pow(Expr::var("x"), Expr::mul(Expr::num(2.0), Expr::var("y")))
        );
        assert_eq!(
            Parser::parse(r"2^{x+1}").unwrap(),
            Expr::pow(Expr::num(2.0), Expr::add(Expr::var("x"), Expr::num(1.0)))
        );
    }

    #[test]
    fn tex_left_right_nested() {
        assert_eq!(
            Parser::parse(r"\left( \left( 1 + 2 \right) + 3 \right)").unwrap(),
            Parser::parse("((1+2)+3)").unwrap()
        );
    }

    #[test]
    fn tex_text_command() {
        assert_eq!(
            Parser::parse(r"\text{alpha} + 1").unwrap(),
            Expr::add(Expr::var("alpha"), Expr::num(1.0))
        );
    }

    #[test]
    fn markdown_mixed_delimiters() {
        // All delimiter styles should produce same result
        let plain = Parser::parse(r"\frac{1}{2}").unwrap();
        assert_eq!(Parser::parse(r"$\frac{1}{2}$").unwrap(), plain);
        assert_eq!(Parser::parse(r"$$\frac{1}{2}$$").unwrap(), plain);
        assert_eq!(Parser::parse(r"\(\frac{1}{2}\)").unwrap(), plain);
        assert_eq!(Parser::parse(r"\[\frac{1}{2}\]").unwrap(), plain);
    }

    #[test]
    fn no_delimiter_unchanged() {
        // Plain expressions without delimiters should work as before
        assert_eq!(
            Parser::parse("sin(pi/4) + 2^3").unwrap(),
            Parser::parse("sin(pi/4) + 2^3").unwrap()
        );
    }
}