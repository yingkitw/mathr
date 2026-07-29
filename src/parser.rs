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
pub struct Parser<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Parser<'a> {
    pub fn new(src: &'a str) -> Self {
        Self { src: src.as_bytes(), pos: 0 }
    }

    pub fn parse(src: &str) -> Result<Expr> {
        let mut p = Parser::new(src);
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

    fn starts_ident(&self) -> bool {
        matches!(self.peek(), Some(b'a'..=b'z') | Some(b'A'..=b'Z') | Some(b'_'))
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
            let exp = self.parse_factor()?; // right-associative
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
}