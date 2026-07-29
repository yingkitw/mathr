use std::fmt;

/// A symbolic math expression used by the parser, evaluator, differentiator,
/// simplifier and solver. Numeric coefficients are stored as `f64`; symbolic
/// leaves are stored as variable names.
#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Num(f64),
    Var(String),
    Neg(Box<Expr>),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    /// Function call: name and arguments.
    Func(String, Vec<Expr>),
}

pub use Expr::*;

impl Expr {
    pub fn num(x: f64) -> Expr {
        Num(x)
    }

    pub fn var<S: Into<String>>(name: S) -> Expr {
        Var(name.into())
    }

    pub fn neg(e: Expr) -> Expr {
        Neg(Box::new(e))
    }

    pub fn add(a: Expr, b: Expr) -> Expr {
        Add(Box::new(a), Box::new(b))
    }

    pub fn sub(a: Expr, b: Expr) -> Expr {
        Sub(Box::new(a), Box::new(b))
    }

    pub fn mul(a: Expr, b: Expr) -> Expr {
        Mul(Box::new(a), Box::new(b))
    }

    pub fn div(a: Expr, b: Expr) -> Expr {
        Div(Box::new(a), Box::new(b))
    }

    pub fn pow(a: Expr, b: Expr) -> Expr {
        Pow(Box::new(a), Box::new(b))
    }

    pub fn func<S: Into<String>>(name: S, args: Vec<Expr>) -> Expr {
        Func(name.into(), args)
    }

    /// Collect all free variable names in the expression.
    pub fn variables(&self) -> Vec<String> {
        let mut out = Vec::new();
        self.collect_vars(&mut out);
        let mut seen = std::collections::HashSet::new();
        out.retain(|v| seen.insert(v.clone()));
        out
    }

    fn collect_vars(&self, out: &mut Vec<String>) {
        match self {
            Num(_) => {}
            Var(v) => out.push(v.clone()),
            Neg(e) => e.collect_vars(out),
            Add(a, b) | Sub(a, b) | Mul(a, b) | Div(a, b) | Pow(a, b) => {
                a.collect_vars(out);
                b.collect_vars(out);
            }
            Func(_, args) => args.iter().for_each(|a| a.collect_vars(out)),
        }
    }

    /// Substitute every occurrence of `var` with `replacement`.
    pub fn substitute(&self, var: &str, replacement: &Expr) -> Expr {
        match self {
            Num(_) => self.clone(),
            Var(v) if v == var => replacement.clone(),
            Var(_) => self.clone(),
            Neg(e) => Expr::neg(e.substitute(var, replacement)),
            Add(a, b) => Expr::add(a.substitute(var, replacement), b.substitute(var, replacement)),
            Sub(a, b) => Expr::sub(a.substitute(var, replacement), b.substitute(var, replacement)),
            Mul(a, b) => Expr::mul(a.substitute(var, replacement), b.substitute(var, replacement)),
            Div(a, b) => Expr::div(a.substitute(var, replacement), b.substitute(var, replacement)),
            Pow(a, b) => Expr::pow(a.substitute(var, replacement), b.substitute(var, replacement)),
            Func(name, args) => Expr::func(
                name.clone(),
                args.iter().map(|a| a.substitute(var, replacement)).collect(),
            ),
        }
    }

    /// Check whether this is a constant (no variables).
    pub fn is_constant(&self) -> bool {
        match self {
            Num(_) => true,
            Var(_) => false,
            Neg(e) => e.is_constant(),
            Add(a, b) | Sub(a, b) | Mul(a, b) | Div(a, b) | Pow(a, b) => {
                a.is_constant() && b.is_constant()
            }
            Func(_, args) => args.iter().all(|a| a.is_constant()),
        }
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write(f, 0, false)
    }
}

impl Expr {
    /// Pretty-printer with operator precedence and parenthesisation.
    /// `prec` is the surrounding context's precedence (0 = top-level).
    /// `in_arg` is true when this expression is itself a function argument —
    /// multiplies get parenthesised to avoid `f x y` ambiguity.
    fn write(
        &self,
        f: &mut fmt::Formatter<'_>,
        parent_prec: u8,
        in_arg: bool,
    ) -> fmt::Result {
        let (my_prec, parens) = match self {
            Num(_) | Var(_) | Func(_, _) => (10, false),
            Neg(_) => (8, true),
            Pow(_, _) => (7, false),
            Mul(_, _) | Div(_, _) => (5, true),
            Add(_, _) | Sub(_, _) => (3, true),
        };
        let need_parens = parens && my_prec < parent_prec
            || in_arg && matches!(self, Mul(_, _) | Div(_, _) | Neg(_));
        if need_parens {
            f.write_str("(")?;
        }
        match self {
            Num(n) => {
                if n.is_nan() {
                    write!(f, "NaN")?;
                } else if n.is_infinite() {
                    if *n > 0.0 {
                        write!(f, "inf")?;
                    } else {
                        write!(f, "-inf")?;
                    }
                } else if n.fract() == 0.0 && n.abs() < 1e15 {
                    write!(f, "{}", *n as i64)?;
                } else {
                    write!(f, "{}", n)?;
                }
            }
            Var(v) => f.write_str(v)?,
            Neg(e) => {
                f.write_str("-")?;
                e.write(f, my_prec, in_arg)?;
            }
            Add(a, b) => {
                a.write(f, my_prec, false)?;
                f.write_str(" + ")?;
                b.write(f, my_prec, false)?;
            }
            Sub(a, b) => {
                a.write(f, my_prec, false)?;
                f.write_str(" - ")?;
                b.write(f, my_prec, false)?;
            }
            Mul(a, b) => {
                a.write(f, my_prec, false)?;
                f.write_str("*")?;
                b.write(f, my_prec, false)?;
            }
            Div(a, b) => {
                a.write(f, my_prec, false)?;
                f.write_str("/")?;
                b.write(f, my_prec, false)?;
            }
            Pow(a, b) => {
                a.write(f, my_prec + 1, false)?;
                f.write_str("^")?;
                b.write(f, my_prec + 1, false)?;
            }
            Func(name, args) => {
                write!(f, "{}(", name)?;
                for (i, a) in args.iter().enumerate() {
                    if i > 0 {
                        f.write_str(", ")?;
                    }
                    a.write(f, 0, true)?;
                }
                f.write_str(")")?;
            }
        }
        if need_parens {
            f.write_str(")")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variables_collected() {
        let e = Expr::add(Expr::var("x"), Expr::mul(Expr::var("y"), Expr::var("x")));
        assert_eq!(e.variables(), vec!["x".to_string(), "y".to_string()]);
    }

    #[test]
    fn substitution() {
        let e = Expr::mul(Expr::var("x"), Expr::var("y"));
        let sub = e.substitute("x", &Expr::num(2.0));
        assert_eq!(sub, Expr::mul(Expr::num(2.0), Expr::var("y")));
    }

    #[test]
    fn display_polishes_parens() {
        let e = Expr::add(
            Expr::num(1.0),
            Expr::mul(Expr::num(2.0), Expr::var("x")),
        );
        assert_eq!(format!("{}", e), "1 + 2*x");
    }

    #[test]
    fn display_function_call() {
        let e = Expr::func("sin", vec![Expr::var("x")]);
        assert_eq!(format!("{}", e), "sin(x)");
    }
}