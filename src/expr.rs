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

    /// Produce a canonical form of the expression for equality checking.
    ///
    /// Normalizes:
    /// - `Sub(a, b)` → `Add(a, Neg(b))`
    /// - Flattens nested `Add`/`Mul` into sorted n-ary lists
    /// - Sorts commutative operands by their string representation
    /// - Folds numeric constants in Add/Mul
    pub fn canonicalize(&self) -> Expr {
        match self {
            Num(n) => Num(*n),
            Var(v) => Var(v.clone()),
            Neg(e) => Expr::neg(e.canonicalize()),
            Add(_, _) => {
                let terms = flatten_add(self);
                let canonical_terms: Vec<Expr> = terms.iter().map(|t| t.canonicalize()).collect();
                recombine_add(&canonical_terms)
            }
            Sub(a, b) => {
                let terms = flatten_add(&Expr::sub((**a).clone(), (**b).clone()));
                let canonical_terms: Vec<Expr> = terms.iter().map(|t| t.canonicalize()).collect();
                recombine_add(&canonical_terms)
            }
            Mul(_, _) => {
                let factors = flatten_mul(self);
                let canonical_factors: Vec<Expr> = factors.iter().map(|f| f.canonicalize()).collect();
                recombine_mul(&canonical_factors)
            }
            Div(a, b) => Expr::div(a.canonicalize(), b.canonicalize()),
            Pow(a, b) => Expr::pow(a.canonicalize(), b.canonicalize()),
            Func(name, args) => {
                Expr::func(name.clone(), args.iter().map(|a| a.canonicalize()).collect())
            }
        }
    }

    /// Check structural equality after canonicalization.
    /// Two expressions are equal if they have the same canonical form.
    pub fn equals(&self, other: &Expr) -> bool {
        self.canonicalize() == other.canonicalize()
    }
}

/// Flatten an expression into a list of additive terms.
/// `a + b - c` → `[a, b, Neg(c)]`
fn flatten_add(e: &Expr) -> Vec<Expr> {
    match e {
        Add(a, b) => {
            let mut v = flatten_add(a);
            v.extend(flatten_add(b));
            v
        }
        Sub(a, b) => {
            let mut v = flatten_add(a);
            v.extend(flatten_add(b).into_iter().map(Expr::neg));
            v
        }
        Neg(e) => flatten_add(e).into_iter().map(Expr::neg).collect(),
        _ => vec![e.clone()],
    }
}

/// Flatten an expression into a list of multiplicative factors.
/// `a * b` → `[a, b]`
fn flatten_mul(e: &Expr) -> Vec<Expr> {
    match e {
        Mul(a, b) => {
            let mut v = flatten_mul(a);
            v.extend(flatten_mul(b));
            v
        }
        _ => vec![e.clone()],
    }
}

/// Recombine a list of terms into a canonical Add expression.
/// Sorts terms by their string representation, folds numeric constants.
fn recombine_add(terms: &[Expr]) -> Expr {
    if terms.is_empty() {
        return Expr::num(0.0);
    }
    if terms.len() == 1 {
        return terms[0].clone();
    }

    // Separate numeric and non-numeric terms
    let mut num_sum = 0.0;
    let mut symbolic: Vec<Expr> = Vec::new();
    for t in terms {
        match t {
            Num(n) => num_sum += n,
            Neg(inner) => {
                if let Num(n) = **inner {
                    num_sum -= n;
                } else {
                    symbolic.push(t.clone());
                }
            }
            _ => symbolic.push(t.clone()),
        }
    }

    // Sort symbolic terms by cached string representation
    let mut keyed: Vec<(String, Expr)> = symbolic
        .into_iter()
        .map(|e| (e.to_string(), e))
        .collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0));

    let mut result = if num_sum != 0.0 {
        Expr::num(num_sum)
    } else if !keyed.is_empty() {
        keyed.remove(0).1
    } else {
        return Expr::num(0.0);
    };

    for (_, t) in keyed {
        result = Expr::add(result, t);
    }
    result
}

/// Recombine a list of factors into a canonical Mul expression.
/// Sorts factors by their string representation, folds numeric constants.
fn recombine_mul(factors: &[Expr]) -> Expr {
    if factors.is_empty() {
        return Expr::num(1.0);
    }
    if factors.len() == 1 {
        return factors[0].clone();
    }

    let mut num_prod = 1.0;
    let mut symbolic: Vec<Expr> = Vec::new();
    for f in factors {
        match f {
            Num(n) => num_prod *= n,
            _ => symbolic.push(f.clone()),
        }
    }

    let mut keyed: Vec<(String, Expr)> = symbolic
        .into_iter()
        .map(|e| (e.to_string(), e))
        .collect();
    keyed.sort_by(|a, b| a.0.cmp(&b.0));

    if num_prod == 0.0 {
        return Expr::num(0.0);
    }

    let mut result = if num_prod != 1.0 {
        Expr::num(num_prod)
    } else if !keyed.is_empty() {
        keyed.remove(0).1
    } else {
        return Expr::num(1.0);
    };

    for (_, f) in keyed {
        result = Expr::mul(result, f);
    }
    result
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.write(f, 0, false)
    }
}

/// Map of function names → LaTeX macros.
fn tex_fn_name(name: &str) -> Option<&'static str> {
    match name {
        "sin" => Some(r"\sin"),
        "cos" => Some(r"\cos"),
        "tan" => Some(r"\tan"),
        "asin" => Some(r"\arcsin"),
        "acos" => Some(r"\arccos"),
        "atan" => Some(r"\arctan"),
        "sinh" => Some(r"\sinh"),
        "cosh" => Some(r"\cosh"),
        "tanh" => Some(r"\tanh"),
        "exp" => Some(r"\exp"),
        "ln" => Some(r"\ln"),
        "log" => Some(r"\log"),
        "log2" => Some(r"\log_2"),
        "log10" => Some(r"\log_{10}"),
        "sqrt" => Some(r"\sqrt"),
        "abs" => Some(r"\abs"),
        "floor" => Some(r"\lfloor"),
        "ceil" => Some(r"\lceil"),
        "gamma" => Some(r"\Gamma"),
        "erf" => Some(r"\operatorname{erf}"),
        "erfc" => Some(r"\operatorname{erfc}"),
        "sinc" => Some(r"\operatorname{sinc}"),
        _ => None,
    }
}

impl Expr {
    /// Render this expression as a LaTeX math string.
    pub fn to_tex(&self) -> String {
        let mut s = String::new();
        self.write_tex(&mut s, 0);
        s
    }

    fn write_tex(&self, out: &mut String, parent_prec: u8) {
        let (my_prec, parens) = match self {
            Num(_) | Var(_) | Func(_, _) => (10, false),
            Neg(_) => (8, true),
            Pow(_, _) => (7, false),
            Mul(_, _) | Div(_, _) => (5, true),
            Add(_, _) | Sub(_, _) => (3, true),
        };
        let need_parens = parens && my_prec < parent_prec;

        if need_parens {
            out.push_str(r"\left(");
        }
        match self {
            Num(n) => {
                if n.is_nan() {
                    out.push_str(r"\text{NaN}");
                } else if n.is_infinite() {
                    out.push_str(if *n > 0.0 { r"\infty" } else { r"-\infty" });
                } else if n.fract() == 0.0 && n.abs() < 1e15 {
                    out.push_str(&format!("{}", *n as i64));
                } else {
                    out.push_str(&format!("{}", n));
                }
            }
            Var(v) => {
                match v.as_str() {
                    "pi" => out.push_str(r"\pi"),
                    "tau" => out.push_str(r"\tau"),
                    "inf" => out.push_str(r"\infty"),
                    _ => {
                        // Multi-char variables get subscript treatment for Greek-like names
                        if v.len() > 1 && v.chars().all(|c| c.is_ascii_alphabetic()) {
                            out.push_str(&format!(r"\text{{{}}}", v));
                        } else {
                            out.push_str(v);
                        }
                    }
                }
            }
            Neg(e) => {
                out.push('-');
                e.write_tex(out, my_prec);
            }
            Add(a, b) => {
                a.write_tex(out, my_prec);
                out.push_str(" + ");
                b.write_tex(out, my_prec);
            }
            Sub(a, b) => {
                a.write_tex(out, my_prec);
                out.push_str(" - ");
                b.write_tex(out, my_prec);
            }
            Mul(a, b) => {
                a.write_tex(out, my_prec);
                // Use \cdot for numeric * symbolic, juxtaposition for symbolic * symbolic
                let a_is_num = matches!(**a, Num(_));
                let b_is_num = matches!(**b, Num(_));
                if a_is_num || b_is_num {
                    out.push_str(r" \cdot ");
                } else {
                    out.push(' ');
                }
                b.write_tex(out, my_prec);
            }
            Div(a, b) => {
                // Use \frac for division
                out.push_str(r"\frac{");
                a.write_tex(out, 0);
                out.push_str("}{");
                b.write_tex(out, 0);
                out.push('}');
            }
            Pow(a, b) => {
                a.write_tex(out, my_prec + 1);
                out.push_str("^{");
                b.write_tex(out, 0);
                out.push('}');
            }
            Func(name, args) => {
                match tex_fn_name(name) {
                    Some(r"\sqrt") => {
                        out.push_str(r"\sqrt{");
                        if let Some(a) = args.first() {
                            a.write_tex(out, 0);
                        }
                        out.push('}');
                    }
                    Some(r"\lfloor") => {
                        out.push_str(r"\lfloor ");
                        if let Some(a) = args.first() { a.write_tex(out, 0); }
                        out.push_str(r" \rfloor");
                    }
                    Some(r"\lceil") => {
                        out.push_str(r"\lceil ");
                        if let Some(a) = args.first() { a.write_tex(out, 0); }
                        out.push_str(r" \rceil");
                    }
                    Some(tex_name) => {
                        out.push_str(tex_name);
                        if !args.is_empty() {
                            out.push_str(r"\left(");
                            for (i, a) in args.iter().enumerate() {
                                if i > 0 { out.push_str(", "); }
                                a.write_tex(out, 0);
                            }
                            out.push_str(r"\right)");
                        }
                    }
                    None => {
                        out.push_str(r"\operatorname{");
                        out.push_str(name);
                        out.push_str(r"}\left(");
                        for (i, a) in args.iter().enumerate() {
                            if i > 0 { out.push_str(", "); }
                            a.write_tex(out, 0);
                        }
                        out.push_str(r"\right)");
                    }
                }
            }
        }
        if need_parens {
            out.push_str(r"\right)");
        }
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

    #[test]
    fn equals_add_commutative() {
        // x + y == y + x
        let a = Expr::add(Expr::var("x"), Expr::var("y"));
        let b = Expr::add(Expr::var("y"), Expr::var("x"));
        assert!(a.equals(&b));
    }

    #[test]
    fn equals_mul_commutative() {
        // x * y == y * x
        let a = Expr::mul(Expr::var("x"), Expr::var("y"));
        let b = Expr::mul(Expr::var("y"), Expr::var("x"));
        assert!(a.equals(&b));
    }

    #[test]
    fn equals_sub_vs_add_neg() {
        // a - b == a + (-b)
        let a = Expr::sub(Expr::var("x"), Expr::var("y"));
        let b = Expr::add(Expr::var("x"), Expr::neg(Expr::var("y")));
        assert!(a.equals(&b));
    }

    #[test]
    fn equals_flatten_and_sort() {
        // (a + b) + c == c + (a + b)
        let a = Expr::add(Expr::add(Expr::var("a"), Expr::var("b")), Expr::var("c"));
        let b = Expr::add(Expr::var("c"), Expr::add(Expr::var("a"), Expr::var("b")));
        assert!(a.equals(&b));
    }

    #[test]
    fn equals_not_equal() {
        let a = Expr::add(Expr::var("x"), Expr::num(1.0));
        let b = Expr::add(Expr::var("x"), Expr::num(2.0));
        assert!(!a.equals(&b));
    }

    #[test]
    fn equals_constant_fold() {
        // 2 + 3 + x == 5 + x
        let a = Expr::add(Expr::add(Expr::num(2.0), Expr::num(3.0)), Expr::var("x"));
        let b = Expr::add(Expr::num(5.0), Expr::var("x"));
        assert!(a.equals(&b));
    }

    #[test]
    fn equals_mul_constant_fold() {
        // 2 * 3 * x == 6 * x
        let a = Expr::mul(Expr::mul(Expr::num(2.0), Expr::num(3.0)), Expr::var("x"));
        let b = Expr::mul(Expr::num(6.0), Expr::var("x"));
        assert!(a.equals(&b));
    }
}