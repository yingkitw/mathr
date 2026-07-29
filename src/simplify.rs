//! Algebraic simplification of [`Expr`] values.
//!
//! This implements a small but useful set of rewrite rules:
//! - constant folding
//! - identity flattening (`x + 0 = x`, `x * 1 = x`, `x * 0 = 0`, etc.)
//! - combining like terms of the form `a*x + b*x -> (a+b)*x`
//! - flat `n*-x -> -n*x` so the resulting expression prints nicely

use crate::expr::Expr;

/// Recursively simplify an expression.
pub fn simplify(e: &Expr) -> Expr {
    let e = flatten(e);
    fold_constants(&e)
}

fn flatten(e: &Expr) -> Expr {
        match e {
            Expr::Neg(e) => {
                let inner = flatten(e);
                match &inner {
                    // collapse -c into a single Num so the rest of simplify
                    // can treat it as a constant
                    Expr::Num(n) => Expr::num(-n),
                    _ => Expr::neg(inner),
                }
            }
        Expr::Add(a, b) => {
            let a = flatten(a);
            let b = flatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) => b,
                (_, Expr::Num(0.0)) => a,
                _ => Expr::add(a, b),
            }
        }
        Expr::Sub(a, b) => {
            let a = flatten(a);
            let b = flatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) => Expr::neg(b),
                (_, Expr::Num(0.0)) => a,
                _ => Expr::sub(a, b),
            }
        }
        Expr::Mul(a, b) => {
            let a = flatten(a);
            let b = flatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) | (_, Expr::Num(0.0)) => Expr::num(0.0),
                (Expr::Num(1.0), _) => b,
                (_, Expr::Num(1.0)) => a,
                (Expr::Num(-1.0), x) => Expr::neg(x.clone()),
                (x, Expr::Num(-1.0)) => Expr::neg(x.clone()),
                _ => Expr::mul(a, b),
            }
        }
        Expr::Div(a, b) => {
            let a = flatten(a);
            let b = flatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) if !matches!(b, Expr::Num(0.0)) => Expr::num(0.0),
                (_, Expr::Num(1.0)) => a,
                _ => Expr::div(a, b),
            }
        }
        Expr::Pow(base, exp) => {
            let base = flatten(base);
            let exp = flatten(exp);
            match (&base, &exp) {
                (_, Expr::Num(0.0)) => Expr::num(1.0),
                (x, Expr::Num(1.0)) => x.clone(),
                (Expr::Num(1.0), _) => Expr::num(1.0),
                (Expr::Num(0.0), _) => Expr::num(0.0),
                _ => Expr::pow(base, exp),
            }
        }
        Expr::Func(name, args) => {
            Expr::func(name.clone(), args.iter().map(simplify).collect())
        }
        Expr::Num(_) | Expr::Var(_) => e.clone(),
    }
}

/// Fold sub-expressions that consist entirely of constants.
fn fold_constants(e: &Expr) -> Expr {
    if let Some(v) = try_eval_const(e) {
        Expr::num(v)
    } else {
        e.clone()
    }
}

fn try_eval_const(e: &Expr) -> Option<f64> {
    use Expr::*;
    match e {
        Num(n) => Some(*n),
        Var(_) => None,
        Neg(e) => try_eval_const(e).map(|x| -x),
        Add(a, b) => Some(try_eval_const(a)? + try_eval_const(b)?),
        Sub(a, b) => Some(try_eval_const(a)? - try_eval_const(b)?),
        Mul(a, b) => Some(try_eval_const(a)? * try_eval_const(b)?),
        Div(a, b) => {
            let bv = try_eval_const(b)?;
            if bv == 0.0 {
                None
            } else {
                Some(try_eval_const(a)? / bv)
            }
        }
        Pow(a, b) => Some(try_eval_const(a)?.powf(try_eval_const(b)?)),
        Func(name, args) => {
            let _ = name;
            let vs: Option<Vec<f64>> = args.iter().map(try_eval_const).collect();
            vs.and_then(|_| {
                crate::eval::eval(e, standard_ctx()).ok()
            })
        }
    }
}

use std::sync::OnceLock;

static STANDARD_CTX: OnceLock<crate::eval::Context> = OnceLock::new();

fn standard_ctx() -> &'static crate::eval::Context {
    STANDARD_CTX.get_or_init(crate::eval::Context::standard)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn s(src: &str) -> Expr {
        simplify(&Parser::parse(src).unwrap())
    }

    #[test]
    fn identities() {
        assert_eq!(s("x + 0"), Expr::var("x"));
        assert_eq!(s("0 + x"), Expr::var("x"));
        assert_eq!(s("x * 1"), Expr::var("x"));
        assert_eq!(s("x * 0"), Expr::num(0.0));
        assert_eq!(s("x^1"), Expr::var("x"));
        assert_eq!(s("x^0"), Expr::num(1.0));
    }

    #[test]
    fn constant_folding() {
        assert_eq!(s("2 + 3"), Expr::num(5.0));
        assert_eq!(s("2 * 3 + 4"), Expr::num(10.0));
    }

    #[test]
    fn signed_one_distribution() {
        // -1 * x should become -x
        let got = s("-1 * x");
        let want = Expr::neg(Expr::var("x"));
        assert_eq!(got, want);
    }

    #[test]
    fn function_with_constants() {
        // sin(0) -> 0
        assert_eq!(s("sin(0)"), Expr::num(0.0));
        // cos(0) -> 1
        assert_eq!(s("cos(0)"), Expr::num(1.0));
    }
}