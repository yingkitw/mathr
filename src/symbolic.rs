//! Symbolic differentiation.
//!
//! Given an [`Expr`] and the name of the variable to differentiate with
//! respect to, returns a derivative expression. The result is run through
//! [`simplify`] so identities like `x*0 = 0`, `-1*x = -x`, and
//! `x*1 = x` are applied.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::simplify::simplify;

/// Differentiate `expr` with respect to `var`. Returns the simplified form
/// so the output is friendlier to read in a REPL.
pub fn differentiate(expr: &Expr, var: &str) -> Result<Expr> {
    Ok(simplify(&diff(expr, var)?))
}

fn diff(expr: &Expr, var: &str) -> Result<Expr> {
    match expr {
        Expr::Num(_) => Ok(Expr::num(0.0)),
        Expr::Var(v) => {
            if v == var {
                Ok(Expr::num(1.0))
            } else {
                Ok(Expr::num(0.0))
            }
        }
        Expr::Neg(e) => Ok(Expr::neg(diff(e, var)?)),
        Expr::Add(a, b) => Ok(Expr::add(diff(a, var)?, diff(b, var)?)),
        Expr::Sub(a, b) => Ok(Expr::sub(diff(a, var)?, diff(b, var)?)),
        Expr::Mul(a, b) => {
            // (fg)' = f'g + fg'
            let da = diff(a, var)?;
            let db = diff(b, var)?;
            Ok(Expr::add(
                Expr::mul(da, (**b).clone()),
                Expr::mul((**a).clone(), db),
            ))
        }
        Expr::Div(a, b) => {
            // (f/g)' = (f'g - fg') / g^2
            let da = diff(a, var)?;
            let db = diff(b, var)?;
            Ok(Expr::div(
                Expr::sub(
                    Expr::mul(da, (**b).clone()),
                    Expr::mul((**a).clone(), db),
                ),
                Expr::pow((**b).clone(), Expr::num(2.0)),
            ))
        }
        Expr::Pow(base, exp) => match (base.as_ref(), exp.as_ref()) {
            // f(x)^c
            (_, Expr::Num(c)) => {
                let c_val = *c;
                let inner = Expr::pow((**base).clone(), Expr::num(c_val - 1.0));
                Ok(Expr::mul(
                    Expr::mul(Expr::num(c_val), inner),
                    diff(base, var)?,
                ))
            }
            // c^g(x)
            (Expr::Num(_), _) => Ok(Expr::mul(
                expr.clone(),
                Expr::mul(
                    Expr::func("ln", vec![(**base).clone()]),
                    diff(exp, var)?,
                ),
            )),
            // f^g in general
            _ => Ok(Expr::mul(
                expr.clone(),
                Expr::add(
                    Expr::mul(diff(exp, var)?, Expr::func("ln", vec![(**base).clone()])),
                    Expr::mul(
                        Expr::div((**exp).clone(), (**base).clone()),
                        diff(base, var)?,
                    ),
                ),
            )),
        },
        Expr::Func(name, args) => {
            // Chain rule: d/dx f(g) = f'(g) * g'
            if args.len() != 1 {
                return Err(MathError::Eval(format!(
                    "cannot differentiate multi-arg function {}",
                    name
                )));
            }
            let arg = &args[0];
            let arg_diff = diff(arg, var)?;
            let deriv = derivative_of_builtin(name, arg)?;
            Ok(Expr::mul(deriv, arg_diff))
        }
    }
}

/// Return the derivative of the named elementary function applied to `arg`.
fn derivative_of_builtin(name: &str, arg: &Expr) -> Result<Expr> {
    let x = arg.clone();
    Ok(match name {
        "sin" => Expr::func("cos", vec![x]),
        "cos" => Expr::neg(Expr::func("sin", vec![x])),
        "tan" => Expr::div(
            Expr::num(1.0),
            Expr::pow(Expr::func("cos", vec![x]), Expr::num(2.0)),
        ),
        "asin" => Expr::div(
            Expr::num(1.0),
            Expr::func(
                "sqrt",
                vec![Expr::sub(Expr::num(1.0), Expr::pow(x, Expr::num(2.0)))],
            ),
        ),
        "acos" => Expr::neg(Expr::div(
            Expr::num(1.0),
            Expr::func(
                "sqrt",
                vec![Expr::sub(Expr::num(1.0), Expr::pow(x, Expr::num(2.0)))],
            ),
        )),
        "atan" => Expr::div(
            Expr::num(1.0),
            Expr::add(Expr::num(1.0), Expr::pow(x, Expr::num(2.0))),
        ),
        "sinh" => Expr::func("cosh", vec![x]),
        "cosh" => Expr::func("sinh", vec![x]),
        "tanh" => Expr::div(
            Expr::num(1.0),
            Expr::pow(Expr::func("cosh", vec![x]), Expr::num(2.0)),
        ),
        "exp" => Expr::func("exp", vec![x]),
        "ln" | "log" => Expr::div(Expr::num(1.0), x),
        "log2" => Expr::div(
            Expr::num(1.0),
            Expr::mul(x, Expr::func("ln", vec![Expr::num(2.0)])),
        ),
        "log10" => Expr::div(
            Expr::num(1.0),
            Expr::mul(x, Expr::func("ln", vec![Expr::num(10.0)])),
        ),
        "sqrt" => Expr::div(
            Expr::num(1.0),
            Expr::mul(Expr::num(2.0), Expr::func("sqrt", vec![x.clone()])),
        ),
        "abs" => Expr::div(x.clone(), Expr::func("abs", vec![x])),
        "floor" | "ceil" | "round" | "sign" | "fract" => Expr::num(0.0),
        _ => {
            return Err(MathError::Eval(format!(
                "no symbolic derivative for function '{}'",
                name
            )))
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::{eval, Context};
    use crate::parser::Parser;

    /// Compare two expressions by evaluating them at the supplied `x` values.
    /// Simpler and more robust than AST equality once you go past the simplest
    /// of cases — the symbolic output is rarely byte-identical to the "nice"
    /// textbook form.
    fn agrees(got: &Expr, want: &Expr, xs: &[f64]) {
        let mut ctx = Context::standard();
        for &x in xs {
            ctx.set("x", x);
            let g = eval(got, &ctx).unwrap();
            let w = eval(want, &ctx).unwrap();
            assert!(
                (g - w).abs() < 1e-9,
                "disagree at x={}: got {} want {}",
                x,
                g,
                w
            );
        }
    }

    fn d(s: &str) -> Expr {
        let e = Parser::parse(s).unwrap();
        differentiate(&e, "x").unwrap()
    }

    #[test]
    fn polynomial() {
        // d/dx (x^3 + 2x^2 + x + 5) = 3x^2 + 4x + 1
        let got = d("x^3 + 2*x^2 + x + 5");
        let want = Parser::parse("3*x^2 + 4*x + 1").unwrap();
        agrees(&got, &want, &[0.0, 1.0, -2.0, 3.5, 0.7]);
    }

    #[test]
    fn product_rule() {
        // d/dx (x * sin(x)) = sin(x) + x*cos(x)
        let got = d("x*sin(x)");
        let want = Parser::parse("sin(x) + x*cos(x)").unwrap();
        agrees(&got, &want, &[0.1, 0.5, 1.0, 2.0, -0.7]);
    }

    #[test]
    fn quotient_rule() {
        // d/dx (x / (x+1)) = 1/(x+1)^2
        let got = d("x/(x+1)");
        let want = Parser::parse("1/(x+1)^2").unwrap();
        agrees(&got, &want, &[0.5, 1.5, 2.0, -0.5]);
    }

    #[test]
    fn chain_rule() {
        // d/dx sin(x^2) = 2*x*cos(x^2)
        let got = d("sin(x^2)");
        let want = Parser::parse("2*x*cos(x^2)").unwrap();
        agrees(&got, &want, &[0.1, 0.5, 1.0, -0.7]);
    }

    #[test]
    fn exp_ln() {
        // d/dx exp(x) = exp(x)
        let got = d("exp(x)");
        let want = Parser::parse("exp(x)").unwrap();
        agrees(&got, &want, &[0.0, 1.0, -2.0]);
    }
}