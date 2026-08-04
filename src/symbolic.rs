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

/// Compute the gradient of a multivariate expression.
///
/// Returns a vector of `(variable_name, partial_derivative)` pairs, one for
/// each free variable in `expr` (sorted alphabetically). Each partial
/// derivative is simplified.
pub fn gradient(expr: &Expr) -> Result<Vec<(String, Expr)>> {
    let vars = expr.variables();
    let mut result = Vec::with_capacity(vars.len());
    for v in &vars {
        result.push((v.clone(), differentiate(expr, v)?));
    }
    Ok(result)
}

/// Symbolic indefinite integration for the common elementary rules.
///
/// Handles:
/// - constants: ∫ c dx = c·x
/// - variable: ∫ x dx = x²/2
/// - powers: ∫ x^n dx = x^(n+1)/(n+1) for `n ≠ -1`, otherwise `ln(x)`
/// - sums / differences: linearity
/// - constant multiples: ∫ c·f dx = c·∫f dx
/// - exponentials: ∫ e^x dx = e^x,  ∫ a^x dx = a^x / ln(a)
/// - trig: ∫ sin(x) dx = −cos(x),  ∫ cos(x) dx = sin(x),  ∫ sec²(x) dx = tan(x)
/// - inverses: ∫ 1/x dx = ln(x),  ∫ 1/(1+x²) dx = atan(x),  ∫ 1/√(1−x²) dx = asin(x)
///
/// Returns `Err` for unsupported integrands (e.g., products of non-constants).
pub fn integrate(expr: &Expr, var: &str) -> Result<Expr> {
    Ok(simplify(&int_step(expr, var)?))
}

fn int_step(expr: &Expr, var: &str) -> Result<Expr> {
    match expr {
        Expr::Num(n) => Ok(Expr::mul(Expr::num(*n), Expr::var(var))),
        Expr::Var(v) if v == var => Ok(Expr::div(Expr::pow(Expr::var(var), Expr::num(2.0)), Expr::num(2.0))),
        Expr::Var(_) => Ok(Expr::mul(expr.clone(), Expr::var(var))),
        Expr::Neg(e) => Ok(Expr::neg(int_step(e, var)?)),
        Expr::Add(a, b) => Ok(Expr::add(int_step(a, var)?, int_step(b, var)?)),
        Expr::Sub(a, b) => Ok(Expr::sub(int_step(a, var)?, int_step(b, var)?)),
        Expr::Mul(a, b) => {
            if a.is_constant() {
                Ok(Expr::mul((**a).clone(), int_step(b, var)?))
            } else if b.is_constant() {
                Ok(Expr::mul(int_step(a, var)?, (**b).clone()))
            } else {
                Err(MathError::Eval(format!(
                    "integrate: cannot integrate non-linear product: {}",
                    expr
                )))
            }
        }
        Expr::Div(a, b) => {
            if b.is_constant() {
                Ok(Expr::mul(Expr::div(Expr::num(1.0), (**b).clone()), int_step(a, var)?))
            } else if a.is_constant() {
                integrate_constant_over(a, b, var)
            } else {
                Err(MathError::Eval(format!(
                    "integrate: cannot integrate non-constant numerator over non-constant denominator: {}",
                    expr
                )))
            }
        }
        Expr::Pow(base, exp) => match (base.as_ref(), exp.as_ref()) {
            (Expr::Var(v), Expr::Num(n)) if v == var && *n == -1.0 => {
                Ok(Expr::func("ln", vec![Expr::var(var)]))
            }
            (Expr::Var(v), Expr::Num(n)) if v == var => {
                let np1 = *n + 1.0;
                Ok(Expr::div(
                    Expr::pow(Expr::var(var), Expr::num(np1)),
                    Expr::num(np1),
                ))
            }
            (Expr::Num(_), _) | (_, _) if (**base).is_constant() => {
                Ok(Expr::div(
                    expr.clone(),
                    Expr::func("ln", vec![(**base).clone()]),
                ))
            }
            _ => Err(MathError::Eval(format!(
                "integrate: cannot integrate power: {}",
                expr
            ))),
        },
        Expr::Func(name, args) if args.len() == 1 => {
            let inner = &args[0];
            let inner_is_var = matches!(inner, Expr::Var(v) if v == var);
            match (name.as_str(), inner_is_var) {
                ("exp", true) => Ok(Expr::func("exp", vec![Expr::var(var)])),
                ("sin", true) => Ok(Expr::neg(Expr::func("cos", vec![Expr::var(var)]))),
                ("cos", true) => Ok(Expr::func("sin", vec![Expr::var(var)])),
                ("sec", true) => Ok(Expr::func("ln", vec![Expr::add(
                    Expr::func("sec", vec![Expr::var(var)]),
                    Expr::func("tan", vec![Expr::var(var)]),
                )])),
                ("tan", true) => Ok(Expr::neg(Expr::func("ln", vec![Expr::func("cos", vec![Expr::var(var)])]))),
                _ => Err(MathError::Eval(format!(
                    "integrate: unsupported integrand `{}{}`", name, inner
                ))),
            }
        }
        Expr::Func(name, args) => Err(MathError::Eval(format!(
            "integrate: cannot integrate multi-argument function {} with {} args",
            name,
            args.len()
        ))),
    }
}

/// Handle integrands of the form `c / g(x)` where `c` is constant.
/// Recognises `1/x` and `1/(1+x²)` and `1/√(1−x²)`.
fn integrate_constant_over(num: &Expr, den: &Expr, var: &str) -> Result<Expr> {
    let num_val = if let Expr::Num(n) = num { *n } else { 1.0 };
    let _ = num_val;
    match den {
        Expr::Var(v) if v == var => Ok(Expr::mul(Expr::num(num_val), Expr::func("ln", vec![Expr::var(var)]))),
        Expr::Add(a, b) | Expr::Sub(a, b) => {
            // 1 / (1 ± x²) → atan or -atan
            let is_one = matches!(a.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
            let is_xsq = match b.as_ref() {
                Expr::Pow(p, e) => matches!(p.as_ref(), Expr::Var(v) if v == var)
                    && matches!(e.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12),
                _ => false,
            };
            if is_one && is_xsq {
                let sign = if matches!(den, Expr::Sub(..)) { -1.0 } else { 1.0 };
                Ok(Expr::mul(
                    Expr::num(sign * num_val),
                    Expr::func("atan", vec![Expr::var(var)]),
                ))
            } else {
                Err(MathError::Eval(format!(
                    "integrate: unsupported integrand 1/{}", den
                )))
            }
        }
        Expr::Func(name, args) if name == "sqrt" && args.len() == 1 => {
            // 1 / sqrt(1 - x²) → asin
            if let Expr::Sub(a_inner, b_inner) = &args[0] {
                let is_one = matches!(a_inner.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
                let is_xsq = matches!(
                    b_inner.as_ref(),
                    Expr::Pow(pp, ee)
                        if matches!(pp.as_ref(), Expr::Var(v) if v == var)
                            && matches!(ee.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12)
                );
                if is_one && is_xsq {
                    return Ok(Expr::mul(
                        Expr::num(num_val),
                        Expr::func("asin", vec![Expr::var(var)]),
                    ));
                }
            }
            Err(MathError::Eval(format!(
                "integrate: unsupported integrand 1/{}",
                den
            )))
        }
        Expr::Pow(p, e) => {
            // 1 / (1 - x²)^{1/2} (as Pow rather than sqrt) → asin
            if let Expr::Sub(a_inner, b_inner) = p.as_ref() {
                let is_one = matches!(a_inner.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
                let is_xsq = matches!(
                    b_inner.as_ref(),
                    Expr::Pow(pp, ee)
                        if matches!(pp.as_ref(), Expr::Var(v) if v == var)
                            && matches!(ee.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12)
                );
                let is_half = matches!(e.as_ref(), Expr::Num(n) if (*n - 0.5).abs() < 1e-12);
                if is_one && is_xsq && is_half {
                    return Ok(Expr::mul(
                        Expr::num(num_val),
                        Expr::func("asin", vec![Expr::var(var)]),
                    ));
                }
            }
            Err(MathError::Eval(format!(
                "integrate: unsupported integrand 1/{}",
                den
            )))
        }
        _ => Err(MathError::Eval(format!(
            "integrate: unsupported integrand {}/{}", num, den
        ))),
    }
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

    fn integrate_agrees(src: &str, want: &str, xs: &[f64]) {
        // Differentiate the symbolic integral and verify it matches the
        // original integrand at sample points.
        let e = Parser::parse(src).unwrap();
        let antideriv = integrate(&e, "x").unwrap();
        let derived = differentiate(&antideriv, "x").unwrap();
        let want_e = Parser::parse(want).unwrap();
        agrees(&derived, &want_e, xs);
    }

    #[test]
    fn integrate_constant() {
        // ∫ 3 dx = 3x
        let e = Parser::parse("3").unwrap();
        let result = integrate(&e, "x").unwrap();
        let want = Parser::parse("3*x").unwrap();
        agrees(&result, &want, &[1.0, 2.0, -5.0]);
    }

    #[test]
    fn integrate_polynomial() {
        integrate_agrees("x", "x", &[0.5, 1.0, 2.0]);
        integrate_agrees("x^2", "x^2", &[0.5, 1.0, 2.0]);
        integrate_agrees("x^3 - 2*x + 1", "x^3 - 2*x + 1", &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn integrate_reciprocal() {
        // ∫ 1/x dx = ln(x)
        integrate_agrees("1/x", "1/x", &[0.5, 1.5, 3.0]);
    }

    #[test]
    fn integrate_exp() {
        integrate_agrees("exp(x)", "exp(x)", &[0.5, 1.0, 2.0]);
        integrate_agrees("2^x", "2^x", &[0.0, 1.0, 2.0]);
    }

    #[test]
    fn integrate_trig() {
        integrate_agrees("sin(x)", "sin(x)", &[0.5, 1.0, 2.0]);
        integrate_agrees("cos(x)", "cos(x)", &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn integrate_atan_arcsin() {
        integrate_agrees("1/(1+x^2)", "1/(1+x^2)", &[0.5, 1.0, 2.0]);
        integrate_agrees("1/sqrt(1-x^2)", "1/sqrt(1-x^2)", &[0.0, 0.3, 0.5]);
    }

    #[test]
    fn partial_derivative() {
        // ∂/∂x (x^2 * y + y^3) = 2*x*y
        let e = Parser::parse("x^2 * y + y^3").unwrap();
        let got = differentiate(&e, "x").unwrap();
        let want = Parser::parse("2*x*y").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at x={},y={}: got {} want {}", x, y, g, w);
        }
    }

    #[test]
    fn partial_derivative_other_var() {
        // ∂/∂y (x^2 * y + y^3) = x^2 + 3*y^2
        let e = Parser::parse("x^2 * y + y^3").unwrap();
        let got = differentiate(&e, "y").unwrap();
        let want = Parser::parse("x^2 + 3*y^2").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at x={},y={}: got {} want {}", x, y, g, w);
        }
    }

    #[test]
    fn gradient_multivariate() {
        // ∇(x^2 + x*y + y^2) = [∂/∂x = 2*x + y, ∂/∂y = x + 2*y]
        let e = Parser::parse("x^2 + x*y + y^2").unwrap();
        let grad = gradient(&e).unwrap();
        assert_eq!(grad.len(), 2);
        assert_eq!(grad[0].0, "x");
        assert_eq!(grad[1].0, "y");

        let want_dx = Parser::parse("2*x + y").unwrap();
        let want_dy = Parser::parse("x + 2*y").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let gdx = eval(&grad[0].1, &ctx).unwrap();
            let wdx = eval(&want_dx, &ctx).unwrap();
            assert!((gdx - wdx).abs() < 1e-9, "dx at x={},y={}: got {} want {}", x, y, gdx, wdx);
            let gdy = eval(&grad[1].1, &ctx).unwrap();
            let wdy = eval(&want_dy, &ctx).unwrap();
            assert!((gdy - wdy).abs() < 1e-9, "dy at x={},y={}: got {} want {}", x, y, gdy, wdy);
        }
    }

    #[test]
    fn gradient_three_vars() {
        // ∇(x*y*z) = [∂/∂x = y*z, ∂/∂y = x*z, ∂/∂z = x*y]
        let e = Parser::parse("x*y*z").unwrap();
        let grad = gradient(&e).unwrap();
        assert_eq!(grad.len(), 3);

        let wants = [
            Parser::parse("y*z").unwrap(),
            Parser::parse("x*z").unwrap(),
            Parser::parse("x*y").unwrap(),
        ];
        let mut ctx = Context::standard();
        for &(x, y, z) in &[(1.0, 2.0, 3.0), (0.5, -1.0, 2.0)] {
            ctx.set("x", x);
            ctx.set("y", y);
            ctx.set("z", z);
            for (i, want) in wants.iter().enumerate() {
                let g = eval(&grad[i].1, &ctx).unwrap();
                let w = eval(want, &ctx).unwrap();
                assert!((g - w).abs() < 1e-9, "var {} at ({},{},{}): got {} want {}", grad[i].0, x, y, z, g, w);
            }
        }
    }

    #[test]
    fn gradient_constant() {
        // ∇(42) = [] (no variables)
        let e = Parser::parse("42").unwrap();
        let grad = gradient(&e).unwrap();
        assert!(grad.is_empty());
    }
}