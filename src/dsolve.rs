//! Symbolic solution of first-order ordinary differential equations.
//!
//! Solves the standard first-order patterns **symbolically**, producing an
//! explicit expression for `y` in terms of `x` and the constant of
//! integration `C`:
//!
//! | Pattern | Example | Solution |
//! |---|---|---|
//! | Separable | `y' = x·y` | `y = C·e^(x²/2)` |
//! | Linear | `y' + 2x·y = x` | `y = ½ + C·e^(−x²)` |
//! | Exact derivative | `y' = 2x` | `y = x² + C` |
//!
//! The equation is given as `y' = <rhs>` (or `dy/dx = <rhs>`), with the
//! dependent variable named `y` and the independent one `x` by default.
//!
//! Everything is built from the existing symbolic core: `differentiate`,
//! `integrate` and `simplify`. The solution is verified by differentiating it
//! and comparing with the original right-hand side.
//!
//! Cross-port note: this is the "dsolve" of Qalculate. The `unsolve`
//! TypeScript port does not have it.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::simplify::simplify;
use crate::symbolic::{differentiate, integrate};

/// The result of solving a differential equation.
#[derive(Debug, Clone)]
pub struct Solution {
    /// The general solution `y(x)`, with the constant `C` as a free symbol.
    pub y: Expr,
    /// The pattern that was recognised.
    pub kind: SolutionKind,
    /// The intermediate steps, for the step-by-step view.
    pub steps: Vec<String>,
}

/// Which first-order pattern was recognised.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionKind {
    /// `y' = f'(x)` — the right-hand side is the derivative of something.
    ExactDerivative,
    /// `y' = f(x)·g(y)` — separable.
    Separable,
    /// `y' + P(x)·y = Q(x)` — linear.
    Linear,
}

/// Solve `y' = rhs` for `y`, with `x` the independent variable.
///
/// The right-hand side may mention `y` (then a separable or linear pattern
/// is required) or not (then it is integrated directly).
pub fn solve_ode(rhs: &Expr, x: &str, y: &str) -> Result<Solution> {
    let mut steps: Vec<String> = Vec::new();
    let rhs = simplify(rhs);

    // The right-hand side may mention y.
    if !rhs.variables().iter().any(|v| v == y) {
        // y' = f(x)  →  y = ∫f dx + C
        let antideriv = integrate(&rhs, x)?;
        let y_expr = simplify(&Expr::add(antideriv, Expr::var("C")));
        steps.push(format!("y' = {} (no dependence on {})", rhs, y));
        steps.push(format!("y = ∫{} d{} + C", rhs, x));
        steps.push(format!("y = {}", y_expr));
        return Ok(Solution {
            y: y_expr,
            kind: SolutionKind::ExactDerivative,
            steps,
        });
    }

    // Try the separable pattern first: y' = f(x)·g(y).
    // This must come before the linear test, because `x·y` is *also*
    // affine in y (with coefficient x) — but it is genuinely separable, and
    // the separable solution is the natural one.
    if let Some(sol) = solve_separable(&rhs, x, y, &mut steps)? {
        return Ok(sol);
    }

    // Try the linear pattern: y' + P(x)·y = Q(x).
    // The equation is y' = rhs, so it is linear when rhs can be written as
    // Q(x) − P(x)·y.
    if let Some(sol) = solve_linear(&rhs, x, y, &mut steps)? {
        return Ok(sol);
    }

    Err(MathError::Eval(format!(
        "dsolve: cannot solve y' = {} — recognised the separable and linear first-order \
         patterns only",
        rhs
    )))
}

/// Solve the linear pattern `y' = Q(x) − P(x)·y`.
///
/// Returns `None` when the right-hand side is not of that shape.
fn solve_linear(
    rhs: &Expr,
    x: &str,
    y: &str,
    steps: &mut Vec<String>,
) -> Result<Option<Solution>> {
    // The right-hand side must be affine in y: Q(x) + c(y)·y where c(y)
    // does not mention y.
    let (q, p) = match split_affine_in_y(rhs, x, y) {
        Some(v) => v,
        None => return Ok(None),
    };
    // P(x) = −c(y); the integrating factor is e^(∫P dx).
    let int_p = match integrate(&p, x) {
        Ok(v) => v,
        // The integrating factor has no closed form — decline the pattern.
        Err(_) => return Ok(None),
    };
    let mu = simplify(&Expr::func("exp", vec![int_p.clone()]));
    // y = (1/mu)·(∫ mu·Q dx + C)
    let integrand = simplify(&Expr::mul(mu.clone(), q.clone()));
    let int_mu_q = match integrate(&integrand, x) {
        Ok(v) => v,
        // No closed form for the transformed integrand — decline.
        Err(_) => return Ok(None),
    };
    let y_expr = simplify(&Expr::div(
        Expr::add(int_mu_q, Expr::var("C")),
        mu.clone(),
    ));
    steps.push(format!("linear: y' = {} − ({})·{}", q, p, y));
    steps.push(format!("integrating factor: e^(∫{} d{}) = {}", p, x, mu));
    steps.push(format!("y = {}", y_expr));
    Ok(Some(Solution {
        y: y_expr,
        kind: SolutionKind::Linear,
        steps: steps.clone(),
    }))
}

/// Split `rhs` into `(Q(x), P(x))` when it is affine in `y`.
///
/// Returns `Some((q, p))` when `rhs = q(x) − p(x)·y` for expressions that do
/// not mention `y`.
fn split_affine_in_y(rhs: &Expr, x: &str, y: &str) -> Option<(Expr, Expr)> {
    // The right-hand side may be a sum: separate the part that mentions y
    // from the part that does not.
    let (y_part, q_part): (Expr, Expr) = match rhs {
        Expr::Add(a, b) => {
            if mentions(a, y) && !mentions(b, y) {
                ((**a).clone(), (**b).clone())
            } else if mentions(b, y) && !mentions(a, y) {
                ((**b).clone(), (**a).clone())
            } else {
                return None;
            }
        }
        other => {
            // A single term: it must be a multiple of y.
            (other.clone(), Expr::num(0.0))
        }
    };
    // The y-part must be exactly `y` times a coefficient free of y.
    // The coefficient may be a nested product: `−2·x·y` has the coefficient
    // `−2·x`. A negated y (`−y`) has the coefficient `−1`.
    let coef = match &y_part {
        Expr::Var(v) if v == y => Expr::num(1.0),
        Expr::Neg(inner) => {
            let c = match inner.as_ref() {
                Expr::Var(v) if v == y => Expr::num(1.0),
                _ => {
                    let (c, _) = split_y_factor(inner, y)?;
                    c
                }
            };
            Expr::neg(c)
        }
        Expr::Mul(_, _) => {
            // Pull the `y` factor out of the product.
            let (coef, rest) = split_y_factor(&y_part, y)?;
            if mentions(&rest, y) {
                return None;
            }
            coef
        }
        _ => return None,
    };
    if mentions(&coef, y) {
        return None;
    }
    let _ = x;
    // rhs = q + coef·y  →  Q = q, P = −coef
    Some((simplify(&q_part), simplify(&Expr::neg(coef))))
}

/// Split a product into `(coefficient, rest)` where `rest` is the factor
/// `y` and the coefficient is everything else.
///
/// `−2·x·y` → `(−2·x, y)`, `2·y` → `(2, y)`, `y` → `(1, y)`.
/// Returns `None` when the product does not contain `y` exactly once.
fn split_y_factor(e: &Expr, y: &str) -> Option<(Expr, Expr)> {
    // Collect the product's factors, then find the one that is `y`.
    let mut factors: Vec<&Expr> = Vec::new();
    collect_factors(e, &mut factors);
    let y_count = factors
        .iter()
        .filter(|f| matches!(f, Expr::Var(v) if v == y))
        .count();
    if y_count != 1 {
        return None;
    }
    let mut coef = Expr::num(1.0);
    for f in &factors {
        if matches!(f, Expr::Var(v) if v == y) {
            continue;
        }
        if mentions(f, y) {
            // A factor that contains y in a deeper position (e.g. y²) is not
            // linear.
            return None;
        }
        coef = Expr::mul(coef, (*f).clone());
    }
    Some((coef, Expr::var(y.to_string())))
}

/// Flatten a product into its factors.
fn collect_factors<'a>(e: &'a Expr, out: &mut Vec<&'a Expr>) {
    match e {
        Expr::Mul(a, b) => {
            collect_factors(a, out);
            collect_factors(b, out);
        }
        other => out.push(other),
    }
}

/// Solve the separable pattern `y' = f(x)·g(y)`.
///
/// Returns `None` when the right-hand side does not factor that way.
fn solve_separable(
    rhs: &Expr,
    x: &str,
    y: &str,
    steps: &mut Vec<String>,
) -> Result<Option<Solution>> {
    // The right-hand side must be a product of an x-only part and a y-only
    // part.
    let (f, g) = match split_separable(rhs, x, y) {
        Some(v) => v,
        None => return Ok(None),
    };
    // ∫ dy/g(y) = ∫f(x) dx + C  →  y = g⁻¹(∫f dx + C)
    let int_f = integrate(&f, x)?;
    // The antiderivative of 1/g(y) with respect to y.
    let g_inv = Expr::div(Expr::num(1.0), g.clone());
    let int_g = integrate(&g_inv, y)?;
    // Solve for y: the antiderivative is hard to invert in general, so the
    // answer is stated in the implicit form plus the explicit one when the
    // integrand is invertible by inspection.
    steps.push(format!("separable: y' = ({})·({})", f, g));
    steps.push(format!("∫dy/({}) = ∫{} d{} + C", g, f, x));
    steps.push(format!("∫dy/({}) = {} + C", g, int_f));
    let explicit = invert_antiderivative(&int_g, &int_f, y);
    match explicit {
        Some(y_expr) => {
            steps.push(format!("y = {}", y_expr));
            Ok(Some(Solution {
                y: y_expr,
                kind: SolutionKind::Separable,
                steps: steps.clone(),
            }))
        }
        None => {
            // Fall back to the implicit form.
            let implicit = simplify(&Expr::add(int_g, Expr::neg(int_f)));
            steps.push(format!("(implicit) {} = C", implicit));
            Ok(Some(Solution {
                y: implicit,
                kind: SolutionKind::Separable,
                steps: steps.clone(),
            }))
        }
    }
}

/// Split `rhs` into `(f(x), g(y))` when it factors that way.
fn split_separable(rhs: &Expr, x: &str, y: &str) -> Option<(Expr, Expr)> {
    match rhs {
        Expr::Mul(a, b) => {
            let (fx, gy) = if !mentions(a, y) && mentions(b, y) {
                ((**a).clone(), (**b).clone())
            } else if !mentions(b, y) && mentions(a, y) {
                ((**b).clone(), (**a).clone())
            } else {
                return None;
            };
            // The x-only part must not mention y, and the y-only part must
            // not mention x.
            if mentions(&fx, y) || mentions(&gy, x) {
                return None;
            }
            Some((fx, gy))
        }
        _ => None,
    }
}

/// Try to invert an antiderivative to get `y` explicitly.
///
/// Handles the common case where the antiderivative of `1/g(y)` is a simple
/// function that can be inverted: `ln(y)` → `exp(·)`, `y` → `·`, `1/y` → `1/·`.
fn invert_antiderivative(int_g: &Expr, int_f: &Expr, y: &str) -> Option<Expr> {
    // The antiderivative must be exactly a function of y (plus a constant).
    let inner = strip_constant(int_g, y)?;
    match &inner {
        // ∫dy/y = ln|y| → y = ±exp(·); the sign is absorbed into C.
        Expr::Func(name, args) if name == "ln" && args.len() == 1 => {
            let arg = &args[0];
            if !mentions(arg, y) {
                return None;
            }
            let e = Expr::func("exp", vec![int_f.clone()]);
            Some(simplify(&Expr::mul(e, Expr::var("C"))))
        }
        // ∫dy = y → y = ∫f dx + C
        Expr::Var(v) if v == y => Some(simplify(&Expr::add(int_f.clone(), Expr::var("C")))),
        _ => None,
    }
}

/// Remove an additive constant from an expression, leaving the part that
/// mentions `y`.
///
/// Returns `None` when the expression mixes `y` and a constant additively.
fn strip_constant(e: &Expr, y: &str) -> Option<Expr> {
    match e {
        Expr::Add(a, b) => {
            if !mentions(a, y) && mentions(b, y) {
                Some((**b).clone())
            } else if !mentions(b, y) && mentions(a, y) {
                Some((**a).clone())
            } else {
                None
            }
        }
        Expr::Neg(a) => {
            let inner = strip_constant(a, y)?;
            Some(Expr::neg(inner))
        }
        _ => Some(e.clone()),
    }
}

/// Does this expression mention the given variable?
fn mentions(e: &Expr, var: &str) -> bool {
    e.variables().iter().any(|v| v == var)
}

/// Verify a solution by differentiating it and comparing with the
/// right-hand side.
///
/// The comparison is numerical: the derivative of the solution and the
/// right-hand side are both evaluated at sample points of `x`, with the
/// dependent variable `y` evaluated from the solution itself (so the
/// identity is checked along the solution curve, not at an arbitrary `y`).
///
/// The free constant `C` is bound to a sample value.
pub fn verify(solution: &Expr, rhs: &Expr, x: &str, y: &str) -> bool {
    let Ok(deriv) = differentiate(solution, x) else {
        return false;
    };
    let mut ctx = crate::eval::Context::standard();
    ctx.set("C", 1.0);
    for &xv in &[0.3f64, 1.1, 2.7] {
        ctx.set(x, xv);
        // The right-hand side may mention y; evaluate it at the solution's
        // own value so the identity is checked along the curve.
        if let Ok(yv) = crate::eval::eval(solution, &ctx) {
            ctx.set(y, yv);
        }
        let (Ok(d), Ok(r)) = (
            crate::eval::eval(&deriv, &ctx),
            crate::eval::eval(rhs, &ctx),
        ) else {
            return false;
        };
        if !d.is_finite() || !r.is_finite() {
            continue;
        }
        if (d - r).abs() > 1e-6 * r.abs().max(1.0) {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn sol(rhs_src: &str) -> Solution {
        let rhs = Parser::parse(rhs_src).unwrap();
        solve_ode(&rhs, "x", "y").unwrap()
    }

    #[test]
    fn solves_without_y_dependence() {
        // y' = 2x → y = x² + C
        let s = sol("2*x");
        assert_eq!(s.kind, SolutionKind::ExactDerivative);
        assert!(s.y.to_string().contains("x^2"), "got {}", s.y);
    }

    #[test]
    fn solves_separable() {
        // y' = x·y → y = C·e^(x²/2)
        let s = sol("x*y");
        assert_eq!(s.kind, SolutionKind::Separable);
        let text = s.y.to_string();
        assert!(text.contains("exp"), "got {}", text);
    }

    #[test]
    fn solves_linear() {
        // y' + y = 0  →  y = C·e^(−x)
        let s = sol("-y");
        assert_eq!(s.kind, SolutionKind::Linear);
        let text = s.y.to_string();
        assert!(text.contains("exp"), "got {}", text);
    }

    #[test]
    fn separable_solution_verifies() {
        // Differentiating the solution must reproduce the right-hand side.
        let s = sol("x*y");
        assert!(verify(&s.y, &Parser::parse("x*y").unwrap(), "x", "y"), "got {}", s.y);
    }

    #[test]
    fn linear_solution_verifies() {
        let rhs = Parser::parse("-y").unwrap();
        let s = solve_ode(&rhs, "x", "y").unwrap();
        assert!(verify(&s.y, &rhs, "x", "y"), "got {}", s.y);
    }

    #[test]
    fn exact_derivative_solution_verifies() {
        let rhs = Parser::parse("2*x").unwrap();
        let s = solve_ode(&rhs, "x", "y").unwrap();
        assert!(verify(&s.y, &rhs, "x", "y"), "got {}", s.y);
    }

    #[test]
    fn rejects_an_unsupported_pattern() {
        // y' = y² is separable, but y' = y·x·(x+1) with a mixed
        // non-separable shape is not. Use a genuinely unsupported one.
        let rhs = Parser::parse("y*exp(y)").unwrap();
        let err = solve_ode(&rhs, "x", "y").unwrap_err();
        assert!(err.to_string().contains("dsolve"), "got {}", err);
    }

    #[test]
    fn steps_are_recorded() {
        let s = sol("x*y");
        assert!(s.steps.len() >= 3, "got {:?}", s.steps);
        assert!(s.steps[0].contains("separable"), "got {:?}", s.steps);
    }
}
