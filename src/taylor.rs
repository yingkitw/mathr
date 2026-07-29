//! Symbolic Taylor series expansion from scratch.
//!
//! Given an expression `f(x)`, computes the Taylor series around a point `a`
//! up to `n` terms: `f(a) + f'(a)(x-a) + f''(a)(x-a)^2/2! + ...`

use crate::error::{MathError, Result};
use crate::eval::{eval, Context};
use crate::expr::Expr;
use crate::parser::Parser;
use crate::simplify::simplify;
use crate::symbolic::differentiate;

/// Compute the Taylor series of `f` around `a` up to `order` terms.
///
/// Returns the symbolic expression for the series. Each term is
/// `f^(n)(a) / n! * (x - a)^n`.
pub fn taylor_series(f: &Expr, var: &str, a: f64, order: usize) -> Result<Expr> {
    if order == 0 {
        return Err(MathError::InvalidArgument("taylor_series: order must be > 0".into()));
    }

    let ctx = Context::standard();
    let mut terms: Vec<Expr> = Vec::with_capacity(order);
    let mut current = f.clone();

    for n in 0..order {
        // Evaluate the nth derivative at x = a
        let mut eval_ctx = ctx.clone();
        eval_ctx.set(var, a);
        let coeff = eval(&current, &eval_ctx).unwrap_or(0.0);

        if coeff.abs() < 1e-15 {
            // Skip zero terms but still need to differentiate for next iteration
            if n + 1 < order {
                current = differentiate(&current, var)?;
            }
            continue;
        }

        let factorial_n = factorial(n);
        let normalized = coeff / factorial_n as f64;

        // (x - a)^n
        let x_minus_a = if a == 0.0 {
            Expr::var(var)
        } else {
            Expr::sub(Expr::var(var), Expr::num(a))
        };
        let power = if n == 0 {
            Expr::num(1.0)
        } else if n == 1 {
            x_minus_a
        } else {
            Expr::pow(x_minus_a, Expr::num(n as f64))
        };

        let term = if normalized == 1.0 {
            power
        } else if normalized == -1.0 {
            Expr::neg(power)
        } else {
            Expr::mul(Expr::num(normalized), power)
        };
        terms.push(term);

        // Differentiate for next iteration
        if n + 1 < order {
            current = differentiate(&current, var)?;
        }
    }

    if terms.is_empty() {
        return Ok(Expr::num(0.0));
    }

    // Sum all terms
    let mut result = terms[0].clone();
    for t in &terms[1..] {
        result = Expr::add(result, t.clone());
    }

    Ok(simplify(&result))
}

/// Compute the Taylor series from a string expression.
pub fn taylor_series_str(src: &str, var: &str, a: f64, order: usize) -> Result<Expr> {
    let f = Parser::parse(src)?;
    taylor_series(&f, var, a, order)
}

fn factorial(n: usize) -> u64 {
    let mut r: u64 = 1;
    for i in 2..=n {
        r *= i as u64;
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    fn eval_at(e: &Expr, x: f64) -> f64 {
        let mut ctx = Context::standard();
        ctx.set("x", x);
        eval(e, &ctx).unwrap_or(f64::NAN)
    }

    #[test]
    fn taylor_exp_around_zero() {
        // e^x = 1 + x + x^2/2 + x^3/6 + ...
        let series = taylor_series_str("exp(x)", "x", 0.0, 5).unwrap();
        // Check that it approximates e^x well near 0
        for &x in &[0.0, 0.1, 0.5] {
            let approx = eval_at(&series, x);
            let exact = x.exp();
            assert!(close(approx, exact, 1e-3), "at x={}: got {} want {}", x, approx, exact);
        }
    }

    #[test]
    fn taylor_sin_around_zero() {
        // sin(x) = x - x^3/6 + x^5/120 - ...
        let series = taylor_series_str("sin(x)", "x", 0.0, 7).unwrap();
        for &x in &[0.0, 0.1, 0.5] {
            let approx = eval_at(&series, x);
            let exact = x.sin();
            assert!(close(approx, exact, 1e-4), "at x={}: got {} want {}", x, approx, exact);
        }
    }

    #[test]
    fn taylor_cos_around_zero() {
        // cos(x) = 1 - x^2/2 + x^4/24 - ...
        let series = taylor_series_str("cos(x)", "x", 0.0, 6).unwrap();
        for &x in &[0.0, 0.1, 0.5] {
            let approx = eval_at(&series, x);
            let exact = x.cos();
            assert!(close(approx, exact, 1e-4), "at x={}: got {} want {}", x, approx, exact);
        }
    }

    #[test]
    fn taylor_polynomial_exact() {
        // Taylor of a polynomial should be exact
        // f(x) = x^2 + 3x + 2, around a=1
        let series = taylor_series_str("x^2 + 3*x + 2", "x", 1.0, 4).unwrap();
        for &x in &[-1.0, 0.0, 1.0, 2.0, 5.0] {
            let approx = eval_at(&series, x);
            let exact = x * x + 3.0 * x + 2.0;
            assert!(close(approx, exact, 1e-10), "at x={}: got {} want {}", x, approx, exact);
        }
    }

    #[test]
    fn taylor_around_nonzero() {
        // e^x around a=1: should approximate e^x well near x=1
        let series = taylor_series_str("exp(x)", "x", 1.0, 6).unwrap();
        for &x in &[0.8, 1.0, 1.2] {
            let approx = eval_at(&series, x);
            let exact = x.exp();
            assert!(close(approx, exact, 1e-3), "at x={}: got {} want {}", x, approx, exact);
        }
    }
}
