//! Summation (Σ) and product (∏) over integer bounds.
//!
//! Evaluation strategy, in order:
//! 1. **Symbolic closed forms** (exact, O(1), work for huge ranges):
//!    var-independent expressions (`expr·count`), `var` (arithmetic series),
//!    `var^2`/`var^3` (Faulhaber), `c^var` (geometric).
//! 2. **Exact rational loop**: term-by-term `eval_rational` accumulation when
//!    every term is rational (capped — rational accumulation can overflow).
//! 3. **f64 loop** fallback.
//!
//! Empty ranges follow the standard convention: Σ = 0, ∏ = 1.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::rational::{eval_rational, Rational};

const MAX_TERMS: i64 = 1_000_000;
const EXACT_LOOP_CAP: i64 = 50_000;

fn as_f64(e: &Expr) -> Result<f64> {
    crate::eval::eval(e, &crate::eval::Context::standard())
}

/// Integer result → `Num(n)`; fractional → exact `Div(Num(n), Num(d))`
/// (displayed as `n/d`, re-parsable, and rational-evaluable).
fn rational_expr(r: Rational) -> Expr {
    if r.den() == 1 {
        Expr::num(r.num() as f64)
    } else {
        Expr::div(Expr::num(r.num() as f64), Expr::num(r.den() as f64))
    }
}

fn i128_expr(v: i128) -> Option<Expr> {
    if v.abs() <= 9_007_199_254_740_992 {
        Some(Expr::num(v as f64))
    } else {
        None // beyond f64 integers — caller falls back to f64
    }
}

fn closed_form_sum(expr: &Expr, var: &str, a: i64, b: i64) -> Option<Expr> {
    let n = (b - a + 1) as i128;
    match expr {
        // constant: expr · n
        e if !e.variables().contains(&var.to_string()) => {
            let v = as_f64(e).ok()?;
            let s = v * n as f64;
            if v.fract() == 0.0 {
                i128_expr((v as i128) * n).or(Some(Expr::num(s)))
            } else {
                Some(Expr::num(s))
            }
        }
        // Σ var = (a + b)·n / 2
        Expr::Var(v) if v == var => {
            let s = (a as i128 + b as i128) * n / 2;
            i128_expr(s)
        }
        // Σ var² = b(b+1)(2b+1)/6 − (a−1)a(2a−1)/6
        Expr::Pow(base, exponent)
            if matches!(**exponent, Expr::Num(2.0)) && matches!(**base, Expr::Var(ref v) if v.as_str() == var) =>
        {
            let f = |k: i128| k * (k + 1) * (2 * k + 1) / 6;
            i128_expr(f(b as i128) - f(a as i128 - 1))
        }
        // Σ var³ = (b(b+1)/2)² − ((a−1)a/2)²
        Expr::Pow(base, exponent)
            if matches!(**exponent, Expr::Num(3.0)) && matches!(**base, Expr::Var(ref v) if v.as_str() == var) =>
        {
            let f = |k: i128| {
                let t = k * (k + 1) / 2;
                t * t
            };
            i128_expr(f(b as i128) - f(a as i128 - 1))
        }
        // Σ c^var = (c^a − c^(b+1)) / (1 − c), c ≠ 1
        Expr::Pow(base, exponent)
            if matches!(**exponent, Expr::Var(ref v) if v.as_str() == var) && matches!(**base, Expr::Num(_)) =>
        {
            let Expr::Num(c) = &**base else { return None };
            if *c == 1.0 {
                return i128_expr(n);
            }
            let s = c.powi(a as i32) - c.powi(b as i32 + 1);
            let denom = 1.0 - c;
            Some(Expr::num(s / denom))
        }
        _ => None,
    }
}

fn closed_form_product(expr: &Expr, var: &str, a: i64, b: i64) -> Option<Expr> {
    let n = (b - a + 1) as i32;
    match expr {
        // constant: expr^n
        e if !e.variables().contains(&var.to_string()) => {
            let v = as_f64(e).ok()?;
            Some(Expr::num(v.powi(n)))
        }
        // c^var over consecutive integers = c^(n·(a+b)/2)
        Expr::Pow(base, exponent)
            if matches!(**exponent, Expr::Var(ref v) if v.as_str() == var) && matches!(**base, Expr::Num(_)) =>
        {
            let Expr::Num(c) = &**base else { return None };
            let total = (a + b) as f64 * n as f64 / 2.0;
            Some(Expr::num(c.powf(total)))
        }
        _ => None,
    }
}

fn generic_loop(expr: &Expr, var: &str, a: i64, b: i64, is_sum: bool) -> Result<Expr> {
    // f64 reference pass (also the fallback result)
    let mut ctx = crate::eval::Context::standard();
    let mut acc = if is_sum { 0.0 } else { 1.0 };
    for i in a..=b {
        ctx.set(var, i as f64);
        let t = crate::eval::eval(expr, &ctx)?;
        acc = if is_sum { acc + t } else { acc * t };
    }
    // exact rational accumulation when every term is rational. i64-backed
    // Rationals can silently overflow on long/denominator-heavy sums, so the
    // result must agree with the f64 pass — otherwise trust f64.
    if b - a + 1 <= EXACT_LOOP_CAP {
        let mut exact = Rational::from_int(if is_sum { 0 } else { 1 });
        let mut clean = true;
        for i in a..=b {
            let term_e = expr.substitute(var, &Expr::num(i as f64));
            match eval_rational(&term_e) {
                Some(r) => exact = if is_sum { exact + r } else { exact * r },
                None => {
                    clean = false;
                    break;
                }
            }
        }
        if clean {
            let r = exact.to_f64();
            if (r - acc).abs() <= 1e-6 * acc.abs().max(1.0) {
                return Ok(rational_expr(exact));
            }
        }
    }
    Ok(Expr::num(acc))
}

/// Σ `expr` for `var` = a..=b (inclusive).
pub fn summation(expr: &Expr, var: &str, a: i64, b: i64) -> Result<Expr> {
    if b - a + 1 > MAX_TERMS {
        return Err(MathError::InvalidArgument(format!(
            "range too large: {} terms (max {MAX_TERMS})",
            b - a + 1
        )));
    }
    if a > b {
        return Ok(Expr::num(0.0));
    }
    if let Some(e) = closed_form_sum(expr, var, a, b) {
        return Ok(e);
    }
    generic_loop(expr, var, a, b, true)
}

/// ∏ `expr` for `var` = a..=b (inclusive).
pub fn product(expr: &Expr, var: &str, a: i64, b: i64) -> Result<Expr> {
    if b - a + 1 > MAX_TERMS {
        return Err(MathError::InvalidArgument(format!(
            "range too large: {} terms (max {MAX_TERMS})",
            b - a + 1
        )));
    }
    if a > b {
        return Ok(Expr::num(1.0));
    }
    if let Some(e) = closed_form_product(expr, var, a, b) {
        return Ok(e);
    }
    generic_loop(expr, var, a, b, false)
}

/// REPL-facing parse: `<expr> [<var>] <a> <b>` with bounds taken from the end
/// and an optional single-letter variable token (minimize precedent — avoids
/// the trailing-token ambiguity class by design).
pub fn parse_bounds(rest: &str) -> Result<(String, String, i64, i64)> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(MathError::Eval("needs: <expr> [<var>] <a> <b>".into()));
    }
    let b: i64 = tokens[tokens.len() - 1]
        .parse()
        .map_err(|_| MathError::Eval(format!("bad upper bound: {}", tokens[tokens.len() - 1])))?;
    let a: i64 = tokens[tokens.len() - 2]
        .parse()
        .map_err(|_| MathError::Eval(format!("bad lower bound: {}", tokens[tokens.len() - 2])))?;
    let mut var = String::new();
    let expr_end = tokens.len() - 2;
    let expr_end = if expr_end >= 2 && tokens[expr_end - 1].len() == 1 && tokens[expr_end - 1].chars().all(|c| c.is_alphabetic()) {
        var = tokens[expr_end - 1].to_string();
        expr_end - 1
    } else {
        expr_end
    };
    let expr_src = tokens[..expr_end].join(" ");
    Ok((expr_src, var, a, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn value(e: &Expr) -> f64 {
        crate::eval::eval(e, &crate::eval::Context::standard()).unwrap()
    }

    fn sum_of(input: &str) -> Expr {
        let body = input.strip_prefix("sum ").unwrap_or(input);
        let (src, var, a, b) = parse_bounds(body).unwrap();
        let e = Parser::parse(src.trim()).unwrap();
        let var = if var.is_empty() {
            let mut vs = e.variables();
            if vs.len() == 1 { vs.remove(0) } else { String::new() }
        } else {
            var
        };
        summation(&e, &var, a, b).unwrap()
    }

    fn prod_of(input: &str) -> Expr {
        let body = input.strip_prefix("prod ").unwrap_or(input);
        let (src, var, a, b) = parse_bounds(body).unwrap();
        let e = Parser::parse(src.trim()).unwrap();
        let var = if var.is_empty() {
            let mut vs = e.variables();
            if vs.len() == 1 { vs.remove(0) } else { String::new() }
        } else {
            var
        };
        product(&e, &var, a, b).unwrap()
    }

    #[test]
    fn arithmetic_series_closed_form() {
        assert_eq!(value(&sum_of("sum x 1 100")), 5050.0);
        // huge range — exact via formula, no loop
        assert_eq!(value(&sum_of("sum x 1 1000000")), 500000500000.0);
        assert_eq!(value(&sum_of("sum k 5 10")), 45.0);
    }

    #[test]
    fn faulhaber_closed_forms() {
        assert_eq!(value(&sum_of("sum x^2 1 10")), 385.0);
        assert_eq!(value(&sum_of("sum x^3 1 3")), 36.0);
        // 3.33e17 exceeds f64's exact-integer range — compare within float precision
        let got = value(&sum_of("sum x^2 1 1000000"));
        assert!((got - 333333833333500000.0).abs() <= 1e-9 * 333333833333500000.0, "got {got}");
    }

    #[test]
    fn constant_and_geometric() {
        assert_eq!(value(&sum_of("sum 1 1 50")), 50.0);
        assert_eq!(value(&sum_of("sum 2^x 0 10")), 2047.0);
        assert_eq!(value(&prod_of("prod 2 1 10")), 1024.0);
        assert_eq!(value(&prod_of("prod 3^x 0 4")), 3_f64.powi(10));
    }

    #[test]
    fn exact_rational_loop() {
        // Σ 1/k for k = 1..10 = 7381/2520 exactly
        let e = sum_of("sum 1/k 1 10");
        assert_eq!(e.to_string(), "7381/2520", "got {e}");
        assert!((value(&e) - 2.9289682539682538).abs() < 1e-12);
    }

    #[test]
    fn product_loop() {
        assert_eq!(value(&prod_of("prod x 1 5")), 120.0);
        assert_eq!(value(&prod_of("prod x 3 5")), 60.0);
        // empty product = 1
        assert_eq!(value(&prod_of("prod x 5 1")), 1.0);
    }

    #[test]
    fn empty_sum_is_zero() {
        assert_eq!(value(&sum_of("sum x 5 1")), 0.0);
    }

    #[test]
    fn float_fallback() {
        let e = sum_of("sum sin(x) 1 5");
        let expected: f64 = (1..=5).map(|i| (i as f64).sin()).sum();
        assert!((value(&e) - expected).abs() < 1e-10);
    }

    #[test]
    fn explicit_var_token() {
        // "sum k^2 k 1 10" — var k explicitly, bounds from the end
        let e = sum_of("sum k^2 k 1 10");
        assert_eq!(value(&e), 385.0);
        // explicit var that the expression doesn't contain, with unbound x:
        // honest error (NOT silently treated as a constant)
        let body = "x^2 k 1 10";
        let (src, var, a, b) = parse_bounds(body).unwrap();
        let e = Parser::parse(src.trim()).unwrap();
        assert!(summation(&e, &var, a, b).is_err());
    }

    #[test]
    fn range_cap_rejects_huge() {
        let e = Parser::parse("x").unwrap();
        assert!(summation(&e, "x", 1, 2_000_000).is_err());
    }
}
