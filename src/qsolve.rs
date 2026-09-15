//! Symbolic solving of linear and quadratic equations.
//!
//! Pipeline: `poly::to_poly` extracts the coefficients (so `x^2 - 4`,
//! `(x - 1)(x + 2)`, and `x^2 = 2x + 3` all work — the RHS is moved over
//! by the caller), then the quadratic formula classifies by discriminant:
//! perfect-square discriminants give exact integer roots, negative ones
//! give complex-conjugate roots built as `Expr`s with the `i` literal
//! (evaluable by `ceval`), anything else snaps to 12 significant digits.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::poly::{to_poly, Poly, Term};

/// Snap values indistinguishable from integers, then round to 12 significant
/// digits (same convention as `limit.rs::snap_num`).
fn snap(v: f64) -> f64 {
    if !v.is_finite() {
        return v;
    }
    if (v - v.round()).abs() < 1e-9 * v.abs().max(1.0) {
        return v.round();
    }
    let mag = v.abs();
    let digits = 11 - mag.log10().floor() as i32;
    if !(-20..=20).contains(&digits) {
        return v;
    }
    let factor = 10f64.powi(digits);
    (v * factor).round() / factor
}

fn num(v: f64) -> Expr {
    Expr::num(snap(v))
}

/// Extract ascending coefficients `[c, b, a]` of a univariate polynomial in
/// `var` (degree <= 2). Returns `None` if the expression is not polynomial
/// in `var` alone, or has degree > 2.
fn coeffs(poly: &Poly, var: &str) -> Option<[f64; 3]> {
    let mut out = [0.0_f64; 3];
    for Term { coeff, monomial } in poly {
        match monomial.get(var).copied().unwrap_or(0) {
            0 if monomial.is_empty() => out[0] += coeff,
            1 if monomial.len() == 1 => out[1] += coeff,
            2 if monomial.len() == 1 => out[2] += coeff,
            _ => return None,
        }
    }
    Some(out)
}

/// Build `re + im*i` with sign-aware joins (never renders `+ -`).
fn complex_root_expr(re: f64, im: f64) -> Expr {
    let re = snap(re);
    let im = snap(im);
    let imag = |mag: f64| -> Expr {
        if mag == 1.0 {
            Expr::var("i")
        } else {
            Expr::mul(Expr::num(mag), Expr::var("i"))
        }
    };
    if im == 0.0 {
        return Expr::num(re);
    }
    if re == 0.0 {
        return if im < 0.0 { Expr::neg(imag(-im)) } else { imag(im) };
    }
    if im < 0.0 {
        Expr::sub(Expr::num(re), imag(-im))
    } else {
        Expr::add(Expr::num(re), imag(im))
    }
}

/// Solve `expr = 0` symbolically for `var`. Handles linear and quadratic
/// polynomials; returns the root expressions (1 or 2 of them).
pub fn solve_symbolic(expr: &Expr, var: &str) -> Result<Vec<Expr>> {
    let poly = to_poly(&crate::simplify::simplify(expr)).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "'{expr}' is not a polynomial in {var}; use the numeric solve instead"
        ))
    })?;
    let [c, b, a] = coeffs(&poly, var).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "degree > 2 in {var} is not supported symbolically; use poly-roots for numeric roots"
        ))
    })?;

    if a != 0.0 {
        let d = b * b - 4.0 * a * c;
        let d_scale = (b * b).abs().max((4.0 * a * c).abs()).max(1.0);
        let (r1, r2) = if d.abs() < 1e-12 * d_scale {
            // double root -b/(2a)
            let r = -b / (2.0 * a);
            (r, r)
        } else if d > 0.0 {
            let s = d.sqrt();
            ((-b + s) / (2.0 * a), (-b - s) / (2.0 * a))
        } else {
            let im = (-d).sqrt() / (2.0 * a);
            let re = -b / (2.0 * a);
            return Ok(vec![complex_root_expr(re, im), complex_root_expr(re, -im)]);
        };
        Ok(vec![Expr::num(snap(r1)), Expr::num(snap(r2))])
    } else if b != 0.0 {
        Ok(vec![Expr::num(snap(-c / b))])
    } else {
        Err(MathError::InvalidArgument(
            "equation has no variable dependence (constant = 0)".into(),
        ))
    }
}

/// REPL-facing steps: `qsolve <expr> [= rhs]`. The variable is inferred;
/// exactly one free variable is required (avoiding the trailing-token
/// ambiguity class of `fit`/`romberg`).
pub fn qsolve_steps_str(input: &str) -> Result<Vec<String>> {
    let (lhs, rhs) = match input.split_once('=') {
        Some((l, r)) => (l.trim(), Some(r.trim())),
        None => (input.trim(), None),
    };
    let e = crate::parser::Parser::parse(lhs)?;
    let e = match rhs {
        Some(r) => {
            let r = crate::parser::Parser::parse(r)?;
            Expr::sub(e, r)
        }
        None => e,
    };
    let mut vars = e.variables();
    vars.sort();
    vars.dedup();
    if vars.len() != 1 {
        return Err(MathError::InvalidArgument(format!(
            "qsolve needs exactly one variable (found: {}); numeric solve handles systems",
            if vars.is_empty() { "none".into() } else { vars.join(", ") }
        )));
    }
    let var = vars.remove(0);

    let poly = to_poly(&crate::simplify::simplify(&e)).ok_or_else(|| {
        MathError::InvalidArgument(format!("'{lhs}' is not a polynomial in {var}"))
    })?;
    let [c, b, a] = coeffs(&poly, &var).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "degree > 2 in {var} is not supported symbolically; use poly-roots for numeric roots"
        ))
    })?;

    let mut steps = Vec::new();
    let eq_display = match rhs {
        Some(r) => format!("{lhs} = {r}"),
        None => format!("{lhs} = 0"),
    };
    steps.push(format!("qsolve: {eq_display}  (in {var})"));

    let deg = if a != 0.0 {
        2
    } else if b != 0.0 {
        1
    } else {
        return Err(MathError::InvalidArgument(
            "equation has no variable dependence (constant = 0)".into(),
        ));
    };
    let fmt_coeff = |v: f64| -> String { format!("{}", num(v)) };
    if deg == 1 {
        steps.push(format!(
            "linear: {}*{} + {} = 0",
            fmt_coeff(b),
            var,
            fmt_coeff(c)
        ));
        let roots = solve_symbolic(&e, &var)?;
        steps.push(format!("{var} = {}", roots[0]));
        return Ok(steps);
    }
    let d = b * b - 4.0 * a * c;
    let d_scale = (b * b).abs().max((4.0 * a * c).abs()).max(1.0);
    steps.push(format!(
        "coefficients: a = {}, b = {}, c = {}",
        fmt_coeff(a),
        fmt_coeff(b),
        fmt_coeff(c)
    ));
    steps.push(format!(
        "discriminant: d = b^2 - 4ac = {} ({})",
        fmt_coeff(d),
        if d.abs() < 1e-12 * d_scale {
            "= 0 (double root)".to_string()
        } else if d > 0.0 {
            "> 0 (two real roots)".to_string()
        } else {
            "< 0 (complex-conjugate roots)".to_string()
        }
    ));
    let roots = solve_symbolic(&e, &var)?;
    let joined: Vec<String> = roots.iter().map(|r| format!("{var} = {r}")).collect();
    steps.push(joined.join(", "));
    Ok(steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn roots(input: &str) -> Vec<String> {
        qsolve_steps_str(input)
            .unwrap()
            .pop()
            .unwrap()
            .split(", ")
            .map(|s| s.to_string())
            .collect()
    }

    #[test]
    fn two_integer_roots() {
        assert_eq!(roots("x^2 - 4"), vec!["x = 2", "x = -2"]);
        assert_eq!(roots("x^2 = 2x + 3"), vec!["x = 3", "x = -1"]);
        assert_eq!(roots("x^2 - x - 6"), vec!["x = 3", "x = -2"]);
    }

    #[test]
    fn double_root() {
        assert_eq!(roots("(x - 3)^2"), vec!["x = 3", "x = 3"]);
        assert_eq!(roots("x^2"), vec!["x = 0", "x = 0"]);
    }

    #[test]
    fn complex_conjugate_roots() {
        assert_eq!(roots("x^2 + 1"), vec!["x = i", "x = -i"]);
        assert_eq!(roots("x^2 - 2x + 5"), vec!["x = 1 + 2*i", "x = 1 - 2*i"]);
        assert_eq!(roots("x^2 + 4x + 5"), vec!["x = -2 + i", "x = -2 - i"]);
    }

    #[test]
    fn linear() {
        assert_eq!(roots("2x - 6"), vec!["x = 3"]);
        assert_eq!(roots("5 - x"), vec!["x = 5"]);
    }

    #[test]
    fn irrational_discriminant_gives_decimals() {
        let r = roots("x^2 - 2");
        let x1: f64 = r[0].strip_prefix("x = ").unwrap().parse().unwrap();
        assert!((x1 - std::f64::consts::SQRT_2).abs() < 1e-10, "got {r:?}");
    }

    #[test]
    fn complex_roots_satisfy_equation_via_ceval() {
        // cross-module check: substitute 1 + 2i into x^2 - 2x + 5 -> 0
        let z = crate::ceval::eval_complex_str("(1 + 2*i)^2 - 2*(1 + 2*i) + 5").unwrap();
        assert!(z.re.abs() < 1e-9 && z.im.abs() < 1e-9, "got {z}");
    }

    #[test]
    fn rejects_non_polynomial_and_high_degree() {
        assert!(qsolve_steps_str("sin(x) = 0").is_err());
        assert!(qsolve_steps_str("x^3 - 1").is_err());
        assert!(qsolve_steps_str("x*y = 1").is_err()); // two variables
        assert!(qsolve_steps_str("1/x = 2").is_err()); // variable denominator
    }

    #[test]
    fn api_level_solve_symbolic() {
        let e = Parser::parse("x^2 - 2x + 5").unwrap();
        let rs = solve_symbolic(&e, "x").unwrap();
        assert_eq!(rs.len(), 2);
        assert_eq!(rs[0].to_string(), "1 + 2*i");
    }
}
