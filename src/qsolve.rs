//! Symbolic solving of univariate polynomial equations (degree ≤ 4).
//!
//! Pipeline: `poly::to_poly` extracts coefficients (expanded products and
//! `lhs = rhs` forms work once the RHS is moved left). For degree ≤ 2 the
//! quadratic formula applies (exact rationals when the discriminant is a
//! perfect square). For degree 3–4 with integer coefficients, rational-root
//! candidates are tested and the polynomial is deflated in exact `Rational`
//! arithmetic before finishing at degree ≤ 2. Remaining cases use
//! Durand–Kerner numeric roots, snapped for display and checked against the
//! original polynomial.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::poly::{to_poly, Poly, Term};
use crate::rational::Rational;
use crate::solver::polynomial_roots_complex;

const MAX_QSOLVE_DEG: usize = 4;

fn num(v: f64) -> Expr {
    Rational::snap_coeff_expr(v, crate::rational::SNAP_MAX_DEN_COMPACT)
}

/// Coefficient is an exact integer (polynomial coeffs from `to_poly` are f64).
fn near_int(v: f64) -> Option<i64> {
    if !v.is_finite() {
        return None;
    }
    let r = v.round();
    if (v - r).abs() <= 1e-9 * v.abs().max(1.0) {
        Some(r as i64)
    } else {
        None
    }
}

fn perfect_square_i64(n: i64) -> Option<i64> {
    if n < 0 {
        return None;
    }
    let r = (n as f64).sqrt().round() as i64;
    if r * r == n {
        Some(r)
    } else {
        None
    }
}

fn expr_from_rational(r: Rational) -> Expr {
    r.to_expr()
}

/// When `a`, `b`, `c` are exact integers and `d = b² − 4ac` is a perfect square,
/// return roots as exact rationals (e.g. `6x² − 5x − 1` → `1`, `-1/6`).
/// Cross-port: unsolve `solver.ts` exact quadratic path.
fn try_exact_rational_quadratic_roots(a: f64, b: f64, c: f64) -> Option<Vec<Expr>> {
    let ai = near_int(a)?;
    let bi = near_int(b)?;
    let ci = near_int(c)?;
    if ai == 0 {
        return None;
    }
    let d128 = (bi as i128) * (bi as i128) - 4 * (ai as i128) * (ci as i128);
    if d128 < 0 {
        return None;
    }
    if d128 > i64::MAX as i128 {
        return None;
    }
    let d = d128 as i64;
    let sd = perfect_square_i64(d)?;
    let two_a = 2 * ai;
    if two_a == 0 {
        return None;
    }
    let r1 = Rational::new(-bi + sd, two_a).ok()?;
    let r2 = Rational::new(-bi - sd, two_a).ok()?;
    Some(vec![expr_from_rational(r1), expr_from_rational(r2)])
}

/// Ascending coefficients `[c0, c1, …, cn]` for `c0 + c1·var + … + cn·var^n`.
fn asc_coeffs(poly: &Poly, var: &str) -> Option<Vec<f64>> {
    let mut out = vec![0.0_f64; MAX_QSOLVE_DEG + 1];
    for Term { coeff, monomial } in poly {
        if monomial.is_empty() {
            out[0] += coeff;
            continue;
        }
        if monomial.len() != 1 {
            return None;
        }
        let (name, exp) = monomial.iter().next()?;
        if name != var {
            return None;
        }
        let exp = *exp as usize;
        if exp > MAX_QSOLVE_DEG {
            return None;
        }
        out[exp] += coeff;
    }
    while out.len() > 1 && out.last().copied().unwrap_or(0.0).abs() == 0.0 {
        out.pop();
    }
    if out.is_empty() {
        out.push(0.0);
    }
    Some(out)
}

fn poly_degree(asc: &[f64]) -> usize {
    asc.len().saturating_sub(1)
}

fn divisors_pos(n: i64) -> Vec<i64> {
    let n = n.abs().max(1);
    let mut out = Vec::new();
    let lim = (n as f64).sqrt().floor() as i64;
    for d in 1..=lim.max(1) {
        if n % d == 0 {
            out.push(d);
            if d != n / d {
                out.push(n / d);
            }
        }
    }
    out.sort_unstable();
    out
}

fn rational_root_candidates(lead: i64, constant: i64) -> Vec<Rational> {
    let ps = divisors_pos(constant);
    let qs = divisors_pos(lead);
    let mut out: Vec<Rational> = Vec::new();
    for &p in &ps {
        for &q in &qs {
            for s in [1i64, -1] {
                if let Ok(r) = Rational::new(s * p, q) {
                    if !out.iter().any(|x| (x.to_f64() - r.to_f64()).abs() < 1e-15) {
                        out.push(r);
                    }
                }
            }
        }
    }
    out.sort_by(|a, b| a.to_f64().partial_cmp(&b.to_f64()).unwrap());
    out
}

fn eval_asc_rat(c: &[Rational], x: &Rational) -> Rational {
    let mut acc = Rational::from_int(0);
    let mut pow = Rational::from_int(1);
    for coeff in c {
        acc = acc + *coeff * pow;
        pow = pow * *x;
    }
    acc
}

fn eval_asc_f64(c: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    let mut pow = 1.0;
    for &coeff in c {
        acc += coeff * pow;
        pow *= x;
    }
    acc
}

/// Synthetic division on descending coefficients; divisor `(x - root)`.
fn synthetic_divide_desc(desc: &[Rational], root: Rational) -> (Vec<Rational>, Rational) {
    let mut b = desc[0];
    let mut qs = vec![b];
    for i in 1..desc.len() {
        let next = desc[i] + root * b;
        if i < desc.len() - 1 {
            b = next;
            qs.push(b);
        } else {
            return (qs, next);
        }
    }
    (qs, Rational::from_int(0))
}

fn deflate_rational_asc(asc: &[Rational], root: Rational) -> Vec<Rational> {
    let desc: Vec<Rational> = asc.iter().rev().cloned().collect();
    let (q_desc, _rem) = synthetic_divide_desc(&desc, root);
    q_desc.into_iter().rev().collect()
}

fn asc_to_rational_int(asc: &[f64]) -> Option<Vec<Rational>> {
    asc.iter()
        .map(|&v| Some(Rational::from_int(near_int(v)?)))
        .collect()
}

fn solve_quadratic_f64(a: f64, b: f64, c: f64) -> Result<Vec<Expr>> {
    if a != 0.0 {
        if let Some(roots) = try_exact_rational_quadratic_roots(a, b, c) {
            return Ok(roots);
        }
        let d = b * b - 4.0 * a * c;
        let d_scale = (b * b).abs().max((4.0 * a * c).abs()).max(1.0);
        let (r1, r2) = if d.abs() < 1e-12 * d_scale {
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
        Ok(vec![num(r1), num(r2)])
    } else if b != 0.0 {
        Ok(vec![num(-c / b)])
    } else if c == 0.0 {
        Ok(vec![num(0.0)])
    } else {
        Err(MathError::InvalidArgument(
            "equation has no variable dependence (constant ≠ 0)".into(),
        ))
    }
}

fn try_exact_rational_deflation(asc: &[f64]) -> Option<Vec<Expr>> {
    let mut work = asc_to_rational_int(asc)?;
    let mut roots: Vec<Expr> = Vec::new();
    loop {
        while work.len() > 1 && work.last().map(|t| t.num()) == Some(0) {
            work.pop();
        }
        let deg = work.len().saturating_sub(1);
        if deg == 0 {
            return Some(roots);
        }
        if deg == 1 {
            let c0 = work[0];
            let c1 = work[1];
            if c1.num() == 0 {
                return if c0.num() == 0 { Some(roots) } else { None };
            }
            roots.push(expr_from_rational(-c0 / c1));
            return Some(roots);
        }
        if deg == 2 {
            let mut quad = solve_quadratic_f64(
                work[2].to_f64(),
                work[1].to_f64(),
                work[0].to_f64(),
            )
            .ok()?;
            roots.append(&mut quad);
            return Some(roots);
        }
        let lead = work[work.len() - 1].num();
        let constant = work[0].num();
        if lead == 0 {
            continue;
        }
        let cands = rational_root_candidates(lead, constant);
        let r = cands
            .into_iter()
            .find(|cand| eval_asc_rat(&work, cand).num() == 0)?;
        loop {
            roots.push(expr_from_rational(r));
            work = deflate_rational_asc(&work, r);
            while work.len() > 1 && work.last().map(|t| t.num()) == Some(0) {
                work.pop();
            }
            if work.len() <= 1 {
                return Some(roots);
            }
            if eval_asc_rat(&work, &r).num() != 0 {
                break;
            }
        }
    }
}

fn solve_numeric_asc(asc: &[f64]) -> Result<Vec<Expr>> {
    let mut desc: Vec<f64> = asc.iter().rev().cloned().collect();
    while desc.len() > 1 && desc[0].abs() == 0.0 {
        desc.remove(0);
    }
    if desc.len() <= 1 {
        return Err(MathError::InvalidArgument(
            "polynomial has no variable dependence".into(),
        ));
    }
    let z = polynomial_roots_complex(&desc)?;
    Ok(z
        .into_iter()
        .map(|(re, im)| {
            if im.abs() < 1e-7 {
                num(re)
            } else {
                complex_root_expr(re, im)
            }
        })
        .collect())
}

fn certify_roots_f64(asc: &[f64], roots: &[Expr]) -> Result<()> {
    let scale = asc.iter().map(|c| c.abs()).sum::<f64>().max(1.0);
    for r in roots {
        if let Some(x) = root_to_f64(r) {
            let val = eval_asc_f64(asc, x);
            if val.abs() > 1e-4 * scale {
                return Err(MathError::Eval(format!(
                    "qsolve: numeric root {r} failed certification (p(x) ≈ {val})"
                )));
            }
        }
    }
    Ok(())
}

fn root_to_f64(e: &Expr) -> Option<f64> {
    use Expr::*;
    match e {
        Num(v) => Some(*v),
        Div(n, d) => {
            if let (Num(a), Num(b)) = (n.as_ref(), d.as_ref()) {
                if *b != 0.0 {
                    return Some(a / b);
                }
            }
            None
        }
        _ => None,
    }
}

fn solve_from_asc(asc: &[f64]) -> Result<Vec<Expr>> {
    let deg = poly_degree(asc);
    if deg == 0 {
        return Err(MathError::InvalidArgument(
            "equation has no variable dependence (constant = 0)".into(),
        ));
    }
    if deg > MAX_QSOLVE_DEG {
        return Err(MathError::InvalidArgument(format!(
            "degree > {MAX_QSOLVE_DEG} is not supported symbolically; use poly-roots for numeric roots"
        )));
    }
    if deg <= 2 {
        let a = *asc.get(2).unwrap_or(&0.0);
        let b = asc.get(1).copied().unwrap_or(0.0);
        let c = asc[0];
        return solve_quadratic_f64(a, b, c);
    }
    if let Some(roots) = try_exact_rational_deflation(asc) {
        return Ok(roots);
    }
    let roots = solve_numeric_asc(asc)?;
    certify_roots_f64(asc, &roots)?;
    Ok(roots)
}

/// Build `re + im*i` with sign-aware joins (never renders `+ -`).
fn complex_root_expr(re: f64, im: f64) -> Expr {
    let re_e = num(re);
    let im = Rational::snap_rational(im, crate::rational::SNAP_MAX_DEN_COMPACT);
    let imag = |mag: f64| -> Expr {
        let e = num(mag);
        if matches!(e, Expr::Num(v) if (v - 1.0).abs() < 1e-15) {
            Expr::var("i")
        } else {
            Expr::mul(e, Expr::var("i"))
        }
    };
    if im == 0.0 {
        return re_e;
    }
    if matches!(re_e, Expr::Num(v) if v == 0.0) {
        return if im < 0.0 { Expr::neg(imag(-im)) } else { imag(im) };
    }
    if im < 0.0 {
        Expr::sub(re_e, imag(-im))
    } else {
        Expr::add(re_e, imag(im))
    }
}

/// Solve `expr = 0` symbolically for `var` (polynomial degree ≤ 4).
pub fn solve_symbolic(expr: &Expr, var: &str) -> Result<Vec<Expr>> {
    let poly = to_poly(&crate::simplify::simplify(expr)).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "'{expr}' is not a polynomial in {var}; use the numeric solve instead"
        ))
    })?;
    let asc = asc_coeffs(&poly, var).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "degree > {MAX_QSOLVE_DEG} in {var} is not supported symbolically; use poly-roots for numeric roots"
        ))
    })?;
    solve_from_asc(&asc)
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
    let asc = asc_coeffs(&poly, &var).ok_or_else(|| {
        MathError::InvalidArgument(format!(
            "degree > {MAX_QSOLVE_DEG} in {var} is not supported symbolically; use poly-roots for numeric roots"
        ))
    })?;

    let mut steps = Vec::new();
    let eq_display = match rhs {
        Some(r) => format!("{lhs} = {r}"),
        None => format!("{lhs} = 0"),
    };
    steps.push(format!("qsolve: {eq_display}  (in {var})"));

    let deg = poly_degree(&asc);
    if deg == 0 {
        return Err(MathError::InvalidArgument(
            "equation has no variable dependence (constant = 0)".into(),
        ));
    }
    let fmt_coeff = |v: f64| -> String { format!("{}", num(v)) };
    if deg == 1 {
        let b = asc[1];
        let c = asc[0];
        steps.push(format!(
            "linear: {}*{} + {} = 0",
            fmt_coeff(b),
            var,
            fmt_coeff(c)
        ));
        let roots = solve_from_asc(&asc)?;
        steps.push(format!("{var} = {}", fmt_root(&roots[0])));
        return Ok(steps);
    }
    if deg == 2 {
        let a = asc[2];
        let b = asc[1];
        let c = asc[0];
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
    } else {
        steps.push(format!("degree {deg}: rational-root test (integer coefficients) + numeric roots if needed"));
    }
    let roots = solve_from_asc(&asc)?;
    let joined: Vec<String> = roots.iter().map(|r| format!("{var} = {}", fmt_root(r))).collect();
    steps.push(joined.join(", "));
    Ok(steps)
}

/// Render a root expression with rational display where possible.
/// `Expr::Num(v)` becomes `"5/3"` instead of `"1.6666666667"` when a
/// clean fraction fits within `max_den = 10_000`.  Complex roots built
/// by `complex_root_expr` get their real part and imaginary coefficient
/// rationalized too (so `x² + x + 1` prints `-1/2 ± 0.866…i`).  Falls
/// back to the default `Expr` Display for anything else.
///
/// Cross-port pattern from unsolve `solver.ts:399-406` (exact complex
/// real parts) and `limit.ts:537-562` (rational pole display).
fn fmt_root(e: &Expr) -> String {
    use Expr::*;
    match e {
        Num(v) => Rational::snap_to_string(*v, crate::rational::SNAP_MAX_DEN_COMPACT, None)
            .unwrap_or_else(|| format!("{}", v)),
        Div(num, den) => {
            if let (Num(n), Num(d)) = (num.as_ref(), den.as_ref()) {
                if *d != 0.0 {
                    if let (Some(ni), Some(di)) = (near_int(*n), near_int(*d)) {
                        if let Ok(r) = Rational::new(ni, di) {
                            return r.to_string();
                        }
                    }
                }
            }
            format!("{e}")
        }
        // ±im*i or ±re ± im*i built by complex_root_expr
        Add(re_box, im_box) | Sub(re_box, im_box) => {
            if let (Num(re), Mul(im_lhs, im_rhs)) = (re_box.as_ref(), im_box.as_ref()) {
                if let (Num(im), Var(v)) = (im_lhs.as_ref(), im_rhs.as_ref()) {
                    if v == "i" {
                        let re_s = crate::rational::Rational::snap_to_string(*re, crate::rational::SNAP_MAX_DEN_COMPACT, None)
                            .unwrap_or_else(|| format!("{}", re));
                        let im_s = crate::rational::Rational::snap_to_string(*im, crate::rational::SNAP_MAX_DEN_COMPACT, None)
                            .unwrap_or_else(|| format!("{}", im));
                        let sign = if matches!(e, Sub(..)) { " - " } else { " + " };
                        let im_part = if im_s == "1" { "i".to_string() } else { format!("{im_s}*i") };
                        return format!("{re_s}{sign}{im_part}");
                    }
                }
            }
            format!("{e}")
        }
        // ±im*i with no real part (Neg(Mul(Num(im), Var(i))))
        Neg(inner_box) => {
            if let Mul(im_lhs, im_rhs) = inner_box.as_ref() {
                if let (Num(im), Var(v)) = (im_lhs.as_ref(), im_rhs.as_ref()) {
                    if v == "i" {
                        let im_s = crate::rational::Rational::snap_to_string(-*im, crate::rational::SNAP_MAX_DEN_COMPACT, None)
                            .unwrap_or_else(|| format!("{}", -im));
                        let im_part = if im_s == "1" { "i".to_string() } else { format!("{im_s}*i") };
                        return format!("-{im_part}");
                    }
                }
            }
            format!("{e}")
        }
        // bare im*i
        Mul(im_lhs, im_rhs) => {
            if let (Num(im), Var(v)) = (im_lhs.as_ref(), im_rhs.as_ref()) {
                if v == "i" {
                    let im_s = crate::rational::Rational::snap_to_string(*im, crate::rational::SNAP_MAX_DEN_COMPACT, None)
                        .unwrap_or_else(|| format!("{}", im));
                    return if im_s == "1" { "i".to_string() } else { format!("{im_s}*i") };
                }
            }
            format!("{e}")
        }
        _ => format!("{e}"),
    }
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
        assert!(qsolve_steps_str("x^5 - 1").is_err());
        assert!(qsolve_steps_str("x*y = 1").is_err()); // two variables
        assert!(qsolve_steps_str("1/x = 2").is_err()); // variable denominator
    }

    #[test]
    fn cubic_integer_factor() {
        let r = roots("x^3 - 6*x^2 + 11*x - 6");
        assert!(r.contains(&"x = 1".to_string()), "{r:?}");
        assert!(r.contains(&"x = 2".to_string()), "{r:?}");
        assert!(r.contains(&"x = 3".to_string()), "{r:?}");
    }

    #[test]
    fn cubic_one_rational_root_then_quadratic() {
        let r = roots("x^3 - 1");
        assert!(r.iter().any(|s| s.contains("x = 1")), "{r:?}");
        assert_eq!(r.len(), 3);
    }

    #[test]
    fn api_level_solve_symbolic() {
        let e = Parser::parse("x^2 - 2x + 5").unwrap();
        let rs = solve_symbolic(&e, "x").unwrap();
        assert_eq!(rs.len(), 2);
        assert_eq!(rs[0].to_string(), "1 + 2*i");
    }

    // ---- Cross-port rational display (unsolve `solver.ts:380-451`) ----

    #[test]
    fn rational_display_linear_root() {
        // 2x + 7 = 0 → x = -7/2 (not -3.5)
        assert_eq!(roots("2*x + 7"), vec!["x = -7/2"]);
    }

    #[test]
    fn rational_display_quadratic_with_rational_roots() {
        // 6x² - 5x - 1 = 0 → roots are exactly 1 and -1/6.
        let r = roots("6*x^2 - 5*x - 1");
        assert!(r.contains(&"x = 1".to_string()), "got {r:?}");
        assert!(r.contains(&"x = -1/6".to_string()), "got {r:?}");
    }

    #[test]
    fn rational_display_complex_real_part() {
        // x² + x + 1 = 0 → re = -1/2, not -0.5
        let r = roots("x^2 + x + 1");
        for s in &r {
            assert!(s.contains("-1/2"), "expected exact -1/2 in: {s}");
            assert!(!s.contains("-0.5"), "expected no decimal -0.5 in: {s}");
        }
    }

    #[test]
    fn rational_display_simple_integer_roots() {
        assert_eq!(roots("x^2 - 5x + 6"), vec!["x = 3", "x = 2"]);
    }
}
