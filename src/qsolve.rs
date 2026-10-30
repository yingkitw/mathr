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

/// Read a coefficient as an exact rational, if it is one.
///
/// A coefficient is exact when it is an integer, or a *small* rational that
/// the value round-trips to (so `1/2` read back as `0.5` is recognised as
/// `1/2`, but `0.3333333333` is not mistaken for `1/3`).
fn exact_coeff(v: f64) -> Option<Rational> {
    if !v.is_finite() {
        return None;
    }
    if let Some(i) = near_int(v) {
        return Some(Rational::from_int(i));
    }
    // A fraction with a small denominator: verify by round-tripping.
    let r = Rational::from_float_snap(v, crate::rational::SNAP_MAX_DEN_COMPACT)?;
    if (r.to_f64() - v).abs() > 1e-12 * v.abs().max(1.0) {
        return None;
    }
    Some(r)
}

fn expr_from_rational(r: Rational) -> Expr {
    r.to_expr()
}

/// Exact quadratic roots for **rational** coefficients, expressed with exact
/// radicals when the discriminant is not a perfect square.
///
/// The roots of `a·x² + b·x + c` are `(-b ± √(b² − 4ac)) / 2a`. When the
/// coefficients are rationals, that whole expression can be kept exact:
///
/// - `(1/2)x² − (3/4)x + 1/8` → `(3 + sqrt(5))/4` and `(3 − sqrt(5))/4`
/// - `(2/3)x² + (1/3)x − 1` → `1` and `-3/2` (the discriminant is a
///   perfect square, so the radical disappears)
///
/// The radical is built with [`Expr::func`]("sqrt", …) over a *rational*
/// discriminant, so no floating-point value is ever introduced.
fn exact_rational_quadratic_radicals(a: f64, b: f64, c: f64) -> Option<Vec<Expr>> {
    let ar = exact_coeff(a)?;
    let br = exact_coeff(b)?;
    let cr = exact_coeff(c)?;
    if ar.num() == 0 {
        return None;
    }
    // Discriminant: b² − 4ac, computed exactly.
    let disc = br * br - Rational::from_int(4) * ar * cr;
    if disc.num() < 0 {
        // Negative discriminant — the roots are complex, not real radicals.
        return None;
    }
    let two_a = ar + ar;
    if two_a.num() == 0 {
        return None;
    }
    // A perfect-square discriminant gives plain rationals.
    if let Some(sd) = disc.exact_sqrt() {
        let r1 = (-br + sd) / two_a;
        let r2 = (-br - sd) / two_a;
        return Some(vec![expr_from_rational(r1), expr_from_rational(r2)]);
    }
    // Otherwise keep the radical symbolic: (-b ± √disc) / 2a.
    //
    // The radical is built in *simplified* form: the perfect-square factor
    // of the discriminant is pulled out, so `√8` becomes `2·√2` and `√(5/16)`
    // becomes `¼·√5`. Every coefficient stays a rational, so no
    // floating-point value is ever introduced.
    let (rad_coeff, rad_inner) = split_radical(disc);
    let sqrt_disc = Expr::func("sqrt", vec![expr_from_rational(rad_inner)]);
    Some(vec![
        build_radical_root(&br, &rad_coeff, &sqrt_disc, &two_a, false),
        build_radical_root(&br, &rad_coeff, &sqrt_disc, &two_a, true),
    ])
}

/// Build one root of the quadratic in tidied radical form.
///
/// The root is `(-b ± c·√d) / e` with `b`, `c`, `e` rational and `d`
/// square-free. The result is written as `p + q·√d` (or `q·√d`) with the
/// denominator already folded in, so `2·√2 / 2` reads `√2` and
/// `(3/4 + ¼·√5) / 1` reads `3/4 + ¼·√5`.
fn build_radical_root(
    b: &Rational,
    rad_coeff: &Rational,
    sqrt_disc: &Expr,
    den: &Rational,
    minus: bool,
) -> Expr {
    // The rational part: -b / e.
    let rat = (-*b) / *den;
    // The radical coefficient: ±c / e.
    let mut coef = *rad_coeff / *den;
    if minus {
        coef = -coef;
    }
    // A rational part of zero is dropped.
    let rat_part = if rat.num() == 0 {
        None
    } else {
        Some(expr_from_rational(rat))
    };
    // The radical term: coef·√d. A coefficient with a denominator is
    // written as a division (`√5/4` rather than `1/4*√5`), which reads
    // better and matches how a rational is displayed elsewhere.
    let radical_part = if coef.num() == 0 {
        None
    } else if coef.num() == 1 && coef.den() == 1 {
        Some(sqrt_disc.clone())
    } else if coef.num() == -1 && coef.den() == 1 {
        Some(Expr::neg(sqrt_disc.clone()))
    } else if coef.den() != 1 {
        // √d · num / den — keep the sign on the numerator.
        let signed = if coef.num() < 0 {
            Expr::neg(sqrt_disc.clone())
        } else {
            sqrt_disc.clone()
        };
        Some(Expr::div(signed, Expr::num(coef.den() as f64)))
    } else {
        Some(Expr::mul(Expr::num(coef.num() as f64), sqrt_disc.clone()))
    };
    match (rat_part, radical_part) {
        (None, None) => Expr::num(0.0),
        (Some(r), None) => r,
        (None, Some(t)) => t,
        (Some(r), Some(t)) => {
            // A negative radical term reads better as a subtraction:
            // `3/4 - sqrt(5)/4` rather than `3/4 + -sqrt(5)/4`.
            match &t {
                Expr::Neg(inner) => Expr::sub(r, (**inner).clone()),
                // `-√d / n` → `r - √d/n`.
                Expr::Div(n, d) => match n.as_ref() {
                    Expr::Neg(inner) => Expr::sub(r, Expr::div((**inner).clone(), (**d).clone())),
                    _ => Expr::add(r, t),
                },
                _ => Expr::add(r, t),
            }
        }
    }
}

/// Split a positive rational into `coefficient · √(inner)` with `inner`
/// square-free.
///
/// `8` → `(2, 2)`, `5/16` → `(1/4, 5)`, `12` → `(2, 3)`. The coefficient is
/// returned as a [`Rational`] so no floating-point value is introduced.
fn split_radical(r: Rational) -> (Rational, Rational) {
    // Pull the largest square *factor* out of the numerator and the
    // denominator. `8 = 4·2` → the square factor is 4, so √8 = 2·√2.
    let (n_root, inner_num) = split_square_factor(r.num().unsigned_abs());
    let (d_root, inner_den) = split_square_factor(r.den().unsigned_abs());
    let inner = Rational::new(inner_num as i64, inner_den as i64).unwrap_or(r);
    let coeff = Rational::new(n_root as i64, d_root as i64).unwrap_or_else(|_| Rational::from_int(1));
    (coeff, inner)
}

/// Split `n` into `root² · rest` with `rest` square-free.
///
/// `8` → `(2, 2)`, `12` → `(2, 3)`, `16` → `(4, 1)`, `5` → `(1, 5)`,
/// `18` → `(3, 2)`, `50` → `(5, 2)`.
fn split_square_factor(n: u64) -> (u64, u64) {
    if n == 0 {
        return (0, 0);
    }
    let mut rest = n;
    let mut root: u64 = 1;
    // Trial division by small primes is enough to find the square factor.
    let mut p: u64 = 2;
    while p * p <= rest {
        let mut count = 0u32;
        while rest.is_multiple_of(p) {
            rest /= p;
            count += 1;
        }
        // Every *pair* of the prime contributes to the square root; a lone
        // leftover one stays in the square-free part.
        root *= p.pow(count / 2);
        if count % 2 == 1 {
            rest *= p;
        }
        p += if p == 2 { 1 } else { 2 };
    }
    // Whatever remains is either 1 or a prime (already square-free).
    (root, rest)
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
        // Rational coefficients: exact rationals when the discriminant is a
        // perfect square, exact radicals otherwise.
        if let Some(roots) = exact_rational_quadratic_radicals(a, b, c) {
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
    fn irrational_discriminant_gives_exact_radicals() {
        // x^2 - 2 has roots ±√2 — kept exact rather than collapsed to
        // decimals.
        let r = roots("x^2 - 2");
        assert_eq!(r[0], "x = sqrt(2)", "got {r:?}");
        assert_eq!(r[1], "x = -sqrt(2)", "got {r:?}");
        // The value must still be right.
        let e = crate::parser::Parser::parse("sqrt(2)").unwrap();
        let v = crate::eval::eval(&e, &crate::eval::Context::standard()).unwrap();
        assert!((v - std::f64::consts::SQRT_2).abs() < 1e-10, "got {v}");
    }

    // ---- exact radical roots ------------------------------------------------

    #[test]
    fn radical_pulls_out_square_factors() {
        // √8 = 2√2, √12 = 2√3, √50 = 5√2.
        assert_eq!(roots("x^2 - 8")[0], "x = 2*sqrt(2)", "got {:?}", roots("x^2 - 8"));
        assert_eq!(roots("x^2 - 12")[0], "x = 2*sqrt(3)", "got {:?}", roots("x^2 - 12"));
        assert_eq!(roots("x^2 - 50")[0], "x = 5*sqrt(2)", "got {:?}", roots("x^2 - 50"));
    }

    #[test]
    fn radical_roots_are_signed() {
        let r = roots("x^2 - 8");
        assert_eq!(r[0], "x = 2*sqrt(2)", "got {r:?}");
        assert_eq!(r[1], "x = -2*sqrt(2)", "got {r:?}");
    }

    #[test]
    fn rational_coefficients_give_exact_radicals() {
        // (1/2)x² − (3/4)x + 1/8: discriminant 5/16, roots (3 ± √5)/4.
        let r = roots("(1/2)x^2-(3/4)x+1/8");
        assert_eq!(r[0], "x = 3/4 + sqrt(5)/4", "got {r:?}");
        assert_eq!(r[1], "x = 3/4 - sqrt(5)/4", "got {r:?}");
    }

    #[test]
    fn rational_coefficients_with_square_discriminant_give_rationals() {
        // (2/3)x² + (1/3)x − 1 has discriminant 1/9 → rational roots.
        let r = roots("(2/3)x^2+(1/3)x-1");
        assert_eq!(r[0], "x = 1", "got {r:?}");
        assert_eq!(r[1], "x = -3/2", "got {r:?}");
    }

    #[test]
    fn radical_roots_evaluate_correctly() {
        // The numeric value must match the true root.
        for (expr, want) in [("x^2-2", 2.0f64.sqrt()), ("x^2-8", 8.0f64.sqrt())] {
            let r = roots(expr);
            let e = crate::parser::Parser::parse(r[0].strip_prefix("x = ").unwrap()).unwrap();
            let v = crate::eval::eval(&e, &crate::eval::Context::standard()).unwrap();
            assert!((v - want).abs() < 1e-12, "{}: got {} want {}", expr, v, want);
        }
    }

    #[test]
    fn radical_roots_satisfy_the_polynomial() {
        // Substituting a radical root back into the equation gives ~0.
        for expr in ["x^2-2", "x^2-8", "(1/2)x^2-(3/4)x+1/8"] {
            let r = roots(expr);
            let e = crate::parser::Parser::parse(expr).unwrap();
            for root in &r {
                let root_expr =
                    crate::parser::Parser::parse(root.strip_prefix("x = ").unwrap()).unwrap();
                let mut ctx = crate::eval::Context::standard();
                ctx.set("x", 0.0);
                // Evaluate the polynomial with the root substituted.
                let val = crate::eval::eval(&e.substitute("x", &root_expr), &ctx).unwrap();
                assert!(val.abs() < 1e-9, "{} root {} gives {}", expr, root, val);
            }
        }
    }

    #[test]
    fn integer_quadratic_roots_unchanged() {
        // The pre-existing exact-integer path must not regress.
        let r = roots("6x^2-5x-1");
        assert_eq!(r[0], "x = 1", "got {r:?}");
        assert_eq!(r[1], "x = -1/6", "got {r:?}");
    }

    #[test]
    fn complex_roots_still_use_floats() {
        // A negative discriminant has no real radical — the complex path
        // is unchanged.
        let r = roots("x^2+x+1");
        assert!(r[0].contains('i'), "got {r:?}");
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
