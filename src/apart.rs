//! Partial fraction decomposition ("apart") for rational functions.
//!
//! Decomposes `N(x)/D(x)` (univariate, polynomial denominator) into
//!
//! ```text
//! N/D = Q + Σ  cᵢₖ / Fᵢ^k
//! ```
//!
//! where `Q` is the polynomial quotient from long division and the `Fᵢ`
//! are the linear `(x − r)` and irreducible quadratic
//! `(x² + p·x + q)` factors of `D` (roots found numerically via
//! Durand–Kerner, conjugate complex roots paired into quadratics). The
//! unknown coefficients form a square linear system (multiply through by
//! `D` and equate coefficients) solved with Gaussian elimination.
//!
//! This is a numeric CAS decomposition: coefficients are `f64` and
//! repeated roots carry the usual numerical fuzz of multiple-root
//! finding. Coefficients near integers are snapped for display.
//!
//! # Example
//!
//! ```
//! use mathr::apart::apart;
//! use mathr::parser::Parser;
//!
//! let e = Parser::parse("1/(x*(x+1))").unwrap();
//! let out = apart(&e, "x").unwrap().to_string();
//! // 1/(x(x+1)) = 1/x − 1/(x+1)
//! assert!(out.contains("1/x"), "got: {}", out);
//! assert!(out.contains("1/(x + 1)"), "got: {}", out);
//! ```

use crate::error::{MathError, Result};
use crate::expr::Expr::{self, *};
use crate::poly::{to_poly, poly_to_expr};
use crate::simplify::simplify;

/// Safety cap on denominator degree.
const MAX_DEGREE: usize = 32;

/// One denominator factor.
#[derive(Debug, Clone)]
enum Factor {
    /// `(x − root)^mult`
    Linear { root: f64, mult: u32 },
    /// `(x² + p·x + q)^mult` with conjugate-complex roots
    Quadratic { p: f64, q: f64, mult: u32 },
}

/// Decompose `expr` (a quotient of univariate polynomials in `var`) into
/// partial fractions. Returns the simplified result.
pub fn apart(expr: &Expr, var: &str) -> Result<Expr> {
    let mut steps = Vec::new();
    let result = apart_inner(expr, var, &mut steps)?;
    Ok(simplify(&result))
}

/// Decompose with human-readable intermediate steps.
pub fn apart_steps(expr: &Expr, var: &str) -> Result<Vec<String>> {
    let mut steps = vec![format!("apart: {}", simplify(expr))];
    let result = apart_inner(expr, var, &mut steps)?;
    steps.push(format!("= {}", simplify(&result)));
    Ok(steps)
}

/// Convenience wrapper for the REPL step dispatch: parse `"<expr> [var]"`
/// (variable inferred when omitted) and return decomposition steps.
pub fn apart_steps_str(src: &str) -> Result<Vec<String>> {
    let tokens: Vec<&str> = src.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(MathError::Eval("apart needs: <expr> [<var>]".into()));
    }
    let (e, var) = if tokens.len() >= 2 {
        let var = tokens[tokens.len() - 1].to_string();
        let expr_src = tokens[..tokens.len() - 1].join(" ");
        (crate::parser::Parser::parse(&expr_src)?, var)
    } else {
        let e = crate::parser::Parser::parse(tokens[0])?;
        let mut vars = e.variables();
        if vars.len() != 1 {
            return Err(MathError::Eval(
                "apart needs: <expr> <var> (expression has multiple or no variables)".into(),
            ));
        }
        (e, vars.remove(0))
    };
    apart_steps(&e, &var)
}

fn apart_inner(expr: &Expr, var: &str, steps: &mut Vec<String>) -> Result<Expr> {
    let s = simplify(expr);

    // The expression must be a quotient; a bare polynomial is already
    // "decomposed".
    let Expr::Div(num, den) = &s else {
        if to_poly(&s).is_some() {
            return Ok(s);
        }
        return Err(MathError::Eval(
            "apart needs a quotient of polynomials, e.g. 1/(x^2 - 1)".into(),
        ));
    };

    let num_sparse = to_poly(num).ok_or_else(|| {
        MathError::Eval("apart: numerator must be a polynomial".into())
    })?;
    let den_sparse = to_poly(den).ok_or_else(|| {
        MathError::Eval("apart: denominator must be a polynomial".into())
    })?;

    let n = dense(&num_sparse, var)?;
    let d = dense(&den_sparse, var)?;

    if d.len() < 2 {
        // Constant denominator: nothing to decompose.
        return Ok(s);
    }
    if d.len() - 1 > MAX_DEGREE {
        return Err(MathError::Eval(format!(
            "apart: denominator degree {} exceeds the maximum of {}",
            d.len() - 1,
            MAX_DEGREE
        )));
    }

    // 1. Polynomial long division: N = Q·D + R with deg(R) < deg(D).
    let (q, r) = divrem(&n, &d);
    if !q.is_empty() {
        steps.push(format!(
            "long division: {} = {}·({}) + {}",
            poly_str(&n, var),
            poly_str(&q, var),
            poly_str(&d, var),
            poly_str(&r, var)
        ));
    }

    // 2. Factor the denominator from its numeric roots.
    let factors = factor_denominator(&d)?;
    steps.push(format!(
        "denominator {} = {}",
        poly_str(&d, var),
        factors
            .iter()
            .map(factor_str)
            .collect::<Vec<_>>()
            .join(" · ")
    ));

    // 3. Build the linear system: R = Σ columns·coefficients.
    // One column per unknown — a linear factor power k contributes
    // D/(x−r)^k (constant numerator); a quadratic factor power k
    // contributes TWO columns, x·D/F^k and D/F^k (numerator B·x + C).
    #[derive(Debug, Clone, Copy)]
    enum ColKind {
        Linear,
        QuadX,
        QuadConst,
    }
    let tol = d.iter().fold(1.0f64, |m, &c| m.max(c.abs())) * 1e-6;
    let mut columns: Vec<(Vec<f64>, ColKind, usize)> = Vec::new();
    for (fi, f) in factors.iter().enumerate() {
        for k in 1..=f.mult() {
            let fk = factor_pow(f, k);
            let (mut col, rem) = divrem(&d, &fk);
            if rem.iter().any(|c| c.abs() > tol) {
                return Err(MathError::Eval(
                    "apart: denominator factoring is inconsistent (numerical)".into(),
                ));
            }
            match f {
                Factor::Linear { .. } => columns.push((col, ColKind::Linear, fi)),
                Factor::Quadratic { .. } => {
                    // Multiplying by x in descending representation appends 0.
                    let quad_const = col.clone();
                    col.push(0.0);
                    columns.push((col, ColKind::QuadX, fi));
                    columns.push((quad_const, ColKind::QuadConst, fi));
                }
            }
        }
    }
    let k_dim = columns.len();
    let deg_r = d.len() - 1;
    if k_dim != deg_r {
        return Err(MathError::Eval(format!(
            "apart: factor degrees do not sum to the denominator degree ({} vs {})",
            k_dim,
            deg_r
        )));
    }

    let mut rhs = vec![0.0f64; k_dim];
    // r is stored descending; reverse into ascending order.
    for (i, c) in r.iter().rev().enumerate() {
        rhs[i] = *c;
    }

    let mut data = vec![0.0f64; k_dim * k_dim];
    for (col, c) in columns.iter().enumerate() {
        for (row, &v) in c.0.iter().rev().enumerate() {
            data[row * k_dim + col] = v;
        }
    }
    let m = crate::matrix::Matrix::from_row_major(k_dim, k_dim, data)?;
    let coeffs = m.solve(&rhs)?;

    // 4. Assemble: Q + Σ numerator / factor^k, grouping the two quadratic
    // columns of each power into one B·x + C numerator.
    let mut terms: Vec<Expr> = Vec::new();
    if !q.is_empty() {
        let q_sparse = sparse(&q, var);
        terms.push(poly_to_expr(&q_sparse));
    }
    // Walk the columns in order, consuming one Linear column per power of
    // a linear factor and the QuadX+QuadConst pair per power of a
    // quadratic factor. Terms are joined with +/− by leading sign so
    // output reads `-1/(x + 1) + 1/x` rather than `+ -…`.
    let mut groups: Vec<(usize, u32, Vec<usize>)> = Vec::new();
    {
        let mut next = 0usize;
        for (fi, f) in factors.iter().enumerate() {
            for k in 1..=f.mult() {
                match f {
                    Factor::Linear { .. } => {
                        groups.push((fi, k, vec![next]));
                        next += 1;
                    }
                    Factor::Quadratic { .. } => {
                        groups.push((fi, k, vec![next, next + 1]));
                        next += 2;
                    }
                }
            }
        }
    }
    let mut out: Option<Expr> = None;
    // Polynomial quotient comes first.
    if !q.is_empty() {
        let q_sparse = sparse(&q, var);
        out = Some(poly_to_expr(&q_sparse));
    }
    for (fi, k, ids) in &groups {
        let f = &factors[*fi];
        let den_e = factor_expr(f, *k, var);
        let (num_mag, neg) = match f {
            Factor::Linear { .. } => {
                let c = snap(coeffs[ids[0]]);
                if c == 0.0 {
                    continue;
                }
                (Num(c.abs()), c < 0.0)
            }
            Factor::Quadratic { .. } => {
                let b = snap(coeffs[ids[0]]);
                let c = snap(coeffs[ids[1]]);
                if b == 0.0 && c == 0.0 {
                    continue;
                }
                (quad_num_mag(b, c, var), b < 0.0)
            }
        };
        let term = Div(Box::new(num_mag), Box::new(den_e));
        out = Some(match out {
            None => {
                if neg {
                    Neg(Box::new(term))
                } else {
                    term
                }
            }
            Some(acc) => {
                if neg {
                    Sub(Box::new(acc), Box::new(term))
                } else {
                    Add(Box::new(acc), Box::new(term))
                }
            }
        });
    }
    let out = out.unwrap_or(Num(0.0));
    Ok(out)
}

/// Magnitude form of a quadratic numerator `B·x + C` (leading coefficient
/// positive). Returns the caller the sign separately.
fn quad_num_mag(b: f64, c: f64, var: &str) -> Expr {
    let bx = |mag: f64| Mul(Box::new(Num(mag)), Box::new(Var(var.to_string())));
    if b.abs() < 1e-12 {
        return Num(c.abs());
    }
    let lhs = bx(b.abs());
    if c.abs() < 1e-12 {
        return lhs;
    }
    if (b < 0.0) == (c < 0.0) {
        Add(Box::new(lhs), Box::new(Num(c.abs())))
    } else {
        Sub(Box::new(lhs), Box::new(Num(c.abs())))
    }
}

/// Dense descending coefficient vector from a sparse polynomial.
fn dense(p: &crate::poly::Poly, var: &str) -> Result<Vec<f64>> {
    let mut max_deg: usize = 0;
    for t in p {
        for (name, e) in &t.monomial {
            if name != var {
                return Err(MathError::Eval(format!(
                    "apart: expression must be univariate in `{}` (found `{}`)",
                    var, name
                )));
            }
            max_deg = max_deg.max(*e as usize);
        }
    }
    let mut coeffs = vec![0.0f64; max_deg + 1];
    for t in p {
        let e = t.monomial.get(var).copied().unwrap_or(0) as usize;
        coeffs[e] += t.coeff;
    }
    while coeffs.len() > 1 && coeffs.last() == Some(&0.0) {
        coeffs.pop();
    }
    coeffs.reverse(); // descending
    Ok(coeffs)
}

/// Sparse ascending polynomial from a dense descending vector.
fn sparse(d: &[f64], var: &str) -> crate::poly::Poly {
    let n = d.len();
    d.iter()
        .enumerate()
        .filter(|(_, &c)| c != 0.0)
        .map(|(i, &c)| {
            let exp = (n - 1 - i) as u32;
            crate::poly::Term {
                coeff: c,
                monomial: if exp == 0 {
                    Default::default()
                } else {
                    [(var.to_string(), exp)].into_iter().collect()
                },
            }
        })
        .collect()
}

/// Polynomial long division: returns (quotient, remainder), both descending.
fn divrem(a: &[f64], b: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut r: Vec<f64> = a.to_vec();
    let db = b.len() - 1;
    if r.len() < b.len() {
        return (Vec::new(), r);
    }
    let dq = r.len() - db - 1;
    let mut q = vec![0.0f64; dq + 1];
    for k in 0..=dq {
        let lead = r[k];
        if b[0] == 0.0 {
            break;
        }
        q[k] = lead / b[0];
        for j in 0..=db {
            r[k + j] -= q[k] * b[j];
        }
    }
    let rest = r[dq + 1..].to_vec();
    let mut rem = rest;
    while rem.len() > 1 && rem[0].abs() < 1e-12 {
        rem.remove(0);
    }
    (q, rem)
}

/// Factor a denominator polynomial into linear and irreducible quadratic
/// factors using numeric roots, clustering repeated roots.
fn factor_denominator(d: &[f64]) -> Result<Vec<Factor>> {
    let roots = crate::solver::polynomial_roots_complex(d)?;
    let deg = d.len() - 1;

    let mut reals: Vec<f64> = Vec::new();
    let mut complex: Vec<(f64, f64)> = Vec::new(); // positive-imag only
    for &(re, im) in &roots {
        let scale = re.abs().max(1.0);
        if im.abs() < 1e-6 * scale {
            reals.push(re);
        } else if im > 0.0 {
            complex.push((re, im));
        }
    }

    // Cluster real roots (multiplicity = count within tolerance).
    reals.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut factors: Vec<Factor> = Vec::new();
    let mut i = 0;
    while i < reals.len() {
        let mut j = i;
        while j + 1 < reals.len() && (reals[j + 1] - reals[i]).abs() <= 1e-5 * reals[i].abs().max(1.0)
        {
            j += 1;
        }
        let avg = reals[i..=j].iter().sum::<f64>() / (j - i + 1) as f64;
        factors.push(Factor::Linear {
            root: snap(avg),
            mult: (j - i + 1) as u32,
        });
        i = j + 1;
    }

    // Pair complex roots into quadratics, clustering conjugates with the
    // same real part and magnitude.
    complex.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then(a.1.partial_cmp(&b.1).unwrap()));
    let mut i = 0;
    while i < complex.len() {
        let (re, im) = complex[i];
        let mut j = i;
        while j + 1 < complex.len()
            && (complex[j + 1].0 - re).abs() <= 1e-5 * re.abs().max(1.0)
            && (complex[j + 1].1 - im).abs() <= 1e-5 * im.abs().max(1.0)
        {
            j += 1;
        }
        let count = j - i + 1;
        // Each occurrence in the list represents a conjugate pair.
        let p = -2.0 * snap(re);
        let q = snap(re * re + im * im);
        factors.push(Factor::Quadratic { p, q, mult: count as u32 });
        i = j + 1;
    }

    // Sanity: Σ multiplicities (quadratics count twice — conjugate pair)
    // must equal the degree.
    let total: usize = factors
        .iter()
        .map(|f| match f {
            Factor::Linear { mult, .. } => *mult as usize,
            Factor::Quadratic { mult, .. } => 2 * *mult as usize,
        })
        .sum();
    if total != deg {
        return Err(MathError::Eval(format!(
            "apart: could not factor the denominator into linear/quadratic factors (covered {} of {} degrees)",
            total, deg
        )));
    }
    Ok(factors)
}

impl Factor {
    fn mult(&self) -> u32 {
        match self {
            Factor::Linear { mult, .. } | Factor::Quadratic { mult, .. } => *mult,
        }
    }
}

fn factor_str(f: &Factor) -> String {
    match f {
        Factor::Linear { root, mult } => {
            let base = format!("(x - {})", root)
                .replace("- -", "+ ")
                .replace("(x - 0)", "x");
            if *mult == 1 {
                base
            } else {
                format!("{}^{}", base, mult)
            }
        }
        Factor::Quadratic { p, q, mult } => {
            let base = format!("(x^2 + {}*x + {})", p, q);
            if *mult == 1 {
                base
            } else {
                format!("{}^{}", base, mult)
            }
        }
    }
}

/// The factor raised to power `k`, as dense descending coefficients.
fn factor_pow(f: &Factor, k: u32) -> Vec<f64> {
    let base: Vec<f64> = match f {
        Factor::Linear { root, .. } => vec![1.0, -root],
        Factor::Quadratic { p, q, .. } => vec![1.0, *p, *q],
    };
    let mut out = vec![1.0f64];
    for _ in 0..k {
        out = poly_dense_mul(&out, &base);
    }
    out
}

fn poly_dense_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

/// The factor expression raised to power `k` (denominator of one term),
/// built with sign-aware joins so `x + 1` renders instead of `x - -1`.
fn factor_expr(f: &Factor, k: u32, var: &str) -> Expr {
    let x = || Var(var.to_string());
    let base = match f {
        Factor::Linear { root, .. } => {
            if root.abs() < 1e-12 {
                x()
            } else if *root < 0.0 {
                Add(Box::new(x()), Box::new(Num(-root)))
            } else {
                Sub(Box::new(x()), Box::new(Num(*root)))
            }
        }
        Factor::Quadratic { p, q, .. } => {
            let mut acc = Pow(Box::new(x()), Box::new(Num(2.0)));
            if p.abs() >= 1e-12 {
                acc = if *p < 0.0 {
                    Sub(Box::new(acc), Box::new(Mul(Box::new(Num(-p)), Box::new(x()))))
                } else {
                    Add(Box::new(acc), Box::new(Mul(Box::new(Num(*p)), Box::new(x()))))
                };
            }
            if q.abs() >= 1e-12 {
                acc = if *q < 0.0 {
                    Sub(Box::new(acc), Box::new(Num(-q)))
                } else {
                    Add(Box::new(acc), Box::new(Num(*q)))
                };
            }
            acc
        }
    };
    if k == 1 {
        base
    } else {
        Pow(Box::new(base), Box::new(Num(k as f64)))
    }
}

/// Snap to integers within 1e-9 relative, else round to 12 significant
/// digits (mirrors limit.rs snapping).
fn snap(v: f64) -> f64 {
    if !v.is_finite() {
        return v;
    }
    if (v - v.round()).abs() < 1e-9 * v.abs().max(1.0) {
        return v.round();
    }
    let digits = 11 - v.abs().log10().floor() as i32;
    if !(-20..=20).contains(&digits) {
        return v;
    }
    let factor = 10f64.powi(digits);
    (v * factor).round() / factor
}

/// Human-readable dense polynomial (descending powers).
fn poly_str(d: &[f64], var: &str) -> String {
    let p = sparse(d, var);
    poly_to_expr(&p).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    /// Check the decomposition numerically: original and result agree at
    /// sample points away from the poles.
    fn agrees(src: &str, var: &str, points: &[f64]) -> bool {
        let e = Parser::parse(src).unwrap();
        let out = apart(&e, var).unwrap();
        for &x in points {
            let mut ctx = crate::eval::Context::standard();
            ctx.set(var, x);
            let a = crate::eval::eval(&e, &ctx).unwrap();
            let b = crate::eval::eval(&out, &ctx).unwrap();
            if (a - b).abs() > 1e-6 * a.abs().max(1.0) {
                return false;
            }
        }
        true
    }

    #[test]
    fn apart_simple_linear_factors() {
        // 1/(x²−1) = 1/2·1/(x−1) − 1/2·1/(x+1)
        let e = Parser::parse("1/(x^2 - 1)").unwrap();
        let out = apart(&e, "x").unwrap();
        let s = out.to_string();
        assert!(s.contains("1/2") || s.contains("0.5"), "got: {}", s);
        assert!(agrees("1/(x^2 - 1)", "x", &[0.3, 0.9, 2.2, -3.1]));
    }

    #[test]
    fn apart_repeated_linear_factor() {
        assert!(agrees("1/(x*(x+1)^2)", "x", &[0.5, 2.0, -3.0, -0.4]));
        assert!(agrees("x/(x+1)^2", "x", &[1.0, 3.0, -0.5]));
    }

    #[test]
    fn apart_improper_fraction_long_division() {
        // (x²+1)/(x−1) = x + 1 + 2/(x−1)
        let e = Parser::parse("(x^2 + 1)/(x - 1)").unwrap();
        let s = apart(&e, "x").unwrap().to_string();
        assert!(s.contains("x"), "got: {}", s);
        assert!(s.contains("2"), "got: {}", s);
        assert!(agrees("(x^2 + 1)/(x - 1)", "x", &[0.5, 2.0, -1.5]));
    }

    #[test]
    fn apart_irreducible_quadratic() {
        assert!(agrees("1/(x^2 + 1)", "x", &[0.5, 2.0, -1.0]));
        assert!(agrees("(2*x + 1)/(x^2 + 1)", "x", &[0.5, 2.0, -1.0]));
        // quadratic + linear mix
        assert!(agrees("(3*x + 5)/((x^2 + x + 1)*(x - 2))", "x", &[0.5, 1.0, -2.0, 3.0]));
    }

    #[test]
    fn apart_trivial_polynomial_input() {
        let e = Parser::parse("x^2 + 1").unwrap();
        let out = apart(&e, "x").unwrap();
        assert!(out.equals(&Parser::parse("x^2 + 1").unwrap()));
    }

    #[test]
    fn apart_shifted_linear_factor() {
        assert!(agrees("1/(x - 3)", "x", &[0.5, 1.0, 2.9, 5.0]));
        assert!(agrees("2/(x + 4)", "x", &[0.5, 1.0, -3.9, 5.0]));
    }

    #[test]
    fn apart_multivariate_errors() {
        let e = Parser::parse("1/(x*y - 1)").unwrap();
        assert!(apart(&e, "x").is_err());
    }

    #[test]
    fn apart_non_quotient_with_functions_errors() {
        let e = Parser::parse("sin(x)").unwrap();
        assert!(apart(&e, "x").is_err());
    }

    #[test]
    fn apart_steps_structure() {
        let e = Parser::parse("1/(x*(x+1))").unwrap();
        let steps = apart_steps(&e, "x").unwrap();
        assert!(steps[0].contains("apart"), "got: {:?}", steps);
        assert!(steps.last().unwrap().contains("1/x"), "got: {:?}", steps);
    }

    #[test]
    fn snap_matches_limit_convention() {
        assert_eq!(snap(0.9999999999999983), 1.0);
        assert_eq!(snap(2.0000000000000004), 2.0);
        assert_eq!(snap(0.5), 0.5);
    }

    // ===== Additional coverage =====

    #[test]
    fn apart_repeated_quadratic() {
        // 1/(x²+1)² — irreducible quadratic with multiplicity 2
        let src = "1/(x^2 + 1)^2";
        assert!(agrees(src, "x", &[0.0, 0.5, 1.5, 2.0, -1.0, -0.3]));
    }

    #[test]
    fn apart_mixed_linear_and_quadratic() {
        // x/((x-1)*(x^2+1))
        let src = "x/((x - 1)*(x^2 + 1))";
        assert!(agrees(src, "x", &[0.0, 2.0, 3.0, -1.0, 0.5, -2.0]));
    }

    #[test]
    fn apart_non_monic_leading_coefficient() {
        // (3x² + 2)/(2x - 1) — non-monic denominator, improper fraction
        let src = "(3*x^2 + 2)/(2*x - 1)";
        assert!(agrees(src, "x", &[0.0, 3.0, 5.0, -1.0, 2.0, 10.0]));
    }

    #[test]
    fn apart_constant_numerator() {
        // 5/(x^2 - 9) = 5/((x-3)(x+3))
        let src = "5/(x^2 - 9)";
        assert!(agrees(src, "x", &[0.0, 1.0, 2.0, -1.0, -2.0, 4.0]));
    }

    #[test]
    fn apart_high_degree() {
        // (x + 1)/(x^3 - 1) = (x+1)/((x-1)(x²+x+1))
        // Cubic denominator with linear and irreducible quadratic factors
        let src = "(x + 1)/(x^3 - 1)";
        assert!(agrees(src, "x", &[2.0, 3.0, -2.0, 0.5, 4.0, -3.0]));
    }

    #[test]
    fn apart_negative_coefficients() {
        // (x - 3)/(x^2 + 2x - 3) = (x-3)/((x+3)(x-1))
        let src = "(x - 3)/(x^2 + 2*x - 3)";
        assert!(agrees(src, "x", &[0.0, 2.0, 4.0, -2.0, 5.0, -4.0]));
    }

    #[test]
    fn apart_steps_str_returns_string() {
        let steps = apart_steps_str("1/(x^2 - 1) x").unwrap();
        let combined = steps.join(" ");
        assert!(combined.contains("apart"), "output should mention apart: {}", combined);
    }

    #[test]
    fn apart_steps_str_empty_input_errors() {
        // Empty input should error
        assert!(apart_steps_str("").is_err());
    }

    #[test]
    fn apart_output_is_simplifiable() {
        // The output of apart should be simplifiable and still equal the input
        let e = Parser::parse("1/(x^2 - 1)").unwrap();
        let out = apart(&e, "x").unwrap();
        // Verify the output has partial fraction structure (contains division)
        let s = out.to_string();
        assert!(s.contains('/') || s.contains("frac"), "output should contain fractions: {}", s);
    }

    #[test]
    fn apart_cubic_denominator() {
        // 1/(x^3 - 1) = 1/((x-1)(x²+x+1))
        let src = "1/(x^3 - 1)";
        assert!(agrees(src, "x", &[0.0, 2.0, 3.0, -1.0, 0.5, -2.0]));
    }

    #[test]
    fn apart_quartic_with_repeated_root() {
        // 1/(x^2*(x-1)^2)
        let src = "1/(x^2 * (x - 1)^2)";
        assert!(agrees(src, "x", &[0.5, 2.0, 3.0, -1.0, -0.5, 4.0]));
    }
}
