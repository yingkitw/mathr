//! Exact symbolic solution of linear systems.
//!
//! Solves a system of linear equations in 2–4 variables **exactly**, using
//! rational Gaussian elimination:
//!
//! ```text
//! 2x + y = 5
//! x + 3y = 10      →  x = 1, y = 3
//! ```
//!
//! The equations are parsed from text, so each may be given in either form:
//!
//! - `2x + y = 5` (an explicit `=`), or
//! - `2x + y - 5` (the right-hand side is zero).
//!
//! Coefficients are read as exact rationals, so the solution is exact rather
//! than a float approximation. A system that is singular (no unique
//! solution) is reported as such.
//!
//! This is the "multisolve" of Qalculate. Cross-port note: the `unsolve`
//! TypeScript port does not have this feature — it solves one variable at a
//! time.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::parser::Parser;
use crate::rational::Rational;
use std::collections::BTreeMap;

/// The largest system this solver handles.
pub const MAX_VARS: usize = 4;

/// A monomial: variable → exponent.
type Monomial = BTreeMap<String, u32>;

/// Read an expression as an exact rational polynomial.
///
/// Returns the coefficients per variable and the constant term.
fn to_rational_poly(expr: &Expr, vars: &[String]) -> Option<(Vec<Rational>, Rational)> {
    let mut terms: Vec<(Monomial, Rational)> = Vec::new();
    collect_rational_terms(expr, &mut terms)?;
    let mut coeffs = vec![Rational::from_int(0); vars.len()];
    let mut constant = Rational::from_int(0);
    for (mono, coeff) in terms {
        match mono.len() {
            0 => constant = constant + coeff,
            1 => {
                let (name, exp) = mono.iter().next().unwrap();
                if *exp != 1 {
                    // A squared variable is not linear.
                    return None;
                }
                let idx = vars.iter().position(|v| v == name)?;
                coeffs[idx] = coeffs[idx] + coeff;
            }
            _ => {
                // A product of variables is not linear.
                return None;
            }
        }
    }
    Some((coeffs, constant))
}

/// Collect the rational terms of an expression.
fn collect_rational_terms(
    expr: &Expr,
    out: &mut Vec<(Monomial, Rational)>,
) -> Option<()> {
    match expr {
        Expr::Num(x) => {
            let r = rational_of(*x)?;
            out.push((Monomial::new(), r));
            Some(())
        }
        Expr::Var(name) => {
            let mut m = Monomial::new();
            m.insert(name.clone(), 1);
            out.push((m, Rational::from_int(1)));
            Some(())
        }
        Expr::Neg(a) => {
            let mut inner = Vec::new();
            collect_rational_terms(a, &mut inner)?;
            for (_, r) in inner.iter_mut() {
                *r = -*r;
            }
            out.extend(inner);
            Some(())
        }
        Expr::Add(a, b) | Expr::Sub(a, b) => {
            let mut la = Vec::new();
            let mut lb = Vec::new();
            collect_rational_terms(a, &mut la)?;
            collect_rational_terms(b, &mut lb)?;
            if matches!(expr, Expr::Sub(_, _)) {
                for (_, r) in lb.iter_mut() {
                    *r = -*r;
                }
            }
            out.extend(la);
            out.extend(lb);
            Some(())
        }
        Expr::Mul(a, b) => {
            let mut la = Vec::new();
            let mut lb = Vec::new();
            collect_rational_terms(a, &mut la)?;
            collect_rational_terms(b, &mut lb)?;
            for (ma, ra) in la {
                for (mb, rb) in &lb {
                    let mut m = ma.clone();
                    let mut ok = true;
                    for (name, exp) in mb {
                        let e = m.entry(name.clone()).or_insert(0);
                        *e += exp;
                        if *e > 1 {
                            ok = false;
                        }
                    }
                    if !ok {
                        return None;
                    }
                    out.push((m, ra * *rb));
                }
            }
            Some(())
        }
        Expr::Div(a, b) => {
            // Only a constant denominator keeps it a rational polynomial.
            let mut lb = Vec::new();
            collect_rational_terms(b, &mut lb)?;
            if lb.len() != 1 || !lb[0].0.is_empty() {
                return None;
            }
            let mut la = Vec::new();
            collect_rational_terms(a, &mut la)?;
            for (_, r) in la.iter_mut() {
                *r = *r / lb[0].1;
            }
            out.extend(la);
            Some(())
        }
        Expr::Pow(base, exp) => {
            // Only a non-negative integer exponent of a linear base.
            let (Expr::Num(k),) = (&**exp,) else {
                return None;
            };
            if *k < 0.0 || k.fract() != 0.0 || *k > 1.0 {
                return None;
            }
            let mut lb = Vec::new();
            collect_rational_terms(base, &mut lb)?;
            if lb.len() != 1 || lb[0].0.len() > 1 {
                return None;
            }
            out.extend(lb);
            Some(())
        }
        Expr::Func(_, _) => None,
    }
}

/// Read a float as an exact rational, when it is one.
fn rational_of(v: f64) -> Option<Rational> {
    if !v.is_finite() {
        return None;
    }
    if v == v.trunc() && v.abs() < 9e15 {
        return Some(Rational::from_int(v as i64));
    }
    // A decimal fraction like 0.25 is exactly 1/4.
    let r = Rational::from_float_snap(v, crate::rational::SNAP_MAX_DEN)?;
    if (r.to_f64() - v).abs() > 1e-12 * v.abs().max(1.0) {
        return None;
    }
    Some(r)
}

/// Solve a system of linear equations given as text.
///
/// Each equation is separated by `;` (or given one per call via
/// [`solve_system_lines`]). The variables are inferred from the equations.
pub fn solve_system(text: &str) -> Result<Vec<(String, Rational)>> {
    let equations: Vec<&str> = text.split(';').map(|s| s.trim()).filter(|s| !s.is_empty()).collect();
    if equations.len() < 2 {
        return Err(MathError::InvalidArgument(
            "solve-system needs at least two equations separated by ';'".into(),
        ));
    }
    // Collect the variables across all equations.
    let mut vars: Vec<String> = Vec::new();
    let mut parsed: Vec<(Expr, Rational)> = Vec::new();
    for eq in &equations {
        let (lhs, rhs) = match eq.split_once('=') {
            Some((l, r)) => (l.trim(), r.trim()),
            None => (*eq, "0"),
        };
        let l = Parser::parse(lhs)?;
        let r = Parser::parse(rhs)?;
        // Move the right-hand side to the left: lhs − rhs = 0.
        let e = crate::simplify::simplify(&Expr::sub(l, r));
        for v in e.variables() {
            if !vars.contains(&v) {
                vars.push(v);
            }
        }
        parsed.push((e, Rational::from_int(0)));
    }
    vars.sort();
    if vars.is_empty() {
        return Err(MathError::InvalidArgument(
            "solve-system: the equations have no variables".into(),
        ));
    }
    if vars.len() > MAX_VARS {
        return Err(MathError::InvalidArgument(format!(
            "solve-system: {} variables exceeds the maximum of {}",
            vars.len(),
            MAX_VARS
        )));
    }
    if parsed.len() != vars.len() {
        return Err(MathError::InvalidArgument(format!(
            "solve-system: {} equations for {} variables — a unique solution needs one \
             equation per variable",
            parsed.len(),
            vars.len()
        )));
    }
    // Build the augmented matrix [A | b] over rationals.
    let mut aug: Vec<Vec<Rational>> = Vec::with_capacity(parsed.len());
    for (e, _) in &parsed {
        let (coeffs, constant) = to_rational_poly(e, &vars).ok_or_else(|| {
            MathError::InvalidArgument(format!(
                "solve-system: the equation '{}' is not linear in {}",
                e,
                vars.join(", ")
            ))
        })?;
        // lhs − rhs = 0 → lhs = rhs, so the constant term moves to the right.
        aug.push(
            coeffs
                .into_iter()
                .chain(std::iter::once(-constant))
                .collect(),
        );
    }
    let solution = gaussian_elimination(&mut aug, vars.len())?;
    Ok(vars.into_iter().zip(solution).collect())
}

/// Solve a system given one equation per line.
pub fn solve_system_lines(lines: &[String]) -> Result<Vec<(String, Rational)>> {
    solve_system(&lines.join("; "))
}

/// Exact Gaussian elimination over rationals.
///
/// `aug` is the augmented matrix [A | b], modified in place. Returns the
/// solution vector.
fn gaussian_elimination(aug: &mut [Vec<Rational>], n: usize) -> Result<Vec<Rational>> {
    let rows = aug.len();
    for col in 0..n {
        // Partial pivoting: find a row with a nonzero pivot.
        let mut pivot_row = None;
        for (r, row) in aug.iter().enumerate().skip(col) {
            if row[col].num() != 0 {
                pivot_row = Some(r);
                break;
            }
        }
        let Some(pivot_row) = pivot_row else {
            return Err(MathError::InvalidArgument(
                "the system is singular — it has no unique solution".into(),
            ));
        };
        if pivot_row != col {
            aug.swap(pivot_row, col);
        }
        // Eliminate below.
        let pivot = aug[col][col];
        // The pivot row is cloned so the current row can be updated in place.
        let pivot_row_data = aug[col].clone();
        // The 2-D indexing is clearer here than an iterator form.
        #[allow(clippy::needless_range_loop)]
        for r in (col + 1)..rows {
            let row = &mut aug[r];
            if row[col].num() == 0 {
                continue;
            }
            let factor = row[col] / pivot;
            for (c, cell) in row.iter_mut().enumerate().skip(col) {
                *cell = *cell - factor * pivot_row_data[c];
            }
        }
    }
    // Back-substitute.
    let mut x = vec![Rational::from_int(0); n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum = sum - aug[i][j] * x[j];
        }
        let pivot = aug[i][i];
        if pivot.num() == 0 {
            return Err(MathError::InvalidArgument(
                "the system is singular — it has no unique solution".into(),
            ));
        }
        x[i] = sum / pivot;
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn solve(text: &str) -> Vec<(String, String)> {
        solve_system(text)
            .unwrap()
            .into_iter()
            .map(|(v, r)| (v, r.to_string()))
            .collect()
    }

    #[test]
    fn solves_two_by_two() {
        // 2x + y = 5, x + 3y = 10  →  x = 1, y = 3
        let s = solve("2x + y = 5; x + 3y = 10");
        assert!(s.contains(&("x".to_string(), "1".to_string())), "got {:?}", s);
        assert!(s.contains(&("y".to_string(), "3".to_string())), "got {:?}", s);
    }

    #[test]
    fn solves_with_rational_answers() {
        // 2x + y = 1, x − y = 0  →  x = y, so 3x = 1  →  x = y = 1/3
        let s = solve("2x + y = 1; x - y = 0");
        assert!(s.contains(&("x".to_string(), "1/3".to_string())), "got {:?}", s);
        assert!(s.contains(&("y".to_string(), "1/3".to_string())), "got {:?}", s);
    }

    #[test]
    fn accepts_the_implicit_zero_form() {
        // Without an explicit '=', the right-hand side is zero.
        let s = solve("2x + y - 5; x + 3y - 10");
        assert!(s.contains(&("x".to_string(), "1".to_string())), "got {:?}", s);
    }

    #[test]
    fn solves_three_by_three() {
        // x + y + z = 6, 2x + 3y + z = 10, x + y + 2z = 9  →  x=2, y=1, z=3
        let s = solve("x + y + z = 6; 2x + 3y + z = 10; x + y + 2z = 9");
        assert!(s.contains(&("x".to_string(), "2".to_string())), "got {:?}", s);
        assert!(s.contains(&("y".to_string(), "1".to_string())), "got {:?}", s);
        assert!(s.contains(&("z".to_string(), "3".to_string())), "got {:?}", s);
    }

    #[test]
    fn solves_four_by_four() {
        let s = solve("x+y+z+w = 10; x-y+z-w = 0; x+y-z-w = 0; x-y-z+w = 0");
        // The solution is x=y=z=w=5/2.
        for v in ["x", "y", "z", "w"] {
            assert!(
                s.contains(&(v.to_string(), "5/2".to_string())),
                "got {:?}",
                s
            );
        }
    }

    #[test]
    fn rejects_a_singular_system() {
        // The second equation is a multiple of the first.
        let err = solve_system("x + y = 1; 2x + 2y = 2").unwrap_err();
        assert!(err.to_string().contains("singular"), "got {}", err);
    }

    #[test]
    fn rejects_a_wrong_number_of_equations() {
        let err = solve_system("x + y = 1").unwrap_err();
        assert!(err.to_string().contains("at least two"), "got {}", err);
        let err = solve_system("x + y = 1; x - y = 0; x = 2").unwrap_err();
        assert!(err.to_string().contains("one equation per variable"), "got {}", err);
    }

    #[test]
    fn rejects_nonlinear_equations() {
        // x² is not linear.
        assert!(solve_system("x^2 + y = 1; x + y = 0").is_err());
        // A product of variables is not linear.
        assert!(solve_system("x*y = 1; x + y = 0").is_err());
    }

    #[test]
    fn rejects_too_many_variables() {
        // Five variables — one more than the solver handles. (Note: `e` is
        // Euler's number in mathr, not a variable, so use other names.)
        let err = solve_system(
            "p+q+r+s+t = 1; p+q+r+s = 1; p+q+r = 1; p+q = 1; p = 1",
        )
        .unwrap_err();
        assert!(err.to_string().contains("maximum"), "got {}", err);
    }

    #[test]
    fn rational_coefficients_are_exact() {
        // (1/2)x + y = 1, x + (1/3)y = 1  →  x = 4/5, y = 3/5
        let s = solve("(1/2)x + y = 1; x + (1/3)y = 1");
        assert!(s.contains(&("x".to_string(), "4/5".to_string())), "got {:?}", s);
        assert!(s.contains(&("y".to_string(), "3/5".to_string())), "got {:?}", s);
    }
}
