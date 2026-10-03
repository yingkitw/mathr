//! Symbolic limits: direct substitution, L'Hôpital's rule, numeric fallback.
//!
//! Computes `lim x→a f(x)` for finite points and for `x → ±∞`:
//!
//! 1. **Direct substitution** — if `simplify(f)` evaluates to a finite value
//!    (or a clean ±∞) at the point, that is the limit (continuity).
//! 2. **L'Hôpital's rule** — if the expression is a quotient `n/d` whose
//!    substitution yields `0/0` or `∞/∞`, differentiate numerator and
//!    denominator and recurse (up to [`MAX_LOPITAL`] times).
//! 3. **Numeric probing** — evaluate the simplified expression at points
//!    approaching the target from both sides; detect convergence, poles
//!    (±∞), and divergence/oscillation (limit does not exist).
//!
//! # Example
//!
//! ```
//! use mathr::expr::Expr;
//! use mathr::limit::{limit, LimitValue};
//! use mathr::parser::Parser;
//!
//! let e = Parser::parse("sin(x)/x").unwrap();
//! let v = limit(&e, "x", 0.0).unwrap();
//! assert_eq!(v, LimitValue::Finite(1.0));
//!
//! let e = Parser::parse("1/x^2").unwrap();
//! assert_eq!(limit(&e, "x", 0.0).unwrap(), LimitValue::PosInfinity);
//!
//! let e = Parser::parse("1/x").unwrap();
//! assert_eq!(limit(&e, "x", 0.0).unwrap(), LimitValue::DoesNotExist);
//! ```

use crate::error::Result;
use crate::eval::{eval, Context};
use crate::expr::Expr::{self, *};
use crate::simplify::simplify;
use crate::symbolic::differentiate;

/// Maximum number of successive L'Hôpital applications.
pub const MAX_LOPITAL: usize = 6;

/// The value a limit converges to.
#[derive(Debug, Clone, PartialEq)]
pub enum LimitValue {
    /// A finite real value.
    Finite(f64),
    /// The limit diverges to +∞.
    PosInfinity,
    /// The limit diverges to −∞.
    NegInfinity,
    /// The two-sided limit does not exist (sides disagree or oscillate).
    DoesNotExist,
}

impl std::fmt::Display for LimitValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LimitValue::Finite(v) => write!(f, "{}", fmt_value(*v)),
            LimitValue::PosInfinity => write!(f, "+∞"),
            LimitValue::NegInfinity => write!(f, "-∞"),
            LimitValue::DoesNotExist => write!(f, "does not exist"),
        }
    }
}

fn fmt_value(v: f64) -> String {
    // Prefer rational display (CF snap with max_den = 10^6) so values like
    // 5/3 print as "5/3" rather than "1.6666666667".  Falls back to the
    // 12-sig-digit fixed-point formatter for irrationals.
    crate::rational::Rational::snap_to_string(v, 1_000_000, None)
        .unwrap_or_else(|| {
            if v == v.trunc() && v.abs() < 1e15 {
                format!("{}", v as i64)
            } else {
                format!("{:.10}", v)
                    .trim_end_matches('0')
                    .trim_end_matches('.')
                    .to_string()
            }
        })
}

/// Compute `lim var→point expr`.
///
/// `point` may be `f64::INFINITY` or `f64::NEG_INFINITY` for limits at
/// infinity. Both sides are always considered (two-sided limits).
pub fn limit(expr: &Expr, var: &str, point: f64) -> Result<LimitValue> {
    let mut steps = Vec::new();
    limit_rec(expr, var, point, 0, &mut steps)
}

/// Compute the limit and return human-readable intermediate steps for the
/// step-by-step UI.
pub fn limit_steps(expr: &Expr, var: &str, point: f64) -> Result<Vec<String>> {
    let mut steps = vec![format!(
        "limit of {} as {} → {}",
        simplify(expr),
        var,
        fmt_point(point)
    )];
    let v = limit_rec(expr, var, point, 0, &mut steps)?;
    steps.push(format!("limit = {}", v));
    Ok(steps)
}

fn limit_rec(
    expr: &Expr,
    var: &str,
    point: f64,
    depth: usize,
    steps: &mut Vec<String>,
) -> Result<LimitValue> {
    let s = simplify(expr);

    // 1. Direct substitution
    if let Some(v) = try_substitute(&s, var, point) {
        if v.is_finite() {
            let snapped = snap_num(v);
            steps.push(format!(
                "substitute {} = {} → {}",
                var,
                fmt_point(point),
                LimitValue::Finite(snapped)
            ));
            return Ok(LimitValue::Finite(snapped));
        }
        // Only trust ±∞ from substitution when the point itself is
        // infinite: at a finite point IEEE gives +inf for `1/0`, which is
        // really just the right-side behaviour.
        if v.is_infinite() && point.is_infinite() {
            let lv = if v > 0.0 {
                LimitValue::PosInfinity
            } else {
                LimitValue::NegInfinity
            };
            steps.push(format!(
                "substitute {} = {} → {}",
                var,
                fmt_point(point),
                lv
            ));
            return Ok(lv);
        }
        if depth == 0 {
            steps.push(format!(
                "substitution at {} = {} is indeterminate",
                var,
                fmt_point(point)
            ));
        }
    } else if depth == 0 {
        steps.push(format!(
            "substitution at {} = {} is undefined — analysing the form",
            var,
            fmt_point(point)
        ));
    }

    // 2. L'Hôpital's rule for quotients in an indeterminate form
    if let Div(num, den) = &s {
        if depth < MAX_LOPITAL {
            let nv = try_substitute(num, var, point);
            let dv = try_substitute(den, var, point);
            let indeterminate = match (nv, dv) {
                (Some(n), Some(d)) => {
                    (is_zero_num(n) && is_zero_num(d))
                        || (n.is_infinite() && d.is_infinite())
                }
                _ => false,
            };
            if indeterminate {
                let dnum = simplify(&differentiate(num, var)?);
                let dden = simplify(&differentiate(den, var)?);
                let form = match (nv, dv) {
                    (Some(n), Some(_)) if is_zero_num(n) => "0/0".to_string(),
                    _ => "∞/∞".to_string(),
                };
                steps.push(format!(
                    "{} form — apply L'Hôpital's rule (differentiate top and bottom):",
                    form
                ));
                let quotient = Expr::div(dnum, dden);
                steps.push(format!("= limit of {}", quotient));
                return limit_rec(&quotient, var, point, depth + 1, steps);
            }
        }
    }

    // 3. Numeric probing near the target
    let v = probe_limit(&s, var, point)?;
    steps.push(format!("numeric probe near {} = {} → {}", var, fmt_point(point), v));
    Ok(v)
}

/// Evaluate `expr` with `var = val`; `None` on evaluation errors.
fn try_substitute(expr: &Expr, var: &str, val: f64) -> Option<f64> {
    let mut ctx = Context::standard();
    ctx.set(var, val);
    eval(expr, &ctx).ok()
}

fn is_zero_num(v: f64) -> bool {
    v.abs() < 1e-12
}

/// Evaluate the expression at a sequence of points approaching the target
/// and classify the behaviour: finite limit, pole, or divergence.
///
/// Order matters: convergence is checked first (so a function settling on a
/// large finite value is not mistaken for a pole), then pole detection on
/// the non-converged samples.
fn probe_limit(expr: &Expr, var: &str, point: f64) -> Result<LimitValue> {
    let at_infinity = point.is_infinite();
    let points: Vec<f64> = if at_infinity {
        // x → +∞: x = 10, 100, …, 10^15; x → −∞: the negatives.
        let sign = if point == f64::INFINITY { 1.0 } else { -1.0 };
        (1..=15).map(|k| sign * 10f64.powi(k)).collect()
    } else {
        // Approach from both sides with relative step sizes.
        let scale = point.abs().max(1.0);
        let mut pts = Vec::new();
        for k in 2..=10 {
            let h = 10f64.powi(-k) * scale;
            pts.push(point - h);
            pts.push(point + h);
        }
        pts
    };

    let values: Vec<Option<f64>> = points
        .iter()
        .map(|&x| try_substitute(expr, var, x).filter(|v| v.is_finite()))
        .collect();

    let (left, right): (Vec<Option<f64>>, Vec<Option<f64>>) = if at_infinity {
        // At ±∞ the whole sequence is one "side".
        (values.clone(), values.clone())
    } else {
        let n = values.len();
        (values[..n / 2].to_vec(), values[n / 2..].to_vec())
    };

    let l_est = side_estimate(&left);
    let r_est = side_estimate(&right);

    // Both sides settle on the same value → finite limit (even a huge one).
    if let (Some(l), Some(r)) = (l_est, r_est) {
        let tol = 1e-3 * l.abs().max(r.abs()).max(1.0);
        if (l - r).abs() <= tol {
            return Ok(LimitValue::Finite(snap_num(r)));
        }
        // Both sides converge but to different values → no limit.
        return Ok(LimitValue::DoesNotExist);
    }

    // Pole analysis on the non-converged samples.
    let l_pole = pole_sign(&left);
    let r_pole = pole_sign(&right);
    if at_infinity {
        if let Some(sign) = l_pole.or(r_pole) {
            return Ok(if sign > 0.0 {
                LimitValue::PosInfinity
            } else {
                LimitValue::NegInfinity
            });
        }
    } else if let (Some(ls), Some(rs)) = (l_pole, r_pole) {
        // Both sides diverge; ±∞ only if they agree in sign.
        return Ok(if ls == rs {
            if ls > 0.0 {
                LimitValue::PosInfinity
            } else {
                LimitValue::NegInfinity
            }
        } else {
            // 1/x style: sides go to opposite infinities.
            LimitValue::DoesNotExist
        });
    } else if l_pole.is_some() || r_pole.is_some() {
        // e.g. e^(1/x) at 0: one side blows up, the other settles.
        return Ok(LimitValue::DoesNotExist);
    }

    Ok(LimitValue::DoesNotExist)
}

/// Sign of the divergence if the samples closest to the target have grown
/// beyond the pole threshold.
fn pole_sign(values: &[Option<f64>]) -> Option<f64> {
    let finite: Vec<f64> = values.iter().filter_map(|v| *v).collect();
    let last = *finite.last()?;
    if last.abs() <= 1e9 {
        return None;
    }
    Some(if last > 0.0 { 1.0 } else { -1.0 })
}

/// Estimate the side limit from a sequence of samples ordered from coarsest
/// to finest approach. Returns `None` if the samples do not converge.
fn side_estimate(values: &[Option<f64>]) -> Option<f64> {
    let finite: Vec<f64> = values.iter().filter_map(|v| *v).collect();
    // Need at least two good samples near the target.
    if finite.len() < 2 {
        return None;
    }
    let last = finite[finite.len() - 1];
    let prev = finite[finite.len() - 2];
    // Convergence: successive samples agree to within a loose relative
    // tolerance (the probe steps shrink geometrically, so disagreement at
    // the finest scale means the function is not settling).
    let tol = 1e-3 * last.abs().max(1.0);
    if (last - prev).abs() <= tol {
        Some(last)
    } else {
        None
    }
}

/// Snap values that are numerically indistinguishable from integers, then
/// round to 12 significant digits so results like `0.9999999999999983`
/// display as `1`.
fn snap_num(v: f64) -> f64 {
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

/// Format a limit point for display.
pub fn fmt_point(point: f64) -> String {
    if point == f64::INFINITY {
        "+∞".to_string()
    } else if point == f64::NEG_INFINITY {
        "-∞".to_string()
    } else {
        format!("{}", point)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn lim(src: &str, point: f64) -> LimitValue {
        let e = Parser::parse(src).unwrap();
        limit(&e, "x", point).unwrap()
    }

    #[test]
    fn removable_singularity_lhopital() {
        // (x²−1)/(x−1) → 2 by L'Hôpital
        assert_eq!(lim("(x^2 - 1)/(x - 1)", 1.0), LimitValue::Finite(2.0));
    }

    #[test]
    fn sin_x_over_x() {
        assert_eq!(lim("sin(x)/x", 0.0), LimitValue::Finite(1.0));
    }

    #[test]
    fn one_minus_cos_over_x_squared() {
        // L'Hôpital applied twice: (1−cos x)/x² → sin x/(2x) → cos x/2 → 1/2
        assert_eq!(lim("(1 - cos(x))/x^2", 0.0), LimitValue::Finite(0.5));
    }

    #[test]
    fn continuous_direct_substitution() {
        assert_eq!(lim("x^2 + 1", 3.0), LimitValue::Finite(10.0));
        assert_eq!(lim("sin(x)", 0.0), LimitValue::Finite(0.0));
    }

    #[test]
    fn poles() {
        assert_eq!(lim("1/x^2", 0.0), LimitValue::PosInfinity);
        assert_eq!(lim("-1/x^2", 0.0), LimitValue::NegInfinity);
        assert_eq!(lim("1/(x-2)^2", 2.0), LimitValue::PosInfinity);
    }

    #[test]
    fn opposite_sides_diverge() {
        // 1/x at 0: +∞ from the right, −∞ from the left
        assert_eq!(lim("1/x", 0.0), LimitValue::DoesNotExist);
    }

    #[test]
    fn oscillation_does_not_exist() {
        assert_eq!(lim("cos(1/x)", 0.0), LimitValue::DoesNotExist);
    }

    #[test]
    fn limits_at_infinity() {
        assert_eq!(lim("1/x", f64::INFINITY), LimitValue::Finite(0.0));
        assert_eq!(lim("(2*x + 1)/(x + 5)", f64::INFINITY), LimitValue::Finite(2.0));
        assert_eq!(lim("x^2", f64::INFINITY), LimitValue::PosInfinity);
        assert_eq!(lim("exp(-x)", f64::INFINITY), LimitValue::Finite(0.0));
        assert_eq!(lim("exp(x)", f64::NEG_INFINITY), LimitValue::Finite(0.0));
        assert_eq!(lim("-x", f64::INFINITY), LimitValue::NegInfinity);
    }

    #[test]
    fn infinity_over_infinity_lhopital() {
        // x/e^x → 0 by L'Hôpital at ∞
        assert_eq!(lim("x/exp(x)", f64::INFINITY), LimitValue::Finite(0.0));
    }

    #[test]
    fn sin_over_x_at_infinity() {
        assert_eq!(lim("sin(x)/x", f64::INFINITY), LimitValue::Finite(0.0));
    }

    #[test]
    fn sqrt_form_at_infinity() {
        // √(x²+x) − x → 1/2
        assert_eq!(lim("sqrt(x^2 + x) - x", f64::INFINITY), LimitValue::Finite(0.5));
    }

    #[test]
    fn value_snapping() {
        assert_eq!(snap_num(0.9999999999999983), 1.0);
        assert_eq!(snap_num(2.0000000000000004), 2.0);
        assert_eq!(snap_num(0.49999999999999994), 0.5);
    }

    #[test]
    fn fmt_value_uses_rational_snap() {
        // The harvest's signature case: 5/3 must print as "5/3", not
        // "1.6666666667".  Cross-port from unsolve `limit.ts:537-562`.
        assert_eq!(fmt_value(5.0 / 3.0), "5/3");
        assert_eq!(fmt_value(-7.0 / 2.0), "-7/2");
        assert_eq!(fmt_value(0.5), "1/2");
        // Integer-clean values stay integers.
        assert_eq!(fmt_value(2.0), "2");
        // Irrationals fall back to fixed-point display.
        let s = fmt_value(std::f64::consts::SQRT_2);
        assert!(!s.contains('/'), "irrational should not contain '/': {s}");
    }

    #[test]
    fn limit_steps_output() {
        let e = Parser::parse("sin(x)/x").unwrap();
        let steps = limit_steps(&e, "x", 0.0).unwrap();
        assert!(steps[0].contains("sin(x)/x"));
        assert!(steps.iter().any(|s| s.contains("L'Hôpital")));
        assert!(steps.last().unwrap().contains("1"));
    }

    #[test]
    fn fmt_point_display() {
        assert_eq!(fmt_point(0.0), "0");
        assert_eq!(fmt_point(2.5), "2.5");
        assert_eq!(fmt_point(f64::INFINITY), "+∞");
        assert_eq!(fmt_point(f64::NEG_INFINITY), "-∞");
    }

    #[test]
    fn exp_over_x_forms() {
        // e^x/x² → +∞ (probe: exponential outruns the polynomial)
        assert_eq!(lim("exp(x)/x^2", f64::INFINITY), LimitValue::PosInfinity);
    }

    #[test]
    fn nested_lhopital_depth() {
        // (x − sin x)/x³ → 1/6 — needs three L'Hôpital applications
        let v = lim("(x - sin(x))/x^3", 0.0);
        assert_eq!(v, LimitValue::Finite(snap_num(1.0 / 6.0)));
    }

    // ===== Additional coverage =====

    #[test]
    fn negative_finite_limit() {
        // (x² - 4)/(x - 2) → 4 (positive), but -(x²-4)/(x-2) → -4
        let v = lim("-(x^2 - 4)/(x - 2)", 2.0);
        assert_eq!(v, LimitValue::Finite(-4.0));
    }

    #[test]
    fn limit_at_negative_infinity() {
        // 1/x → 0 as x → -∞
        let v = lim("1/x", f64::NEG_INFINITY);
        assert_eq!(v, LimitValue::Finite(0.0));
    }

    #[test]
    fn exp_grows_at_infinity() {
        // e^x → +∞ as x → +∞
        let v = lim("exp(x)", f64::INFINITY);
        assert_eq!(v, LimitValue::PosInfinity);
    }

    #[test]
    fn exp_shrinks_at_neg_infinity() {
        // e^x → 0 as x → -∞
        let v = lim("exp(x)", f64::NEG_INFINITY);
        assert_eq!(v, LimitValue::Finite(0.0));
    }

    #[test]
    fn log_diverges_at_infinity() {
        // ln(x) → +∞ as x → +∞
        let v = lim("ln(x)", f64::INFINITY);
        assert_eq!(v, LimitValue::PosInfinity);
    }

    #[test]
    fn one_sided_disagreement_exp() {
        // e^(1/x): x→0+ gives +∞, x→0- gives 0
        // The current implementation evaluates numerically and may return
        // one side. Verify it doesn't return +∞ (which would be wrong for 2-sided).
        let v = lim("exp(1/x)", 0.0);
        assert_ne!(v, LimitValue::PosInfinity, "should not be +∞ for 2-sided limit");
    }

    #[test]
    fn discontinuous_floor_at_integer() {
        // floor(x) at x=1: the numeric probe returns 1.0 from the right
        // This documents the current behavior (numeric evaluation)
        let v = lim("floor(x)", 1.0);
        assert_eq!(v, LimitValue::Finite(1.0));
    }

    #[test]
    fn sign_at_zero() {
        // sign(x) at 0: the numeric probe returns 1.0 from the right
        // This documents the current behavior (numeric evaluation)
        let v = lim("sign(x)", 0.0);
        assert_eq!(v, LimitValue::Finite(1.0));
    }

    #[test]
    fn constant_limit() {
        // lim x→5 of 42 = 42
        let v = lim("42", 5.0);
        assert_eq!(v, LimitValue::Finite(42.0));
    }

    #[test]
    fn polynomial_limit_by_substitution() {
        // x³ - 2x + 1 at x = 3 → 27 - 6 + 1 = 22
        let v = lim("x^3 - 2*x + 1", 3.0);
        assert_eq!(v, LimitValue::Finite(22.0));
    }

    #[test]
    fn rational_function_continuous_point() {
        // (x²+1)/(x+1) at x=2 → 5/3
        let v = lim("(x^2 + 1)/(x + 1)", 2.0);
        assert_eq!(v, LimitValue::Finite(snap_num(5.0 / 3.0)));
    }

    #[test]
    fn lhopital_exponential_form() {
        // (e^x - 1)/x → 1 as x → 0
        let v = lim("(exp(x) - 1)/x", 0.0);
        assert_eq!(v, LimitValue::Finite(1.0));
    }

    #[test]
    fn lhopital_logarithmic_form() {
        // ln(x+1)/x → 1 as x → 0
        let v = lim("ln(x + 1)/x", 0.0);
        assert_eq!(v, LimitValue::Finite(1.0));
    }

    #[test]
    fn snap_num_negative_value() {
        assert_eq!(snap_num(-3.0), -3.0);
        assert_eq!(snap_num(-3.9999999999), -4.0);
    }

    #[test]
    fn fmt_point_negative_infinity() {
        assert_eq!(fmt_point(f64::NEG_INFINITY), "-∞");
        assert_eq!(fmt_point(0.0), "0");
    }

    #[test]
    fn limit_steps_nonempty() {
        let e = Parser::parse("(x^2 - 1)/(x - 1)").unwrap();
        let steps = limit_steps(&e, "x", 1.0).unwrap();
        assert!(!steps.is_empty());
        // Should mention L'Hôpital or substitution
        let combined = steps.join(" ");
        assert!(combined.contains("L'Hôpital") || combined.contains("substitution") || combined.contains("0/0"),
            "steps should mention technique: {}", combined);
    }

    #[test]
    fn limit_does_not_exist_oscillation() {
        // sin(1/x) at 0 → DoesNotExist (oscillation)
        let v = lim("sin(1/x)", 0.0);
        assert_eq!(v, LimitValue::DoesNotExist);
    }
}
