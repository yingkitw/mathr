//! Monotone cubic Hermite interpolation (PCHIP, Fritsch–Carlson).
//!
//! Piecewise cubic interpolation that passes exactly through every knot and
//! — unlike natural cubic splines — never overshoots the local data range.
//! On each interval `[x_k, x_{k+1}]` the interpolant is the cubic Hermite
//! `h00·y_k + h10·h·m_k + h01·y_{k+1} + h11·h·m_{k+1}` where the slopes `m`
//! come from the Fritsch–Carlson construction:
//!
//! - secants `δ_k = Δy_k / h_k`,
//! - interior knots with a sign change of neighbouring secants get `m = 0`
//!   (flat spots), otherwise the weighted harmonic mean of the two secants,
//! - endpoint three-point formulas, clipped to preserve monotonicity.
//!
//! Evaluation outside `[x_0, x_n]` clamps to the endpoint values (keeps the
//! monotonicity guarantees; documented deviation from scipy, which
//! extrapolates).
//!
//! # Example
//!
//! ```
//! use mathr::pchip::Pchip;
//!
//! // linear data is reproduced exactly
//! let p = Pchip::new(&[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0]).unwrap();
//! assert!((p.eval(0.5) - 2.0).abs() < 1e-14);
//! assert!((p.eval(1.7) - 4.4).abs() < 1e-14);
//! ```

use crate::error::{MathError, Result};

/// Monotone piecewise cubic Hermite interpolant over sorted knots.
#[derive(Debug, Clone, PartialEq)]
pub struct Pchip {
    xs: Vec<f64>,
    ys: Vec<f64>,
    /// Fritsch–Carlson derivative estimates at each knot.
    slopes: Vec<f64>,
}

impl Pchip {
    /// Build the interpolant from `xs` (strictly increasing) and `ys`.
    pub fn new(xs: &[f64], ys: &[f64]) -> Result<Pchip> {
        if xs.len() != ys.len() {
            return Err(MathError::InvalidArgument(format!(
                "pchip: {} x-values vs {} y-values",
                xs.len(),
                ys.len()
            )));
        }
        if xs.len() < 2 {
            return Err(MathError::InvalidArgument("pchip: need at least 2 points".into()));
        }
        if xs.windows(2).any(|w| w[1] <= w[0]) || xs.iter().any(|x| !x.is_finite()) {
            return Err(MathError::InvalidArgument(
                "pchip: x-values must be finite and strictly increasing".into(),
            ));
        }
        Ok(Pchip {
            xs: xs.to_vec(),
            ys: ys.to_vec(),
            slopes: fritsch_carlson_slopes(xs, ys),
        })
    }

    /// The Fritsch–Carlson slope estimates at the knots.
    pub fn slopes(&self) -> &[f64] {
        &self.slopes
    }

    /// Evaluate at `x`. Values outside `[x_first, x_last]` clamp to the
    /// endpoint y-values.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.xs.len();
        if x <= self.xs[0] {
            return self.ys[0];
        }
        if x >= self.xs[n - 1] {
            return self.ys[n - 1];
        }
        // rightmost interval with xs[k] <= x
        let k = self.xs.partition_point(|&v| v <= x) - 1;
        let h = self.xs[k + 1] - self.xs[k];
        let t = (x - self.xs[k]) / h;
        // cubic Hermite basis
        let h00 = (1.0 + 2.0 * t) * (1.0 - t) * (1.0 - t);
        let h10 = t * (1.0 - t) * (1.0 - t);
        let h01 = t * t * (3.0 - 2.0 * t);
        let h11 = t * t * (t - 1.0);
        h00 * self.ys[k]
            + h10 * h * self.slopes[k]
            + h01 * self.ys[k + 1]
            + h11 * h * self.slopes[k + 1]
    }
}

/// Fritsch–Carlson slope estimates (non-uniform grid variant).
fn fritsch_carlson_slopes(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    let n = xs.len();
    let h: Vec<f64> = (0..n - 1).map(|k| xs[k + 1] - xs[k]).collect();
    let delta: Vec<f64> = (0..n - 1).map(|k| (ys[k + 1] - ys[k]) / h[k]).collect();
    // two points: the cubic degenerates to the exact straight line
    if n == 2 {
        return vec![delta[0], delta[0]];
    }
    let mut m = vec![0.0; n];

    // interior knots
    for k in 1..n - 1 {
        if delta[k - 1] * delta[k] <= 0.0 {
            m[k] = 0.0; // local extremum or flat: slope 0
        } else {
            let w1 = 2.0 * h[k] + h[k - 1];
            let w2 = h[k] + 2.0 * h[k - 1];
            m[k] = (w1 + w2) / (w1 / delta[k - 1] + w2 / delta[k]);
        }
    }

    // left endpoint: three-point formula, clipped for monotonicity
    m[0] = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1]);
    if m[0] * delta[0] <= 0.0 {
        m[0] = 0.0;
    } else if delta[0] * delta[1] <= 0.0 && m[0].abs() > 3.0 * delta[0].abs() {
        m[0] = 3.0 * delta[0];
    }

    // right endpoint (mirror)
    let r = n - 1;
    m[r] = ((2.0 * h[r - 1] + h[r - 2]) * delta[r - 1] - h[r - 1] * delta[r - 2])
        / (h[r - 1] + h[r - 2]);
    if m[r] * delta[r - 1] <= 0.0 {
        m[r] = 0.0;
    } else if delta[r - 1] * delta[r - 2] <= 0.0 && m[r].abs() > 3.0 * delta[r - 1].abs() {
        m[r] = 3.0 * delta[r - 1];
    }

    m
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs().max(a.abs()))
    }

    #[test]
    fn linear_data_reproduced_exactly() {
        // two secants equal -> slopes equal the common slope -> exact line
        let p = Pchip::new(&[0.0, 1.0, 2.0], &[1.0, 3.0, 5.0]).unwrap();
        for &x in &[0.0, 0.25, 0.5, 1.0, 1.3, 1.7, 2.0] {
            assert!(close(p.eval(x), 1.0 + 2.0 * x, 1e-13), "x={x}");
        }
    }

    #[test]
    fn two_points_exact_line() {
        let p = Pchip::new(&[0.0, 4.0], &[2.0, -2.0]).unwrap();
        for &x in &[0.0, 1.0, 2.0, 3.0, 4.0] {
            assert!(close(p.eval(x), 2.0 - x, 1e-13));
        }
        // clamped outside
        assert_eq!(p.eval(-5.0), 2.0);
        assert_eq!(p.eval(9.0), -2.0);
    }

    #[test]
    fn interpolates_through_knots() {
        let xs = [0.0, 0.5, 1.5, 3.0, 7.0];
        let ys = [1.0, -2.0, 4.0, -1.0, 3.0];
        let p = Pchip::new(&xs, &ys).unwrap();
        for (&x, &y) in xs.iter().zip(ys.iter()) {
            assert!(close(p.eval(x), y, 1e-13), "x={x}");
        }
    }

    #[test]
    fn constant_data_is_constant() {
        let p = Pchip::new(&[0.0, 1.0, 2.0, 3.0], &[5.0; 4]).unwrap();
        assert!(p.slopes().iter().all(|&m| m == 0.0));
        for i in 0..60 {
            // Hermite basis weights sum to 1 up to one ulp
            assert!(close(p.eval(i as f64 / 20.0), 5.0, 1e-14));
        }
    }

    #[test]
    fn no_overshoot_on_monotone_data() {
        // THE PCHIP property: the interpolant stays within the local knot
        // range on every interval (natural cubic splines overshoot here)
        let xs = [0.0f64, 1.0, 2.0, 3.0, 4.0];
        let ys = [0.0, 0.0, 0.3, 4.0, 4.2];
        let p = Pchip::new(&xs, &ys).unwrap();
        for k in 0..xs.len() - 1 {
            let (lo, hi) = (ys[k].min(ys[k + 1]), ys[k].max(ys[k + 1]));
            for s in 1..50 {
                let x = xs[k] + (xs[k + 1] - xs[k]) * s as f64 / 50.0;
                let y = p.eval(x);
                assert!(y >= lo - 1e-12 && y <= hi + 1e-12, "overshoot at x={x}: {y}");
            }
        }
        // global monotonicity of the sampled interpolant
        let mut prev = p.eval(0.0);
        for i in 1..=400 {
            let x = 4.0 * i as f64 / 400.0;
            let y = p.eval(x);
            assert!(y >= prev - 1e-12, "not monotone at x={x}");
            prev = y;
        }
    }

    #[test]
    fn slopes_zero_at_local_extrema() {
        // y rises to a peak at knot 1 then falls: knot 1 slope is 0,
        // the falling side stays negative
        let p = Pchip::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.0, 0.5, -1.0]).unwrap();
        let s = p.slopes();
        assert_eq!(s[1], 0.0, "peak slope should be 0");
        assert!(s[0] > 0.0, "rising side slope positive");
        assert!(s[2] < 0.0 && s[3] < 0.0, "falling side slopes negative");
    }

    #[test]
    fn slope_construction_hand_check_uniform() {
        // uniform grid, all secants = 1: interior slope = weighted harmonic
        // mean of equal values = 1; endpoints also 1 (well inside 3x clip)
        let p = Pchip::new(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.0, 2.0, 3.0]).unwrap();
        for &m in p.slopes() {
            assert!(close(m, 1.0, 1e-14));
        }
    }

    #[test]
    fn non_uniform_grid() {
        // x = [0, 1, 4], y linear 0,2,8 (slope 2 everywhere) -> exact
        let p = Pchip::new(&[0.0, 1.0, 4.0], &[0.0, 2.0, 8.0]).unwrap();
        for &x in &[0.5, 2.0, 3.5] {
            assert!(close(p.eval(x), 2.0 * x, 1e-13));
        }
    }

    #[test]
    fn smooth_function_accuracy() {
        // sin over [0, pi]: monotone increasing, error ~ O(h^3) because the
        // Fritsch-Carlson slopes approximate f' (vs exact-slope Hermite's h^4)
        let xs: Vec<f64> = (0..9).map(|i| std::f64::consts::PI * i as f64 / 8.0).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
        let p = Pchip::new(&xs, &ys).unwrap();
        let mut worst = 0.0f64;
        for i in 1..80 {
            let x = std::f64::consts::PI * i as f64 / 80.0;
            worst = worst.max((p.eval(x) - x.sin()).abs());
        }
        assert!(worst < 1e-2, "9-knot error={worst}");
        // doubling the knots cuts the max error ~4x (order 3)
        let xs2: Vec<f64> = (0..17).map(|i| std::f64::consts::PI * i as f64 / 16.0).collect();
        let ys2: Vec<f64> = xs2.iter().map(|x| x.sin()).collect();
        let p2 = Pchip::new(&xs2, &ys2).unwrap();
        let mut worst2 = 0.0f64;
        for i in 1..160 {
            let x = std::f64::consts::PI * i as f64 / 160.0;
            worst2 = worst2.max((p2.eval(x) - x.sin()).abs());
        }
        let ratio = worst / worst2;
        assert!(ratio > 2.5 && ratio < 5.5, "refinement ratio={ratio}");
    }

    #[test]
    fn validation_errors() {
        assert!(Pchip::new(&[1.0], &[1.0]).is_err());
        assert!(Pchip::new(&[1.0, 2.0], &[1.0]).is_err());
        assert!(Pchip::new(&[2.0, 1.0], &[1.0, 2.0]).is_err()); // unsorted
        assert!(Pchip::new(&[1.0, 1.0], &[1.0, 2.0]).is_err()); // duplicate
        assert!(Pchip::new(&[f64::NAN, 1.0], &[1.0, 2.0]).is_err());
    }
}
