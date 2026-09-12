//! B-spline basis, curves, and interpolating splines.
//!
//! Cox–de Boor recursion for B-spline basis functions, a general
//! [`BSpline`] curve evaluated with de Boor's algorithm, and
//! [`cubic_bspline_interp`] which fits a clamped B-spline through data
//! points using de Boor's knot averaging (guaranteeing a banded,
//! nonsingular collocation system). Degree is reduced automatically when
//! fewer than 4 points are given (linear for 2 points, quadratic for 3).

use crate::error::{MathError, Result};
use crate::matrix::Matrix;

/// B-spline basis function `N_{i,p}(t)` via the Cox–de Boor recursion.
///
/// Zero-degree basis functions are right-continuous. At or beyond the final
/// knot the function is evaluated as its left-hand limit (for a clamped
/// vector the limit mass sits on the last coefficient), which keeps
/// evaluation at the right end of the domain well-defined for every degree.
pub fn basis_function(i: usize, p: usize, knots: &[f64], t: f64) -> f64 {
    let m = knots.len() - 1;
    let width = knots[m] - knots[0];
    let t = if width > 0.0 && t >= knots[m] {
        knots[m] - 1e-12 * width
    } else {
        t
    };
    if p == 0 {
        return if t >= knots[i] && t < knots[i + 1] {
            1.0
        } else {
            0.0
        };
    }
    let den1 = knots[i + p] - knots[i];
    let left = if den1 > 0.0 {
        (t - knots[i]) / den1 * basis_function(i, p - 1, knots, t)
    } else {
        0.0
    };
    let den2 = knots[i + p + 1] - knots[i + 1];
    let right = if den2 > 0.0 {
        (knots[i + p + 1] - t) / den2 * basis_function(i + 1, p - 1, knots, t)
    } else {
        0.0
    };
    left + right
}

/// A B-spline curve: control coefficients over a knot vector, evaluated
/// with de Boor's algorithm. Outside the domain evaluation clamps to the
/// endpoints (like `CubicSpline`).
///
/// Evaluation happens at `u = (x - t0) * t_scale`; [`BSpline::new`] uses the
/// identity map (`t0 = 0, t_scale = 1`), while [`cubic_bspline_interp`]
/// maps data abscissae onto the normalized knot range.
#[derive(Debug, Clone)]
pub struct BSpline {
    degree: usize,
    knots: Vec<f64>,
    coeffs: Vec<f64>,
    t0: f64,
    t_scale: f64,
}

impl BSpline {
    /// Build a B-spline from a knot vector and control coefficients.
    ///
    /// Requires `knots.len() == coeffs.len() + degree + 1`, a nondecreasing
    /// knot vector, `degree >= 1`, and at least `degree + 1` coefficients.
    pub fn new(degree: usize, knots: Vec<f64>, coeffs: Vec<f64>) -> Result<Self> {
        if degree < 1 {
            return Err(MathError::InvalidArgument(format!(
                "degree must be at least 1, got {}",
                degree
            )));
        }
        if coeffs.len() < degree + 1 {
            return Err(MathError::InvalidArgument(format!(
                "need at least degree + 1 = {} coefficients, got {}",
                degree + 1,
                coeffs.len()
            )));
        }
        if knots.len() != coeffs.len() + degree + 1 {
            return Err(MathError::InvalidArgument(format!(
                "knot count {} must equal coefficient count {} + degree {} + 1",
                knots.len(),
                coeffs.len(),
                degree
            )));
        }
        for w in knots.windows(2) {
            if w[1] < w[0] {
                return Err(MathError::InvalidArgument(format!(
                    "knot vector must be nondecreasing, found {} > {}",
                    w[0], w[1]
                )));
            }
        }
        Ok(BSpline {
            degree,
            knots,
            coeffs,
            t0: 0.0,
            t_scale: 1.0,
        })
    }

    /// The degree `p` of the piecewise polynomials.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// The knot vector.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }

    /// The control coefficients.
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// The domain `[knots[p], knots[len-1-p]]` on which the curve is defined,
    /// in evaluation (x) coordinates.
    pub fn domain(&self) -> (f64, f64) {
        let m = self.knots.len() - 1;
        (
            self.t0 + self.knots[self.degree] / self.t_scale,
            self.t0 + self.knots[m - self.degree] / self.t_scale,
        )
    }

    /// Evaluate the curve at `x` (clamped into the domain).
    pub fn eval(&self, x: f64) -> f64 {
        let u = (x - self.t0) * self.t_scale;
        let p = self.degree;
        let m = self.knots.len() - 1;
        let u = u.clamp(self.knots[p], self.knots[m - p]);
        // Knot span: largest k in [p, m-p-1] with knots[k] <= u < knots[k+1];
        // u == right end lands on the final span.
        let mut k = p;
        while k < m - p - 1 && u >= self.knots[k + 1] {
            k += 1;
        }
        let mut d: Vec<f64> = (0..=p).map(|j| self.coeffs[j + k - p]).collect();
        for r in 1..=p {
            for j in (r..=p).rev() {
                let i = j + k - p;
                let den = self.knots[i + p - r + 1] - self.knots[i];
                let alpha = if den > 0.0 {
                    (u - self.knots[i]) / den
                } else {
                    0.0
                };
                d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
            }
        }
        d[p]
    }
}

/// Interpolate data points `(x_i, y_i)` with a clamped cubic B-spline.
///
/// Parameters are the normalized abscissae `t_i = (x_i - x_0)/(x_{n-1} - x_0)`;
/// interior knots follow de Boor's averaging formula (the Schoenberg–Whitney
/// condition holds, so the collocation system is banded and nonsingular) and
/// the control coefficients are found with a dense solve. With fewer than 4
/// points the degree drops (2 points → linear, 3 → quadratic). The result
/// passes through every data point and is evaluated on `[x_0, x_{n-1}]`
/// (clamped outside).
pub fn cubic_bspline_interp(xs: &[f64], ys: &[f64]) -> Result<BSpline> {
    let n = xs.len();
    if n != ys.len() {
        return Err(MathError::InvalidArgument(format!(
            "x and y lengths differ: {} vs {}",
            n,
            ys.len()
        )));
    }
    if n < 2 {
        return Err(MathError::InvalidArgument(
            "need at least 2 data points".into(),
        ));
    }
    for w in xs.windows(2) {
        if w[1] <= w[0] {
            return Err(MathError::InvalidArgument(format!(
                "x values must be strictly increasing, found {} then {}",
                w[0], w[1]
            )));
        }
    }
    let x0 = xs[0];
    let span = xs[n - 1] - x0;
    let params: Vec<f64> = xs.iter().map(|x| (x - x0) / span).collect();
    let p = n.min(4) - 1;

    // Clamped knot vector with de Boor averaging:
    // u_{j+p} = (t_j + ... + t_{j+p-1}) / p for j = 1..=n-p-1.
    let m = n + p;
    let mut knots = vec![0.0; m + 1];
    let pf = p as f64;
    for j in 1..=(n - p - 1) {
        let sum: f64 = (0..p).map(|i| params[j + i]).sum();
        knots[j + p] = sum / pf;
    }
    for k in knots.iter_mut().take(m + 1).skip(n) {
        *k = 1.0;
    }

    // Collocation: N(t_i) · c = y_i.
    let mut a = vec![0.0; n * n];
    for row in 0..n {
        for col in 0..n {
            a[row * n + col] = basis_function(col, p, &knots, params[row]);
        }
    }
    let mat = Matrix::from_row_major(n, n, a)?;
    let coeffs = mat.solve(ys)?;
    let mut spline = BSpline::new(p, knots, coeffs)?;
    spline.t0 = x0;
    spline.t_scale = 1.0 / span;
    Ok(spline)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn basis_functions_partition_of_unity() {
        // Clamped cubic knot vector with 5 coefficients.
        let knots = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        let p = 3usize;
        let n = knots.len() - p - 1;
        let ts = [0.0, 0.05, 0.25, 0.5, 0.5 + 1e-9, 0.75, 1.0];
        for &t in &ts {
            let total: f64 = (0..n)
                .map(|i| basis_function(i, p, &knots, t))
                .sum();
            assert!(
                close(total, 1.0, 1e-10),
                "partition of unity fails at t = {}: {}",
                t,
                total
            );
            for i in 0..n {
                let v = basis_function(i, p, &knots, t);
                assert!(v >= 0.0, "basis {} negative at t = {}: {}", i, t, v);
            }
        }
    }

    #[test]
    fn basis_functions_have_correct_support() {
        let knots = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
        // N_{0,3}(t) = (1 - 2t)³ on its support [0, 0.5).
        assert_abs_diff_eq!(basis_function(0, 3, &knots, 0.25), 0.125, epsilon = 1e-12);
        assert_abs_diff_eq!(basis_function(0, 3, &knots, 0.5), 0.0);
        // N_{4,3} peaks at the right end (clamped: left-limit equals 1 at t = 1).
        assert_abs_diff_eq!(basis_function(4, 3, &knots, 1.0), 1.0, epsilon = 1e-9);
        assert_abs_diff_eq!(basis_function(4, 3, &knots, 0.0), 0.0);
    }

    #[test]
    fn bspline_new_validates() {
        assert!(BSpline::new(0, vec![0.0, 1.0], vec![1.0, 1.0]).is_err());
        assert!(BSpline::new(1, vec![0.0, 1.0], vec![1.0]).is_err());
        // knots.len() must be coeffs.len() + degree + 1
        assert!(BSpline::new(1, vec![0.0, 0.5, 1.0], vec![1.0, 1.0]).is_err());
        assert!(BSpline::new(1, vec![1.0, 0.0, 1.0], vec![1.0, 1.0]).is_err());
    }

    #[test]
    fn linear_bspline_is_exact_between_two_points() {
        // n = 2 → degree 1, knots [0, 0, 1, 1]: exact straight line.
        let s = cubic_bspline_interp(&[0.0, 4.0], &[1.0, 3.0]).unwrap();
        assert_eq!(s.degree(), 1);
        for t in [0.0, 1.0, 2.5, 4.0] {
            assert_abs_diff_eq!(s.eval(t), 1.0 + 0.5 * t, epsilon = 1e-10);
        }
        assert_abs_diff_eq!(s.eval(5.0), 3.0, epsilon = 1e-10); // clamped
    }

    #[test]
    fn cubic_bspline_reproduces_cubics_exactly() {
        // Four points of a cubic → single Bezier segment → exact reproduction.
        let xs: Vec<f64> = (0..4).map(|i| i as f64).collect();
        let ys: Vec<f64> = xs.iter().map(|&x| x * x * x + 2.0 * x).collect();
        let s = cubic_bspline_interp(&xs, &ys).unwrap();
        assert_eq!(s.degree(), 3);
        for k in 1..30 {
            let x = k as f64 / 10.0;
            let expected = x * x * x + 2.0 * x;
            assert!(
                close(s.eval(x), expected, 1e-9),
                "eval({}) = {} vs {}",
                x,
                s.eval(x),
                expected
            );
        }
    }

    #[test]
    fn quadratic_bspline_for_three_points() {
        // Three points of a parabola → degree 2, exact reproduction.
        let xs = [0.0, 1.0, 2.0];
        let ys: Vec<f64> = xs.iter().map(|&x| x * x - 3.0 * x + 1.0).collect();
        let s = cubic_bspline_interp(&xs, &ys).unwrap();
        assert_eq!(s.degree(), 2);
        for k in 1..20 {
            let x = k as f64 / 10.0;
            let expected = x * x - 3.0 * x + 1.0;
            assert!(close(s.eval(x), expected, 1e-9));
        }
    }

    #[test]
    fn bspline_interpolant_is_knot_exact_on_sin_data() {
        let xs: Vec<f64> = (0..8).map(|i| i as f64 * 0.5).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
        let s = cubic_bspline_interp(&xs, &ys).unwrap();
        for i in 0..xs.len() {
            assert!(
                close(s.eval(xs[i]), ys[i], 1e-9),
                "not knot-exact at x = {}: {} vs {}",
                xs[i],
                s.eval(xs[i]),
                ys[i]
            );
        }
        // Accuracy between knots: cubic interpolation error ~ h⁴.
        for k in 1..32 {
            let x = k as f64 / 20.0;
            assert!(close(s.eval(x), x.sin(), 1e-3), "at x = {}", x);
        }
        assert_eq!(s.domain(), (0.0, 3.5));
    }

    #[test]
    fn bspline_eval_clamps_outside_domain() {
        let s = cubic_bspline_interp(&[0.0, 1.0, 2.0, 3.0], &[0.0, 1.0, 4.0, 9.0]).unwrap();
        assert_abs_diff_eq!(s.eval(-1.0), s.eval(0.0), epsilon = 1e-12);
        assert_abs_diff_eq!(s.eval(4.0), s.eval(3.0), epsilon = 1e-12);
    }

    #[test]
    fn bspline_interp_validates_input() {
        assert!(cubic_bspline_interp(&[0.0], &[1.0]).is_err());
        assert!(cubic_bspline_interp(&[0.0, 1.0], &[1.0, 2.0, 3.0]).is_err());
        assert!(cubic_bspline_interp(&[1.0, 0.0], &[1.0, 2.0]).is_err());
        assert!(cubic_bspline_interp(&[0.0, 0.0], &[1.0, 2.0]).is_err());
    }

    #[test]
    fn arbitrary_knot_vector_evaluates() {
        // Manual construction: degree-1 tent over knots [0, 1, 2, 3, 4]
        // with coefficients [0, 2, 0]: rises 0→2 on [1, 2], falls 2→0 on [2, 3].
        let s = BSpline::new(1, vec![0.0, 1.0, 2.0, 3.0, 4.0], vec![0.0, 2.0, 0.0]).unwrap();
        assert_abs_diff_eq!(s.eval(1.0), 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(s.eval(1.5), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(s.eval(2.0), 2.0, epsilon = 1e-12);
        assert_abs_diff_eq!(s.eval(2.5), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(s.eval(3.0), 0.0, epsilon = 1e-12);
        assert_eq!(s.domain(), (1.0, 3.0));
        assert_abs_diff_eq!(s.eval(0.0), 0.0, epsilon = 1e-12); // clamped
    }
}
