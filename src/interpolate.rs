//! Polynomial interpolation from scratch.
//!
//! Provides Lagrange interpolation and Newton's divided-difference
//! interpolation for fitting a polynomial through a set of points.

use crate::error::{MathError, Result};

/// Lagrange interpolation: given points `(x_i, y_i)`, returns the
/// interpolated value at `x`.
pub fn lagrange_interp(points: &[(f64, f64)], x: f64) -> Result<f64> {
    if points.is_empty() {
        return Err(MathError::InvalidArgument("lagrange_interp: empty points".into()));
    }
    if points.len() == 1 {
        return Ok(points[0].1);
    }
    let n = points.len();
    let mut result = 0.0;
    for i in 0..n {
        let (xi, yi) = points[i];
        let mut term = yi;
        for j in 0..n {
            if i == j {
                continue;
            }
            let (xj, _) = points[j];
            let denom = xi - xj;
            if denom.abs() < 1e-14 {
                return Err(MathError::InvalidArgument(
                    "lagrange_interp: duplicate x values".into(),
                ));
            }
            term *= (x - xj) / denom;
        }
        result += term;
    }
    Ok(result)
}

/// Newton's divided-difference interpolation.
///
/// Builds the coefficient table once, then allows efficient evaluation
/// at multiple points via Newton's form.
pub struct NewtonInterpolator {
    coeffs: Vec<f64>,
    xs: Vec<f64>,
}

impl NewtonInterpolator {
    /// Build from a set of `(x, y)` points.
    pub fn new(points: &[(f64, f64)]) -> Result<Self> {
        if points.is_empty() {
            return Err(MathError::InvalidArgument("NewtonInterpolator: empty points".into()));
        }
        let n = points.len();
        let xs: Vec<f64> = points.iter().map(|(x, _)| *x).collect();
        let mut table = vec![vec![0.0; n]; n];
        for i in 0..n {
            table[i][0] = points[i].1;
        }
        for j in 1..n {
            for i in 0..(n - j) {
                let denom = xs[i] - xs[i + j];
                if denom.abs() < 1e-14 {
                    return Err(MathError::InvalidArgument(
                        "NewtonInterpolator: duplicate x values".into(),
                    ));
                }
                table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / denom;
            }
        }
        let coeffs: Vec<f64> = (0..n).map(|j| table[0][j]).collect();
        Ok(Self { coeffs, xs })
    }

    /// Evaluate the interpolating polynomial at `x` using Horner-like form.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.coeffs.len();
        let mut result = self.coeffs[n - 1];
        for i in (0..(n - 1)).rev() {
            result = result * (x - self.xs[i]) + self.coeffs[i];
        }
        result
    }

    /// Return the divided-difference coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.coeffs
    }
}

/// Convenience: Newton interpolation at a single point.
pub fn newton_interp(points: &[(f64, f64)], x: f64) -> Result<f64> {
    let interp = NewtonInterpolator::new(points)?;
    Ok(interp.eval(x))
}

/// Linear interpolation between two points (a special case of Lagrange).
pub fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn lagrange_through_points() {
        let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0)];
        // All points should be recovered exactly
        for &(x, y) in &pts {
            assert!(close(lagrange_interp(&pts, x).unwrap(), y, 1e-10));
        }
        // f(x) = x^2 + 1, so f(0.5) = 1.25
        assert!(close(lagrange_interp(&pts, 0.5).unwrap(), 1.25, 1e-10));
    }

    #[test]
    fn newton_through_points() {
        let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (3.0, 10.0)];
        let interp = NewtonInterpolator::new(&pts).unwrap();
        for &(x, y) in &pts {
            assert!(close(interp.eval(x), y, 1e-10));
        }
        // f(x) = x^2 + 1, so f(1.5) = 3.25
        assert!(close(interp.eval(1.5), 3.25, 1e-10));
    }

    #[test]
    fn lagrange_and_newton_agree() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 8.0), (3.0, 27.0)];
        for &x in &[0.5, 1.5, 2.5, 0.1, 2.9] {
            let l = lagrange_interp(&pts, x).unwrap();
            let n = newton_interp(&pts, x).unwrap();
            assert!(close(l, n, 1e-10), "at x={}: lagrange={} newton={}", x, l, n);
        }
    }

    #[test]
    fn single_point() {
        let pts = vec![(5.0, 42.0)];
        assert!(close(lagrange_interp(&pts, 100.0).unwrap(), 42.0, 1e-10));
    }

    #[test]
    fn lerp_basic() {
        assert!(close(lerp(0.0, 0.0, 10.0, 100.0, 5.0), 50.0, 1e-10));
        assert!(close(lerp(0.0, 0.0, 10.0, 100.0, 2.5), 25.0, 1e-10));
    }

    #[test]
    fn duplicate_x_errors() {
        let pts = vec![(1.0, 2.0), (1.0, 3.0)];
        assert!(lagrange_interp(&pts, 0.5).is_err());
        assert!(NewtonInterpolator::new(&pts).is_err());
    }
}
