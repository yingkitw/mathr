//! Levenberg–Marquardt nonlinear least-squares curve fitting.
//!
//! Fits a model `y = f(x, params)` to `(x, y)` observations by minimizing the
//! sum of squared residuals. Uses the LM damped Gauss–Newton algorithm:
//! each step solves `(JᵀJ + λ·diag(JᵀJ))·δ = −Jᵀr` with Marquardt scaling,
//! where the numeric Jacobian `J` is computed by central differences. The
//! damping parameter `λ` is decreased on successful steps and increased on
//! rejected ones, interpolating between Gauss–Newton (fast) and gradient
//! descent (robust).
//!
//! The linear algebra runs on [`crate::matrix::Matrix`]'s public API (solve,
//! transpose, inverse) — this module adds no solver code of its own.
//!
//! # Example
//!
//! ```
//! use mathr::curvefit::{curve_fit, LmOptions};
//!
//! // y = 3 e^{-0.7 x}, noisy-free synthetic fit from a bad starting point
//! let data: Vec<(f64, f64)> = (0..10)
//!     .map(|i| { let x = i as f64; (x, 3.0 * (-0.7 * x).exp()) })
//!     .collect();
//! let fit = curve_fit(
//!     |x, p| Ok(p[0] * (-p[1] * x).exp()),
//!     &data,
//!     &[1.0, 1.0],
//!     &LmOptions::default(),
//! )
//! .unwrap();
//! assert!(fit.converged);
//! assert!((fit.params[0] - 3.0).abs() < 1e-6);
//! assert!((fit.params[1] - 0.7).abs() < 1e-6);
//! ```

use crate::error::{MathError, Result};
use crate::matrix::Matrix;

/// Tuning parameters for [`curve_fit`]. The defaults work well for typical
/// small problems; see field docs for the knobs.
#[derive(Debug, Clone)]
pub struct LmOptions {
    /// Maximum outer iterations (default 200).
    pub max_iter: usize,
    /// Convergence tolerance on relative SSE decrease and step size
    /// (default 1e-12).
    pub tol: f64,
    /// Initial damping parameter (default 1e-3).
    pub lambda0: f64,
    /// Multiplier applied to `lambda` when a step is rejected (default 10).
    pub lambda_up: f64,
    /// Multiplier applied to `lambda` when a step is accepted (default 0.1).
    pub lambda_down: f64,
}

impl Default for LmOptions {
    fn default() -> Self {
        LmOptions {
            max_iter: 200,
            tol: 1e-12,
            lambda0: 1e-3,
            lambda_up: 10.0,
            lambda_down: 0.1,
        }
    }
}

/// Outcome of a successful [`curve_fit`] run (convergence is reported via
/// [`LmFit::converged`], not an error).
#[derive(Debug, Clone, PartialEq)]
pub struct LmFit {
    /// Fitted parameter values.
    pub params: Vec<f64>,
    /// Sum of squared residuals at the fitted parameters.
    pub sse: f64,
    /// Outer iterations performed.
    pub iterations: usize,
    /// True if the convergence criteria were met within `max_iter`.
    pub converged: bool,
    /// Approximate 1-sigma parameter standard errors from
    /// `σ²·(JᵀJ)⁻¹` with `σ² = sse/(m−n)`. Empty when `m ≤ n` (no residual
    /// degrees of freedom).
    pub std_errors: Vec<f64>,
}

/// Fit `model(x, params)` to `data` by Levenberg–Marquardt least squares.
///
/// `init` provides the starting parameter vector (length `n`); the model is
/// called with `m = data.len()` x-values. Returns the best parameters found
/// with `converged = false` if the tolerance was not met within
/// `opts.max_iter`.
pub fn curve_fit<F>(mut model: F, data: &[(f64, f64)], init: &[f64], opts: &LmOptions) -> Result<LmFit>
where
    F: FnMut(f64, &[f64]) -> Result<f64>,
{
    if data.is_empty() {
        return Err(MathError::InvalidArgument("curve_fit: no data points".into()));
    }
    if init.is_empty() {
        return Err(MathError::InvalidArgument("curve_fit: no parameters".into()));
    }
    if data.len() < init.len() {
        return Err(MathError::InvalidArgument(format!(
            "curve_fit: {} data points cannot fit {} parameters (need m >= n)",
            data.len(),
            init.len()
        )));
    }
    let n = init.len();
    let mut params = init.to_vec();
    if params.iter().any(|p| !p.is_finite()) {
        return Err(MathError::InvalidArgument("curve_fit: non-finite initial guess".into()));
    }

    let residuals = |model: &mut F, p: &[f64]| -> Result<Vec<f64>> {
        let mut r = Vec::with_capacity(data.len());
        for &(x, y) in data {
            let f = model(x, p)?;
            if !f.is_finite() {
                return Err(MathError::Domain(format!(
                    "model produced non-finite value at x = {x}"
                )));
            }
            r.push(f - y);
        }
        Ok(r)
    };
    let sse_of = |r: &[f64]| r.iter().map(|v| v * v).sum::<f64>();

    let mut r = residuals(&mut model, &params)?;
    let mut sse = sse_of(&r);
    let mut lambda = opts.lambda0;
    let mut iterations = 0usize;

    while iterations < opts.max_iter {
        iterations += 1;
        let j = numeric_jacobian(&mut model, data, &params)?;
        let a = (&j.transpose() * &j)?;
        let g = j.transpose().mul_vec(&r)?;
        let step_tol = opts.tol * params.iter().fold(1.0f64, |m, p| m.max(p.abs()));

        // Damping sub-loop: raise lambda until a step is accepted or lambda
        // maxes out (then the problem is too ill-conditioned to progress).
        let mut accepted = false;
        let mut new_sse = sse;
        let mut step = vec![0.0; n];
        while lambda <= 1e12 {
            // (JᵀJ + λ·diag(JᵀJ)) δ = −g  with diagonal floored at 1e-12
            let mut damp = Vec::with_capacity(n * n);
            for i in 0..n {
                for jj in 0..n {
                    let a_ij = a[(i, jj)];
                    damp.push(if i == jj { a_ij + lambda * a_ij.max(1e-12) } else { a_ij });
                }
            }
            let rhs: Vec<f64> = g.iter().map(|v| -v).collect();
            let delta = Matrix::from_row_major(n, n, damp)?.solve(&rhs)?;
            let trial: Vec<f64> = params
                .iter()
                .zip(delta.iter())
                .map(|(p, d)| p + d)
                .collect();
            let tr = residuals(&mut model, &trial)?;
            let tsse = sse_of(&tr);
            if tsse < sse {
                accepted = true;
                new_sse = tsse;
                step = delta;
                params = trial;
                r = tr;
                lambda = (lambda * opts.lambda_down).max(1e-12);
                break;
            }
            lambda *= opts.lambda_up;
        }
        if !accepted {
            return Ok(finish(params, sse, iterations, false, &a, data.len(), n));
        }

        // Convergence: negligible SSE decrease or negligible step.
        let sse_drop = sse - new_sse;
        let step_inf = step.iter().fold(0.0f64, |m, d| m.max(d.abs()));
        sse = new_sse;
        let a_final = (&j.transpose() * &j)?;
        if sse_drop <= opts.tol * sse.max(1.0) || step_inf <= step_tol {
            return Ok(finish(params, sse, iterations, true, &a_final, data.len(), n));
        }
    }

    let j = numeric_jacobian(&mut model, data, &params)?;
    let a = (&j.transpose() * &j)?;
    Ok(finish(params, sse, iterations, false, &a, data.len(), n))
}

/// Compute standard errors and assemble the result.
fn finish(
    params: Vec<f64>,
    sse: f64,
    iterations: usize,
    converged: bool,
    jtj: &Matrix,
    m: usize,
    n: usize,
) -> LmFit {
    let mut std_errors = Vec::new();
    if m > n {
        let s2 = sse / (m - n) as f64;
        if let Ok(cov) = jtj.inverse() {
            for idx in 0..n {
                let v = (s2 * cov[(idx, idx)]).max(0.0).sqrt();
                std_errors.push(v);
            }
        }
    }
    LmFit { params, sse, iterations, converged, std_errors }
}

/// Central-difference numeric Jacobian of the residuals w.r.t. parameters.
fn numeric_jacobian<F>(model: &mut F, data: &[(f64, f64)], params: &[f64]) -> Result<Matrix>
where
    F: FnMut(f64, &[f64]) -> Result<f64>,
{
    let m = data.len();
    let n = params.len();
    let mut col_major = Vec::with_capacity(m * n);
    for (j, p0) in params.iter().enumerate() {
        let h = 1e-6 * p0.abs().max(1.0);
        let mut pp = params.to_vec();
        let mut pm = params.to_vec();
        pp[j] = p0 + h;
        pm[j] = p0 - h;
        for &(x, _y) in data {
            let f_plus = model(x, &pp)?;
            let f_minus = model(x, &pm)?;
            if !(f_plus.is_finite() && f_minus.is_finite()) {
                return Err(MathError::Domain(format!(
                    "model produced non-finite value at x = {x}"
                )));
            }
            col_major.push((f_plus - f_minus) / (2.0 * h));
        }
    }
    // col_major is column-by-column; rebuild row-major J (m×n).
    let mut row_major = vec![0.0; m * n];
    for jj in 0..n {
        for ii in 0..m {
            row_major[ii * n + jj] = col_major[jj * m + ii];
        }
    }
    Matrix::from_row_major(m, n, row_major)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fit_linear(data: &[(f64, f64)]) -> LmFit {
        curve_fit(|x, p| Ok(p[0] * x + p[1]), data, &[0.5, 0.5], &LmOptions::default()).unwrap()
    }

    #[test]
    fn exact_linear_recovery() {
        // y = 2x + 1 at x = 0, 1, 2 -> exact fit
        let data = [(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)];
        let fit = fit_linear(&data);
        assert!(fit.converged);
        assert!((fit.params[0] - 2.0).abs() < 1e-10);
        assert!((fit.params[1] - 1.0).abs() < 1e-10);
        assert!(fit.sse < 1e-22, "sse={}", fit.sse);
        // m = n + 1 -> 1 residual dof -> std errors near zero for exact fit
        assert_eq!(fit.std_errors.len(), 2);
        assert!(fit.std_errors.iter().all(|e| *e < 1e-8));
    }

    #[test]
    fn matches_linear_regression() {
        // overdetermined noisy-ish data: LM should equal the closed form
        let data = [(0.0, 1.0), (1.0, 2.1), (2.0, 2.9), (3.0, 4.2)];
        let (slope, intercept) = crate::stats::linear_regression(
            &data.iter().map(|d| d.0).collect::<Vec<f64>>(),
            &data.iter().map(|d| d.1).collect::<Vec<f64>>(),
        )
        .unwrap();
        let fit = fit_linear(&data);
        assert!(fit.converged);
        assert!((fit.params[0] - slope).abs() < 1e-8, "a={} vs {slope}", fit.params[0]);
        assert!((fit.params[1] - intercept).abs() < 1e-8);
        // sse equals the regression residuals
        let sse_ref: f64 = data
            .iter()
            .map(|&(x, y)| (slope * x + intercept - y).powi(2))
            .sum();
        assert!((fit.sse - sse_ref).abs() < 1e-14);
    }

    #[test]
    fn exponential_decay_recovery() {
        // y = 3 e^{-0.7 x}, x = 0..9; init far away [1, 1]
        let data: Vec<(f64, f64)> = (0..10)
            .map(|i| {
                let x = i as f64;
                (x, 3.0 * (-0.7 * x).exp())
            })
            .collect();
        let fit = curve_fit(
            |x, p| Ok(p[0] * (-p[1] * x).exp()),
            &data,
            &[1.0, 1.0],
            &LmOptions::default(),
        )
        .unwrap();
        assert!(fit.converged);
        assert!((fit.params[0] - 3.0).abs() < 1e-5, "a={}", fit.params[0]);
        assert!((fit.params[1] - 0.7).abs() < 1e-5, "b={}", fit.params[1]);
        assert!(fit.sse < 1e-18);
    }

    #[test]
    fn quadratic_recovery_with_zero_params() {
        // y = x^2 exactly; true params [1, 0, 0] exercise a zero parameter
        let data: Vec<(f64, f64)> = (-4..=4).map(|i| (i as f64, (i * i) as f64)).collect();
        let fit = curve_fit(
            |x, p| Ok(p[0] * x * x + p[1] * x + p[2]),
            &data,
            &[0.5, 0.5, 0.5],
            &LmOptions::default(),
        )
        .unwrap();
        assert!(fit.converged);
        assert!((fit.params[0] - 1.0).abs() < 1e-10);
        assert!(fit.params[1].abs() < 1e-8, "b={}", fit.params[1]);
        assert!(fit.params[2].abs() < 1e-8, "c={}", fit.params[2]);
        assert_eq!(fit.std_errors.len(), 3);
        assert!(fit.std_errors.iter().all(|e| e.is_finite() && *e >= 0.0));
    }

    #[test]
    fn michaelis_menten_recovery() {
        // y = p0·x / (p1 + x), true [4, 1]; canonical unique-optimum NLS model
        let data: Vec<(f64, f64)> = (1..=20)
            .map(|i| {
                let x = 0.25 * i as f64;
                (x, 4.0 * x / (1.0 + x))
            })
            .collect();
        let fit = curve_fit(
            |x, p| Ok(p[0] * x / (p[1] + x)),
            &data,
            &[1.0, 1.0],
            &LmOptions::default(),
        )
        .unwrap();
        assert!(fit.converged);
        assert!((fit.params[0] - 4.0).abs() < 1e-5, "v={}", fit.params[0]);
        assert!((fit.params[1] - 1.0).abs() < 1e-5, "k={}", fit.params[1]);
    }

    #[test]
    fn sine_fit_needs_init_near_solution() {
        // Frequency fitting sin(w·x) has a local minimum near w = 0
        // (small-angle regime); LM from init 1.0 lands there. From an init
        // near the true frequency it converges. Documents init sensitivity.
        let data: Vec<(f64, f64)> = (0..20)
            .map(|i| {
                let x = 0.25 * i as f64;
                (x, (2.0 * x).sin())
            })
            .collect();
        let bad = curve_fit(
            |x, p| Ok((p[0] * x).sin()),
            &data,
            &[1.0],
            &LmOptions::default(),
        )
        .unwrap();
        // from w=1 the fit does NOT reach the true frequency 2
        assert!((bad.params[0] - 2.0).abs() > 1e-3, "w={}", bad.params[0]);
        let good = curve_fit(
            |x, p| Ok((p[0] * x).sin()),
            &data,
            &[1.5],
            &LmOptions::default(),
        )
        .unwrap();
        assert!(good.converged);
        assert!((good.params[0] - 2.0).abs() < 1e-6, "w={}", good.params[0]);
    }

    #[test]
    fn max_iter_zero_reports_not_converged() {
        let data = [(0.0, 1.0), (1.0, 2.0)];
        let opts = LmOptions { max_iter: 0, ..LmOptions::default() };
        let fit = curve_fit(|x, p| Ok(p[0] * x + p[1]), &data, &[0.0, 0.0], &opts).unwrap();
        assert!(!fit.converged);
        assert_eq!(fit.iterations, 0);
        assert_eq!(fit.params, vec![0.0, 0.0]);
    }

    #[test]
    fn std_errors_empty_without_residual_dof() {
        // m = n = 2: no residual degrees of freedom -> no std errors
        let data = [(0.0, 1.0), (1.0, 3.0)];
        let fit = fit_linear(&data);
        assert!(fit.converged);
        assert!(fit.std_errors.is_empty());
    }

    #[test]
    fn validation_errors() {
        assert!(curve_fit(|_x, p| Ok(p[0]), &[], &[1.0], &LmOptions::default()).is_err());
        assert!(curve_fit(|_x, p| Ok(p[0]), &[(1.0, 1.0)], &[], &LmOptions::default()).is_err());
        // more parameters than data points
        assert!(
            curve_fit(|_x, p| Ok(p[0] + p[1]), &[(1.0, 1.0)], &[1.0, 1.0], &LmOptions::default())
                .is_err()
        );
        // model error propagates
        let r = curve_fit(
            |_x, _p| Err(MathError::Domain("boom".into())),
            &[(1.0, 1.0)],
            &[1.0],
            &LmOptions::default(),
        );
        assert!(r.is_err());
    }
}
