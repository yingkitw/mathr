//! Logistic regression via iteratively reweighted least squares (IRLS).
//!
//! Fits the binary-response model `P(y=1|x) = sigmoid(β0 + β1·x1 + … + βk·xk)`
//! by maximizing the log-likelihood with the standard IRLS (Newton) update
//! `β ← (XᵀWX)⁻¹XᵀWz`, where `W = diag(p(1−p))` and `z = Xβ + (y−p)/(p(1−p))`
//! is the working response. Each iteration solves one linear system through
//! [`crate::matrix::Matrix`]'s public API — this module adds no solver code.
//!
//! # Example
//!
//! ```
//! use mathr::logit::{logistic_regression, LogitOptions};
//!
//! // balanced two-point design: y rates 1/3 at x=1 and 2/3 at x=2
//! // -> closed form: b1 = 2 ln 2, b0 = -3 ln 2
//! let x = [&[1.0, 1.0, 1.0, 2.0, 2.0, 2.0][..]];
//! let y = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0];
//! let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
//! assert!(fit.converged);
//! assert!((fit.coefficients[1] - 2.0 * std::f64::consts::LN_2).abs() < 1e-8);
//! ```

use crate::error::{MathError, Result};
use crate::matrix::Matrix;

/// Tuning parameters for [`logistic_regression`].
#[derive(Debug, Clone)]
pub struct LogitOptions {
    /// Maximum IRLS iterations (default 100).
    pub max_iter: usize,
    /// Convergence tolerance on the infinity-norm coefficient update,
    /// relative to `1 + ‖β‖∞` (default 1e-10).
    pub tol: f64,
}

impl Default for LogitOptions {
    fn default() -> Self {
        LogitOptions { max_iter: 100, tol: 1e-10 }
    }
}

/// Fitted logistic regression model.
#[derive(Debug, Clone, PartialEq)]
pub struct LogitFit {
    /// `[intercept, b1, …, bk]`.
    pub coefficients: Vec<f64>,
    /// Wald standard errors: `sqrt(diag((XᵀWX)⁻¹))` at convergence.
    pub std_errors: Vec<f64>,
    /// Log-likelihood at the fitted coefficients.
    pub log_likelihood: f64,
    /// IRLS iterations performed.
    pub iterations: usize,
    /// True if the coefficient updates fell below `tol` within `max_iter`.
    /// `false` often indicates (quasi-)complete separation, where the MLE
    /// does not exist and coefficients diverge.
    pub converged: bool,
}

/// Overflow-safe logistic function.
fn sigmoid(eta: f64) -> f64 {
    if eta >= 0.0 {
        1.0 / (1.0 + (-eta).exp())
    } else {
        let e = eta.exp();
        e / (1.0 + e)
    }
}

/// Predicted probability of the positive class for predictor values
/// `features` (length `k`, no intercept) under fitted `coefficients`
/// `[intercept, b1, …, bk]`.
pub fn predict_proba(coefficients: &[f64], features: &[f64]) -> Result<f64> {
    if coefficients.len() != features.len() + 1 {
        return Err(MathError::InvalidArgument(format!(
            "predict_proba: {} coefficients need {} features",
            coefficients.len(),
            coefficients.len() - 1
        )));
    }
    let eta = coefficients[0]
        + coefficients[1..]
            .iter()
            .zip(features.iter())
            .map(|(b, x)| b * x)
            .sum::<f64>();
    Ok(sigmoid(eta))
}

/// Fit logistic regression by IRLS.
///
/// `x` holds the predictor columns (each of length `m = y.len()`); `y` holds
/// the binary responses (exactly 0.0 or 1.0). Returns coefficients
/// `[intercept, b1, …, bk]` with Wald standard errors.
#[allow(clippy::needless_range_loop)] // design-matrix access needs (row, col) indices
pub fn logistic_regression(x: &[&[f64]], y: &[f64], opts: &LogitOptions) -> Result<LogitFit> {
    if y.is_empty() {
        return Err(MathError::InvalidArgument("logistic_regression: no observations".into()));
    }
    if y.iter().any(|&v| v != 0.0 && v != 1.0) {
        return Err(MathError::InvalidArgument(
            "logistic_regression: responses must be exactly 0 or 1".into(),
        ));
    }
    for (j, col) in x.iter().enumerate() {
        if col.len() != y.len() {
            return Err(MathError::InvalidArgument(format!(
                "logistic_regression: predictor {} has {} values, expected {}",
                j + 1,
                col.len(),
                y.len()
            )));
        }
    }
    let m = y.len();
    let np = x.len() + 1; // design columns: intercept + predictors

    // Design row i: [1, x1_i, …, xk_i]
    let design = |i: usize, c: usize| if c == 0 { 1.0 } else { x[c - 1][i] };

    let mut beta = vec![0.0; np];
    let mut iterations = 0usize;
    let mut converged = false;

    // XᵀWX (symmetric) and XᵀWz built per iteration; + 1e-10 ridge on the
    // diagonal keeps the solve well-defined under (quasi-)separation.
    for _ in 0..opts.max_iter {
        iterations += 1;
        // weighted least squares pieces
        let mut xtwx = vec![0.0; np * np];
        let mut xtwz = vec![0.0; np];
        for i in 0..m {
            let mut eta = 0.0;
            for c in 0..np {
                eta += beta[c] * design(i, c);
            }
            let p = sigmoid(eta);
            let w = (p * (1.0 - p)).max(1e-10);
            let z = eta + (y[i] - p) / w;
            for a in 0..np {
                let xa = design(i, a);
                xtwz[a] += xa * w * z;
                for b in a..np {
                    let v = xa * design(i, b) * w;
                    xtwx[a * np + b] += v;
                    if a != b {
                        xtwx[b * np + a] += v;
                    }
                }
            }
        }
        for d in 0..np {
            xtwx[d * np + d] += 1e-10;
        }
        let next = Matrix::from_row_major(np, np, xtwx)?.solve(&xtwz)?;

        let scale = 1.0 + beta.iter().fold(0.0f64, |mx, v| mx.max(v.abs()));
        let step_inf = beta
            .iter()
            .zip(next.iter())
            .fold(0.0f64, |mx, (a, b)| mx.max((a - b).abs()));
        beta = next;
        if step_inf <= opts.tol * scale {
            converged = true;
            break;
        }
    }

    // Final weights for log-likelihood and standard errors.
    let mut log_likelihood = 0.0;
    let mut xtwx = vec![0.0; np * np];
    for i in 0..m {
        let mut eta = 0.0;
        for c in 0..np {
            eta += beta[c] * design(i, c);
        }
        let p = sigmoid(eta);
        // stable log-likelihood: ln p = -ln(1+e^-eta), ln(1-p) = -ln(1+e^eta)
        log_likelihood += if y[i] == 1.0 {
            -((-eta).exp().ln_1p())
        } else {
            -(eta.exp().ln_1p())
        };
        let w = (p * (1.0 - p)).max(1e-10);
        for a in 0..np {
            for b in a..np {
                let v = design(i, a) * design(i, b) * w;
                xtwx[a * np + b] += v;
                if a != b {
                    xtwx[b * np + a] += v;
                }
            }
        }
    }
    let mut std_errors = Vec::with_capacity(np);
    if let Ok(cov) = Matrix::from_row_major(np, np, xtwx)?.inverse() {
        for d in 0..np {
            std_errors.push((cov[(d, d)]).max(0.0).sqrt());
        }
    } else {
        std_errors = vec![f64::NAN; np];
    }

    Ok(LogitFit {
        coefficients: beta,
        std_errors,
        log_likelihood,
        iterations,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs().max(a.abs()))
    }

    #[test]
    fn balanced_two_point_closed_form() {
        // y rates 1/3 at x=1, 2/3 at x=2 -> empirical-logit closed form:
        // b1 = logit(2/3) - logit(1/3) = 2 ln 2, b0 = -3 ln 2, SE(b1) = sqrt(3)
        let x = [&[1.0, 1.0, 1.0, 2.0, 2.0, 2.0][..]];
        let y = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0];
        let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
        assert!(fit.converged);
        assert!(close(fit.coefficients[0], -3.0 * std::f64::consts::LN_2, 1e-8));
        assert!(close(fit.coefficients[1], 2.0 * std::f64::consts::LN_2, 1e-8));
        // grouped-information SE: 1/(m p(1-p)) summed over the two groups = 3
        assert!(close(fit.std_errors[1], 3.0f64.sqrt(), 1e-6), "se={}", fit.std_errors[1]);
        // fitted probabilities reproduce the empirical rates
        assert!(close(predict_proba(&fit.coefficients, &[1.0]).unwrap(), 1.0 / 3.0, 1e-8));
        assert!(close(predict_proba(&fit.coefficients, &[2.0]).unwrap(), 2.0 / 3.0, 1e-8));
        // closed-form log-likelihood: 4 ln(2/3) + 2 ln(1/3) = 4 ln2 - 6 ln3
        let expect_ll = 4.0 * std::f64::consts::LN_2 - 6.0 * (3.0f64).ln();
        assert!(close(fit.log_likelihood, expect_ll, 1e-8), "ll={}", fit.log_likelihood);
    }

    #[test]
    fn balanced_overlap_gives_null_model() {
        // x=1 and x=2 both have rate 1/2 -> MLE is beta = 0 (p = 0.5 for all)
        let x = [&[1.0, 1.0, 2.0, 2.0][..]];
        let y = [0.0, 1.0, 0.0, 1.0];
        let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
        assert!(fit.converged);
        for c in &fit.coefficients {
            assert!(c.abs() < 1e-8, "coef={c}");
        }
    }

    #[test]
    fn score_equation_zero_at_convergence() {
        // gradient of the log-likelihood, Xᵀ(y − p), vanishes at the MLE
        let x = [
            &[0.5, 1.2, 2.0, 2.7, 3.1, 4.0, 4.8, 5.5][..],
            &[1.0, 0.3, 1.7, 0.9, 2.2, 1.1, 0.5, 2.0][..],
        ];
        let y = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
        assert!(fit.converged);
        // Σ(y_i − p_i) and Σx_j(y_i − p_i) over the design rows
        let mut grad = vec![0.0; fit.coefficients.len()];
        for i in 0..y.len() {
            grad[0] += y[i];
            for (j, col) in x.iter().enumerate() {
                grad[j + 1] += y[i] * col[i];
            }
        }
        for i in 0..y.len() {
            let features: Vec<f64> = x.iter().map(|col| col[i]).collect();
            let p = predict_proba(&fit.coefficients, &features).unwrap();
            grad[0] -= p;
            for (j, col) in x.iter().enumerate() {
                grad[j + 1] -= p * col[i];
            }
        }
        for (j, g) in grad.iter().enumerate() {
            assert!(g.abs() < 1e-7, "grad[{j}]={g}");
        }
    }

    #[test]
    fn coefficient_sign_follows_direction() {
        let up = [&[1.0, 2.0, 3.0, 4.0][..]];
        let y_up = [0.0, 0.0, 1.0, 1.0];
        let fit = logistic_regression(&up, &y_up, &LogitOptions::default()).unwrap();
        assert!(fit.coefficients[1] > 0.0);
        let y_down = [1.0, 1.0, 0.0, 0.0];
        let fit2 = logistic_regression(&up, &y_down, &LogitOptions::default()).unwrap();
        assert!(fit2.coefficients[1] < 0.0);
    }

    #[test]
    fn complete_separation_diverges_but_predicts_correct_side() {
        // perfectly separable data: the MLE does not exist (|beta| -> inf),
        // IRLS stops once steps fall below tolerance with a large slope
        let x = [&[1.0, 2.0, 3.0, 4.0][..]];
        let y = [0.0, 0.0, 1.0, 1.0];
        let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
        assert!(
            fit.coefficients[1] > 10.0,
            "separation should push the slope out: {:?}",
            fit.coefficients
        );
        // yet the direction is still correct
        assert!(fit.coefficients[1] > 0.0);
        assert!(predict_proba(&fit.coefficients, &[1.0]).unwrap() < 0.5);
        assert!(predict_proba(&fit.coefficients, &[4.0]).unwrap() > 0.5);
    }

    #[test]
    fn multivariate_gradient_and_prediction() {
        // duplicated design points with mixed labels -> provably
        // non-separable; grouped rates 1/3 (x2=0) vs 2/3 (x2=1) and no x1
        // effect give the closed form beta = [-ln 2, 0, 2 ln 2]
        let x = [
            &[0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0][..],
            &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0][..],
        ];
        let y = [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0];
        let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();
        assert!(fit.converged);
        assert!(close(fit.coefficients[0], -std::f64::consts::LN_2, 1e-8));
        assert!(close(fit.coefficients[1], 0.0, 1e-6));
        assert!(close(fit.coefficients[2], 2.0 * std::f64::consts::LN_2, 1e-8));
        assert_eq!(fit.coefficients.len(), 3);
        assert_eq!(fit.std_errors.len(), 3);
        // training-set probabilities are self-consistent
        for (&a, &b) in x[0].iter().zip(x[1].iter()) {
            let p = predict_proba(&fit.coefficients, &[a, b]).unwrap();
            assert!((0.0..=1.0).contains(&p));
        }
        // non-degenerate fit: log-likelihood strictly negative and finite
        assert!(fit.log_likelihood.is_finite() && fit.log_likelihood < 0.0);
        // score equations hold at the MLE
        let mut grad = [0.0; 3];
        for i in 0..y.len() {
            let p = predict_proba(&fit.coefficients, &[x[0][i], x[1][i]]).unwrap();
            grad[0] += y[i] - p;
            grad[1] += (y[i] - p) * x[0][i];
            grad[2] += (y[i] - p) * x[1][i];
        }
        for (j, g) in grad.iter().enumerate() {
            assert!(g.abs() < 1e-6, "grad[{j}]={g}");
        }
    }

    #[test]
    fn validation_errors() {
        let x = [&[1.0, 2.0][..]];
        assert!(logistic_regression(&x, &[], &LogitOptions::default()).is_err());
        // non-binary response
        assert!(logistic_regression(&x, &[0.0, 0.5], &LogitOptions::default()).is_err());
        // length mismatch
        assert!(logistic_regression(&[&[1.0, 2.0, 3.0]], &[0.0, 1.0], &LogitOptions::default())
            .is_err());
        // coefficient/feature mismatch in predict_proba
        assert!(predict_proba(&[0.0, 1.0], &[1.0, 2.0]).is_err());
    }
}
