//! Laurent series expansion from scratch.
//!
//! Given an expression `f(x)` with a pole of order `k` at `x = a`, computes
//! the Laurent series:
//!
//! ```text
//! f(x) = Σ_{n=-k}^{N} c_n (x - a)^n
//! ```
//!
//! The principal part (n = -k, ..., -1) is obtained by Taylor-expanding
//! `g(x) = (x - a)^k · f(x)` (which is analytic at `a`) and dividing by
//! `(x - a)^k`. The analytic part (n = 0, ..., N) is the Taylor expansion
//! of `f(x) - principal_part`.

use crate::error::{MathError, Result};
use crate::eval::{eval, Context};
use crate::expr::Expr;
use crate::parser::Parser;

/// A Laurent series expansion: `Σ_{n=-k}^{N} c_n (x - a)^n`.
#[derive(Debug, Clone)]
pub struct LaurentSeries {
    /// Coefficients indexed from `n = -k` to `n = N`.
    /// `coeffs[0]` is `c_{-k}`, `coeffs[i]` is `c_{i - k}`.
    pub coeffs: Vec<f64>,
    /// Pole order (number of negative-power terms).
    pub pole_order: usize,
    /// Expansion point.
    pub center: f64,
    /// Variable name.
    pub var: String,
}

impl LaurentSeries {
    /// Get the coefficient for power `(x - a)^n`.
    pub fn coeff(&self, n: isize) -> f64 {
        let idx = (n + self.pole_order as isize) as usize;
        if idx < self.coeffs.len() {
            self.coeffs[idx]
        } else {
            0.0
        }
    }

    /// Number of terms (including zero coefficients).
    pub fn n_terms(&self) -> usize {
        self.coeffs.len()
    }

    /// Evaluate the series at point `x`.
    pub fn eval(&self, x: f64) -> f64 {
        let mut sum = 0.0;
        let k = self.pole_order as isize;
        for (i, &c) in self.coeffs.iter().enumerate() {
            let n = i as isize - k;
            let term = if n >= 0 {
                c * (x - self.center).powi(n as i32)
            } else {
                c / (x - self.center).powi((-n) as i32)
            };
            sum += term;
        }
        sum
    }

    /// Render as a human-readable string.
    pub fn to_string(&self) -> String {
        let k = self.pole_order as isize;
        let var = &self.var;
        let a = self.center;
        let mut parts: Vec<String> = Vec::new();

        for (i, &c) in self.coeffs.iter().enumerate() {
            if c.abs() < 1e-15 {
                continue;
            }
            let n = i as isize - k;
            let coeff_str = if c == 1.0 {
                "".to_string()
            } else if c == -1.0 {
                "-".to_string()
            } else {
                format!("{}*", c)
            };

            let power_str = if a == 0.0 {
                if n == 0 {
                    "1".to_string()
                } else if n == 1 {
                    var.clone()
                } else if n == -1 {
                    format!("1/{}", var)
                } else if n > 0 {
                    format!("{}^{}", var, n)
                } else {
                    format!("1/{}^{}", var, -n)
                }
            } else {
                let base = format!("({} - {})", var, a);
                if n == 0 {
                    "1".to_string()
                } else if n == 1 {
                    base
                } else if n == -1 {
                    format!("1/{}", base)
                } else if n > 0 {
                    format!("{}^{}", base, n)
                } else {
                    format!("1/{}^{}", base, -n)
                }
            };

            let term = format!("{}{}", coeff_str, power_str);
            if parts.is_empty() {
                parts.push(if c < 0.0 && c == -1.0 {
                    format!("-{}", power_str)
                } else {
                    term
                });
            } else if c < 0.0 {
                parts.push(format!(" - {}", term.trim_start_matches('-')));
            } else {
                parts.push(format!(" + {}", term));
            }
        }

        if parts.is_empty() {
            "0".to_string()
        } else {
            parts.join("")
        }
    }
}

/// Evaluate `e` at `var = x_val`.
fn eval_at(e: &Expr, var: &str, x_val: f64, ctx: &Context) -> f64 {
    let mut eval_ctx = ctx.clone();
    eval_ctx.set(var, x_val);
    eval(e, &eval_ctx).unwrap_or(f64::NAN)
}

/// Compute Taylor coefficients of `e` (viewed as a function of `var`) around
/// `var = a` by polynomial fitting. Evaluates `e` at `n_terms` points near `a`
/// (avoiding `a` itself) and solves the resulting Vandermonde system to get
/// coefficients `d_0, d_1, ..., d_{n_terms-1}` where
/// `e(a + t) ≈ d_0 + d_1*t + d_2*t^2 + ...`.
///
/// This avoids issues with removable singularities (e.g. `g(x) = x*(1/x)`)
/// because we never evaluate at `a` itself.
fn taylor_coeffs_by_fit(
    e: &Expr,
    var: &str,
    a: f64,
    n_terms: usize,
    ctx: &Context,
) -> Vec<f64> {
    // Use Chebyshev nodes on [-h, h] around a (none at 0).
    // Fit in scaled variable u = t/h, then convert back.
    let h = 1e-2;
    let n = n_terms;
    let mut nodes: Vec<f64> = Vec::with_capacity(n); // u values in [-1, 1]
    let mut values: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let u = ((2 * i + 1) as f64 / (2 * n) as f64 * std::f64::consts::PI).cos();
        let t = h * u;
        let x = a + t;
        let val = eval_at(e, var, x, ctx);
        if !val.is_finite() {
            // Fallback: use a non-Chebyshev point
            let u2 = (i as f64 + 1.0) / n as f64 * 2.0 - 1.0;
            let val2 = eval_at(e, var, a + h * u2, ctx);
            nodes.push(u2);
            values.push(val2);
        } else {
            nodes.push(u);
            values.push(val);
        }
    }

    // Solve Vandermonde in scaled variable u: V[i][j] = nodes[i]^j
    // Then d_n (Taylor coeff of t^n) = c_n / h^n where c_n is coeff of u^n.
    let mut aug = vec![vec![0.0f64; n + 1]; n];
    for i in 0..n {
        let mut pwr = 1.0;
        for j in 0..n {
            aug[i][j] = pwr;
            pwr *= nodes[i];
        }
        aug[i][n] = values[i];
    }

    // Gaussian elimination with partial pivoting
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in (col + 1)..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-300 {
            continue;
        }
        aug.swap(col, max_row);
        for row in (col + 1)..n {
            let factor = aug[row][col] / aug[col][col];
            for j in col..=n {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    // Back-substitution → c_n (coefficients in u)
    let mut c = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for j in (i + 1)..n {
            sum -= aug[i][j] * c[j];
        }
        if aug[i][i].abs() > 1e-300 {
            c[i] = sum / aug[i][i];
        }
    }

    // Convert: d_n = c_n / h^n (Taylor coefficient of (x-a)^n)
    let mut coeffs = vec![0.0f64; n];
    let mut h_pwr = 1.0;
    for i in 0..n {
        coeffs[i] = c[i] / h_pwr;
        h_pwr *= h;
    }
    coeffs
}

/// Compute the Laurent series of `f` around `x = a` with pole order `k`
/// and `n_positive` non-negative-power terms.
///
/// - `pole_order`: the order of the pole at `x = a` (0 = no singularity,
///   1 = simple pole, 2 = double pole, etc.)
/// - `n_positive`: number of non-negative-power terms to compute (like Taylor order)
pub fn laurent_series(
    f: &Expr,
    var: &str,
    a: f64,
    pole_order: usize,
    n_positive: usize,
) -> Result<LaurentSeries> {
    if pole_order == 0 && n_positive == 0 {
        return Err(MathError::InvalidArgument(
            "laurent_series: need at least one term".into(),
        ));
    }

    let ctx = Context::standard();

    // Principal part: expand g(x) = (x - a)^k * f(x) as Taylor series.
    // The first k+1 Taylor coefficients of g give c_{-k}, ..., c_0.
    // g(x) = Σ d_n (x-a)^n, so f(x) = g(x) / (x-a)^k = Σ d_n (x-a)^{n-k}
    // => c_{n-k} = d_n for n = 0, ..., k  (gives c_{-k}, ..., c_0)
    // Then the analytic part c_1, c_2, ... comes from Taylor of g beyond index k.

    let total_terms = pole_order + n_positive;
    if pole_order == 0 {
        // No pole — coefficients are Taylor coefficients of f
        let coeffs = taylor_coeffs_by_fit(f, var, a, total_terms, &ctx);
        return Ok(LaurentSeries {
            coeffs,
            pole_order: 0,
            center: a,
            var: var.to_string(),
        });
    }

    // Build g(x) = (x - a)^k * f(x) symbolically, but evaluate numerically.
    // g is analytic at a (removable singularity), so polynomial fitting works.
    let x_minus_a = if a == 0.0 {
        Expr::var(var)
    } else {
        Expr::sub(Expr::var(var), Expr::num(a))
    };
    let g = if pole_order == 1 {
        Expr::mul(x_minus_a, f.clone())
    } else {
        Expr::mul(Expr::pow(x_minus_a, Expr::num(pole_order as f64)), f.clone())
    };

    // Compute Taylor coefficients of g by polynomial fitting
    let d_coeffs = taylor_coeffs_by_fit(&g, var, a, total_terms, &ctx);

    // Laurent coefficients: c_{n-k} = d_n for n = 0, ..., total_terms-1
    // So coeffs[0] = c_{-k} = d_0, coeffs[1] = c_{-k+1} = d_1, etc.
    Ok(LaurentSeries {
        coeffs: d_coeffs,
        pole_order,
        center: a,
        var: var.to_string(),
    })
}

/// Compute the Laurent series from a string expression.
pub fn laurent_series_str(
    src: &str,
    var: &str,
    a: f64,
    pole_order: usize,
    n_positive: usize,
) -> Result<LaurentSeries> {
    let f = Parser::parse(src)?;
    laurent_series(&f, var, a, pole_order, n_positive)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn laurent_simple_pole() {
        // f(x) = 1/x, pole of order 1 at x=0
        // Laurent series: 1/x (exact, single term)
        let ls = laurent_series_str("1/x", "x", 0.0, 1, 3).unwrap();
        assert_eq!(ls.pole_order, 1);
        // c_{-1} = 1, c_0 = 0, c_1 = 0, c_2 = 0
        assert!(close(ls.coeff(-1), 1.0, 1e-6), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 0.0, 1e-6), "c_0 = {}", ls.coeff(0));
        assert!(close(ls.coeff(1), 0.0, 1e-6), "c_1 = {}", ls.coeff(1));
    }

    #[test]
    fn laurent_double_pole() {
        // f(x) = 1/x^2, pole of order 2 at x=0
        let ls = laurent_series_str("1/x^2", "x", 0.0, 2, 3).unwrap();
        assert_eq!(ls.pole_order, 2);
        assert!(close(ls.coeff(-2), 1.0, 1e-4), "c_-2 = {}", ls.coeff(-2));
        assert!(close(ls.coeff(-1), 0.0, 1e-4), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 0.0, 1e-4), "c_0 = {}", ls.coeff(0));
    }

    #[test]
    fn laurent_simple_pole_with_analytic() {
        // f(x) = 1/x + 1 + x, pole of order 1 at x=0
        let ls = laurent_series_str("1/x + 1 + x", "x", 0.0, 1, 3).unwrap();
        assert!(close(ls.coeff(-1), 1.0, 1e-6), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 1.0, 1e-6), "c_0 = {}", ls.coeff(0));
        assert!(close(ls.coeff(1), 1.0, 1e-6), "c_1 = {}", ls.coeff(1));
    }

    #[test]
    fn laurent_eval_simple_pole() {
        // f(x) = 2/x + 3, eval at x=2 should give 2/2 + 3 = 4
        let ls = laurent_series_str("2/x + 3", "x", 0.0, 1, 2).unwrap();
        let val = ls.eval(2.0);
        assert!(close(val, 4.0, 1e-6), "eval at 2: got {} want 4", val);
    }

    #[test]
    fn laurent_eval_double_pole() {
        // f(x) = 1/x^2 + 1/x + x^2
        // eval at x=2: 1/4 + 1/2 + 4 = 4.75
        let ls = laurent_series_str("1/x^2 + 1/x + x^2", "x", 0.0, 2, 3).unwrap();
        let val = ls.eval(2.0);
        assert!(close(val, 4.75, 1e-4), "eval at 2: got {} want 4.75", val);
    }

    #[test]
    fn laurent_no_pole() {
        // f(x) = exp(x), no pole → should match Taylor
        let ls = laurent_series_str("exp(x)", "x", 0.0, 0, 5).unwrap();
        assert_eq!(ls.pole_order, 0);
        // c_0 = 1, c_1 = 1, c_2 = 1/2
        assert!(close(ls.coeff(0), 1.0, 1e-6), "c_0 = {}", ls.coeff(0));
        assert!(close(ls.coeff(1), 1.0, 1e-6), "c_1 = {}", ls.coeff(1));
        assert!(close(ls.coeff(2), 0.5, 1e-6), "c_2 = {}", ls.coeff(2));
    }

    #[test]
    fn laurent_eval_no_pole() {
        // f(x) = sin(x), no pole, 7 terms
        let ls = laurent_series_str("sin(x)", "x", 0.0, 0, 7).unwrap();
        for &x in &[0.1, 0.5, 1.0] {
            let approx = ls.eval(x);
            let exact = x.sin();
            assert!(close(approx, exact, 1e-3), "at x={}: got {} want {}", x, approx, exact);
        }
    }

    #[test]
    fn laurent_exp_over_x() {
        // f(x) = exp(x)/x = (1/x)(1 + x + x^2/2 + ...)
        // = 1/x + 1 + x/2 + x^2/6 + ...
        let ls = laurent_series_str("exp(x)/x", "x", 0.0, 1, 5).unwrap();
        assert!(close(ls.coeff(-1), 1.0, 1e-6), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 1.0, 1e-6), "c_0 = {}", ls.coeff(0));
        assert!(close(ls.coeff(1), 0.5, 1e-6), "c_1 = {}", ls.coeff(1));
        assert!(close(ls.coeff(2), 1.0 / 6.0, 1e-6), "c_2 = {}", ls.coeff(2));
    }

    #[test]
    fn laurent_eval_exp_over_x() {
        // f(x) = exp(x)/x, eval at x=1: e/1 = e ≈ 2.71828
        // Use 6 positive terms (7×7 Vandermonde is well-conditioned)
        let ls = laurent_series_str("exp(x)/x", "x", 0.0, 1, 6).unwrap();
        let val = ls.eval(1.0);
        assert!(close(val, std::f64::consts::E, 1e-2), "eval at 1: got {} want e", val);
    }

    #[test]
    fn laurent_around_nonzero() {
        // f(x) = 1/(x-2), pole of order 1 at x=2
        let ls = laurent_series_str("1/(x-2)", "x", 2.0, 1, 3).unwrap();
        assert!(close(ls.coeff(-1), 1.0, 1e-6), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 0.0, 1e-6), "c_0 = {}", ls.coeff(0));
        // eval at x=3: 1/(3-2) = 1
        let val = ls.eval(3.0);
        assert!(close(val, 1.0, 1e-6), "eval at 3: got {} want 1", val);
    }

    #[test]
    fn laurent_to_string() {
        let ls = laurent_series_str("1/x + 2 + x", "x", 0.0, 1, 2).unwrap();
        let s = ls.to_string();
        assert!(s.contains("1/x"), "string should contain 1/x: {}", s);
    }

    #[test]
    fn laurent_invalid_args() {
        assert!(laurent_series_str("1/x", "x", 0.0, 0, 0).is_err());
    }

    #[test]
    fn laurent_rational_function() {
        // f(x) = 1/(x(1-x)) = 1/x + 1 + x + x^2 + ...  (partial fractions: 1/x + 1/(1-x))
        // Pole of order 1 at x=0
        let ls = laurent_series_str("1/(x*(1-x))", "x", 0.0, 1, 5).unwrap();
        assert!(close(ls.coeff(-1), 1.0, 1e-6), "c_-1 = {}", ls.coeff(-1));
        assert!(close(ls.coeff(0), 1.0, 1e-6), "c_0 = {}", ls.coeff(0));
        assert!(close(ls.coeff(1), 1.0, 1e-6), "c_1 = {}", ls.coeff(1));
        assert!(close(ls.coeff(2), 1.0, 1e-6), "c_2 = {}", ls.coeff(2));
    }
}
