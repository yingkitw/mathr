//! Equation solvers: Newton–Raphson, bisection, secant, and polynomial
//! root finding via the Durand–Kerner method.

use crate::error::{MathError, Result};

/// Configuration shared by all iterative solvers.
#[derive(Debug, Clone)]
pub struct SolveOptions {
    pub max_iter: usize,
    pub tol: f64,
    pub h: f64, // finite-difference step
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-10,
            h: 1e-6,
        }
    }
}

/// Bisection method: finds a root of `f` on `[a, b]` assuming `f(a)*f(b) < 0`.
pub fn bisect<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut lo = a;
    let mut hi = b;
    let fa = f(lo);
    let fb = f(hi);
    if fa == 0.0 {
        return Ok((lo, 0.0));
    }
    if fb == 0.0 {
        return Ok((hi, 0.0));
    }
    if fa.signum() == fb.signum() {
        return Err(MathError::InvalidArgument(format!(
            "bisect: f(a) and f(b) must have opposite signs (got f({})={}, f({})={})",
            lo, fa, hi, fb
        )));
    }
    let mut mid;
    for _ in 0..opts.max_iter {
        mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if fm.abs() < opts.tol || (hi - lo) * 0.5 < opts.tol {
            return Ok((mid, fm));
        }
        if fm.signum() == fa.signum() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(MathError::NotConvergent(format!(
        "bisect did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Newton–Raphson method. Requires a derivative `df`; if you don't have one,
/// use [`newton_central`] which approximates it with a centred difference.
pub fn newton<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(
    f: F,
    df: D,
    x0: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut x = x0;
    let mut fx = f(x);
    for i in 0..opts.max_iter {
        if fx.abs() < opts.tol {
            return Ok((x, fx));
        }
        let dfx = df(x);
        if dfx == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "newton: derivative vanished at iteration {}",
                i
            )));
        }
        let step = fx / dfx;
        x -= step;
        fx = f(x);
        if step.abs() < opts.tol {
            return Ok((x, fx));
        }
    }
    Err(MathError::NotConvergent(format!(
        "newton did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Newton method with a numerically-estimated derivative (central difference).
pub fn newton_central<F>(f: F, x0: f64, opts: SolveOptions) -> Result<(f64, f64)>
where
    F: Fn(f64) -> f64 + Clone,
{
    let h = opts.h;
    let f2 = f.clone();
    let df = move |x: f64| (f2(x + h) - f2(x - h)) / (2.0 * h);
    newton(f, df, x0, opts)
}

/// Secant method: derivative-free, but starts with two initial guesses.
pub fn secant<F: Fn(f64) -> f64>(
    f: F,
    x0: f64,
    x1: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut x_prev = x0;
    let mut x_curr = x1;
    let mut f_prev = f(x_prev);
    let mut f_curr = f(x_curr);
    for _ in 0..opts.max_iter {
        if f_curr.abs() < opts.tol {
            return Ok((x_curr, f_curr));
        }
        let denom = f_curr - f_prev;
        if denom.abs() < 1e-30 {
            return Err(MathError::NotConvergent(
                "secant: zero denominator".into(),
            ));
        }
        let x_next = x_curr - f_curr * (x_curr - x_prev) / denom;
        if (x_next - x_curr).abs() < opts.tol {
            let fn_ = f(x_next);
            return Ok((x_next, fn_));
        }
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
        f_curr = f(x_next);
    }
    Err(MathError::NotConvergent(format!(
        "secant did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Find all (real) roots of a polynomial with real coefficients given by
/// `coeffs` in descending order: `coeffs[0]*x^n + ... + coeffs[n]`.
///
/// Implemented with the Durand–Kerner method, which converges quickly for
/// well-separated roots and is robust enough for moderate degrees.
pub fn polynomial_roots(coeffs: &[f64]) -> Result<Vec<(f64, f64)>> {
    let n = coeffs.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    // normalize so leading coefficient is 1
    let lead = coeffs[0];
    if lead == 0.0 {
        return Err(MathError::InvalidArgument("polynomial leading coefficient is 0".into()));
    }
    let a: Vec<f64> = coeffs.iter().map(|c| c / lead).collect();

    // initial guesses: complex roots spread around the unit circle
    use crate::complex::Complex;
    let mut z: Vec<Complex<f64>> = (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / n as f64;
            Complex::from_polar(0.4 + 0.9 * (n as f64 / 10.0), theta)
        })
        .collect();

    for _ in 0..200 {
        let mut converged = true;
        for i in 0..n {
            let mut p = Complex::new(a[0], 0.0);
            for j in 1..a.len() {
                p = p * z[i] + Complex::new(a[j], 0.0);
            }
            let mut denom = Complex::new(1.0, 0.0);
            for j in 0..n {
                if i != j {
                    denom = denom * (z[i] - z[j]);
                }
            }
            let step = if denom.abs() < 1e-30 {
                Complex::new(0.0, 0.0)
            } else {
                p / denom
            };
            z[i] = z[i] - step;
            if step.abs() > 1e-10 {
                converged = false;
            }
        }
        if converged {
            break;
        }
    }

    // collect real roots (within tolerance) and pair them with f(root)
    let mut out = Vec::new();
    for zi in z {
        if zi.im.abs() < 1e-6 {
            let x = zi.re;
            let fx = polynomial_eval(coeffs, x);
            out.push((x, fx));
        }
    }
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    Ok(out)
}

fn polynomial_eval(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    for c in coeffs {
        acc = acc * x + c;
    }
    acc
}

/// Solve a system of `n` nonlinear equations in `n` unknowns using
/// Newton's method with a finite-difference Jacobian.
///
/// `system` is a function that maps the current guess `x ∈ ℝⁿ` to
/// the residual vector `f(x) ∈ ℝⁿ`. The algorithm iterates
/// `x ← x − J(x)⁻¹ f(x)` until either `‖f(x)‖ < tol` or `max_iter` is hit.
///
/// On non-convergence or singular Jacobian, returns [`MathError::NotConvergent`].
pub fn newton_system<F: Fn(&[f64]) -> Vec<f64>>(
    system: F,
    x0: &[f64],
    opts: SolveOptions,
) -> Result<Vec<f64>> {
    let n = x0.len();
    if n == 0 {
        return Err(MathError::InvalidArgument("newton_system: empty initial guess".into()));
    }
    let h = if opts.h <= 0.0 { 1e-6 } else { opts.h };
    let mut x = x0.to_vec();
    for _ in 0..opts.max_iter {
        let f0 = system(&x);
        if f0.len() != n {
            return Err(MathError::InvalidArgument(format!(
                "newton_system: system returned {} values, expected {}",
                f0.len(), n
            )));
        }
        let norm: f64 = f0.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < opts.tol {
            return Ok(x);
        }
        // Build Jacobian by central finite differences.
        let mut jac = vec![0.0_f64; n * n];
        for j in 0..n {
            let mut xp = x.clone();
            let mut xm = x.clone();
            xp[j] += h;
            xm[j] -= h;
            let fp = system(&xp);
            let fm = system(&xm);
            for i in 0..n {
                jac[i * n + j] = (fp[i] - fm[i]) / (2.0 * h);
            }
        }
        // Solve J · delta = -f0 for delta.
        let neg_f0: Vec<f64> = f0.iter().map(|v| -v).collect();
        let m = crate::matrix::Matrix::from_row_major(n, n, jac.clone())
            .map_err(|e| MathError::InvalidArgument(format!("newton_system: {}", e)))?;
        let delta = m
            .solve(&neg_f0)
            .map_err(|_| MathError::NotConvergent("newton_system: singular Jacobian".into()))?;
        for i in 0..n {
            x[i] += delta[i];
        }
    }
    Err(MathError::NotConvergent(format!(
        "newton_system: failed to converge in {} iterations",
        opts.max_iter
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn bisect_root() {
        // x^2 - 4 = 0 on [0, 5]
        let f = |x: f64| x * x - 4.0;
        let (r, _) = bisect(f, 0.0, 5.0, SolveOptions::default()).unwrap();
        assert!((r - 2.0).abs() < 1e-6);
    }

    #[test]
    fn newton_root() {
        // cos(x) = x  -> Newton's with f=cos(x)-x, df=-sin(x)-1
        let f = |x: f64| x.cos() - x;
        let df = |x: f64| -x.sin() - 1.0;
        let (r, _) = newton(f, df, 0.5, SolveOptions::default()).unwrap();
        // known fixed point ≈ 0.7390851332151607
        assert!((r - 0.7390851332151607).abs() < 1e-8);
    }

    #[test]
    fn newton_central_root() {
        // Newton from 3 on sin(x) converges to the closest root, which is π.
        let f = |x: f64| x.sin();
        let (r, _) = newton_central(f, 3.0, SolveOptions::default()).unwrap();
        assert!((r - PI).abs() < 1e-7, "expected root at π, got {}", r);
    }

    #[test]
    fn secant_root() {
        let f = |x: f64| x.powi(3) - 2.0;
        let (r, _) = secant(f, 1.0, 2.0, SolveOptions::default()).unwrap();
        assert!((r - 2f64.cbrt()).abs() < 1e-8);
    }

    #[test]
    fn polynomial_two_roots() {
        // x^2 - 5x + 6 = (x-2)(x-3)
        let r = polynomial_roots(&[1.0, -5.0, 6.0]).unwrap();
        let xs: Vec<f64> = r.iter().map(|(x, _)| *x).collect();
        assert_eq!(xs.len(), 2);
        let mut sorted = xs.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 2.0).abs() < 1e-4);
        assert!((sorted[1] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn newton_system_linear() {
        // 2x + y = 5; x + 3y = 7  →  (x, y) = (1.6, 1.8).
        let system = |x: &[f64]| vec![2.0 * x[0] + x[1] - 5.0, x[0] + 3.0 * x[1] - 7.0];
        let sol = newton_system(system, &[0.0, 0.0], SolveOptions::default()).unwrap();
        assert!((sol[0] - 1.6).abs() < 1e-6, "x = {}", sol[0]);
        assert!((sol[1] - 1.8).abs() < 1e-6, "y = {}", sol[1]);
    }

    #[test]
    fn newton_system_nonlinear() {
        // Intersection of x² + y² = 1 and y = x² − 0.5.
        let system = |x: &[f64]| {
            vec![
                x[0] * x[0] + x[1] * x[1] - 1.0,
                x[1] - x[0] * x[0] + 0.5,
            ]
        };
        let sol = newton_system(system, &[0.9, 0.3], SolveOptions::default()).unwrap();
        let r = (sol[0] * sol[0] + sol[1] * sol[1] - 1.0).abs();
        assert!(r < 1e-6, "residual = {}", r);
        let s = (sol[1] - sol[0] * sol[0] + 0.5).abs();
        assert!(s < 1e-6, "residual = {}", s);
    }

    #[test]
    fn newton_system_3d() {
        // 3×3 linear: x + y + z = 6, 2x - y + z = 3, x + 2y - z = 2.
        // Solution: x = 1, y = 2, z = 3.
        let system = |v: &[f64]| {
            vec![
                v[0] + v[1] + v[2] - 6.0,
                2.0 * v[0] - v[1] + v[2] - 3.0,
                v[0] + 2.0 * v[1] - v[2] - 2.0,
            ]
        };
        let sol = newton_system(system, &[0.0, 0.0, 0.0], SolveOptions::default()).unwrap();
        assert!((sol[0] - 1.0).abs() < 1e-6, "x = {}", sol[0]);
        assert!((sol[1] - 2.0).abs() < 1e-6, "y = {}", sol[1]);
        assert!((sol[2] - 3.0).abs() < 1e-6, "z = {}", sol[2]);
    }
}