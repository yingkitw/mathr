//! Numerical calculus: differentiation, integration, vector calculus helpers.

use crate::error::{MathError, Result};

/// Numerical derivative of `f` at `x` using the centred five-point stencil:
///
/// `f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)`
///
/// The five-point stencil is `O(h^4)` accurate, much better than the
/// naive `O(h)` forward difference.
pub fn derivative<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    let h = if h <= 0.0 { 1e-5 } else { h };
    let twoh = 2.0 * h;
    (-f(x + twoh) + 8.0 * f(x + h) - 8.0 * f(x - h) + f(x - twoh)) / (12.0 * h)
}

/// Second derivative using the five-point stencil:
/// `f''(x) ≈ (-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) - f(x-2h)) / (12h^2)`
pub fn second_derivative<F: Fn(f64) -> f64>(f: F, x: f64, h: f64) -> f64 {
    let h = if h <= 0.0 { 1e-3 } else { h };
    let twoh = 2.0 * h;
    let h2 = 12.0 * h * h;
    (-f(x + twoh) + 16.0 * f(x + h) - 30.0 * f(x) + 16.0 * f(x - h) - f(x - twoh)) / h2
}

/// Trapezoidal-rule definite integral of `f` over `[a, b]` with `n` panels.
pub fn integrate_trap<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += f(x);
    }
    sum * h
}

/// Simpson's 1/3 rule for definite integrals over `[a, b]` using `n` panels.
/// `n` must be even.
pub fn integrate_simpson<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> Result<f64> {
    if n == 0 || n % 2 != 0 {
        return Err(MathError::InvalidArgument(
            "Simpson's rule requires n > 0 and even".into(),
        ));
    }
    let h = (b - a) / n as f64;
    let mut sum = f(a) + f(b);
    for i in 1..n {
        let x = a + i as f64 * h;
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(x);
    }
    Ok(sum * h / 3.0)
}

/// Adaptive Simpson's quadrature: recursively subdivides the interval until
/// the estimate meets the requested absolute/relative tolerance.
pub fn integrate_adaptive<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    max_depth: u32,
) -> Result<f64> {
    if !(a.is_finite() && b.is_finite() && a < b) {
        return Err(MathError::InvalidArgument(format!(
            "integrate_adaptive: bad interval [{}, {}]",
            a, b
        )));
    }

    fn rec<F2: Fn(f64) -> f64>(
        f: &F2,
        a: f64,
        b: f64,
        fa: f64,
        fb: f64,
        whole: f64,
        tol: f64,
        depth: u32,
        max_depth: u32,
    ) -> f64 {
        let m = 0.5 * (a + b);
        let fm = f(m);
        let lm = (a + m) * 0.5;
        let rm = (m + b) * 0.5;
        let flm = f(lm);
        let frm = f(rm);
        let left = (fa + 4.0 * flm + fm) * (m - a) / 6.0;
        let right = (fm + 4.0 * frm + fb) * (b - m) / 6.0;
        let refined = left + right;
        let err = (refined - whole).abs();
        if depth >= max_depth || err < 15.0 * tol {
            refined + (refined - whole) / 15.0
        } else {
            rec(f, a, m, fa, fm, left, tol * 0.5, depth + 1, max_depth)
                + rec(f, m, b, fm, fb, right, tol * 0.5, depth + 1, max_depth)
        }
    }

    let fa = f(a);
    let fb = f(b);
    let m = 0.5 * (a + b);
    let fm = f(m);
    let whole = (fa + 4.0 * fm + fb) * (b - a) / 6.0;
    Ok(rec(&f, a, b, fa, fb, whole, tol, 0, max_depth))
}

/// Numerically compute a partial derivative with respect to one variable,
/// holding other variables fixed at the supplied values.
pub fn partial<F: Fn(&[f64]) -> f64>(f: F, point: &[f64], var: usize, h: f64) -> Result<f64> {
    if var >= point.len() {
        return Err(MathError::InvalidArgument(format!(
            "partial: variable index {} out of range",
            var
        )));
    }
    let h = if h <= 0.0 { 1e-5 } else { h };
    let twoh = 2.0 * h;
    let mut plus = point.to_vec();
    plus[var] += twoh;
    let mut p1 = plus.clone();
    p1[var] -= h; // actually we'll do two stencils separately below
    let _ = (plus, p1);

    // 5-point centred stencil in the requested coordinate.
    let mut p = point.to_vec();
    let offsets = [(-2.0 * h, 1.0), (-h, -8.0), (h, 8.0), (2.0 * h, -1.0)];
    let mut acc = 0.0;
    for (dx, coef) in offsets {
        p[var] = point[var] + dx;
        acc += coef * f(&p);
    }
    Ok(acc / (12.0 * h))
}

/// Numerical gradient of `f` at `point`.
pub fn gradient<F: Fn(&[f64]) -> f64>(f: F, point: &[f64], h: f64) -> Result<Vec<f64>> {
    let mut out = Vec::with_capacity(point.len());
    for i in 0..point.len() {
        out.push(partial(&f, point, i, h)?);
    }
    Ok(out)
}

/// Romberg integration: starts from the trapezoidal rule with `n = 1` panel
/// and applies Richardson extrapolation to build a triangular table whose
/// entries have successively higher orders of accuracy in `h²`.
///
/// Returns the best (highest-level) estimate. `levels` controls the table
/// size: a value of `k` gives an integrator that is exact for polynomials
/// up to degree `2^k - 1` (i.e. `O(h^(2k))` for smooth integrands).
pub fn integrate_romberg<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, levels: usize) -> Result<f64> {
    if !(a.is_finite() && b.is_finite() && a < b) {
        return Err(MathError::InvalidArgument(format!(
            "integrate_romberg: bad interval [{}, {}]",
            a, b
        )));
    }
    if levels == 0 {
        return Err(MathError::InvalidArgument(
            "integrate_romberg: levels must be >= 1".into(),
        ));
    }
    let levels = levels.min(20);
    // Compute the first column: trapezoidal rule with increasing panel counts.
    // We track the trapezoidal sum `s_n` such that `T_n = h_n · s_n`.  The
    // doubling recurrence for the sum (NOT pre-multiplied by h) is just
    //   s_{2n} = s_n + Σ_{k=0}^{n-1} f(a + (2k+1) · h_new),
    // where `h_new = (b-a)/(2n)` is the new panel width.
    let mut col: Vec<f64> = Vec::with_capacity(levels);
    let mut n = 1usize;
    let h = (b - a) / n as f64;
    let mut sum = 0.5 * (f(a) + f(b));
    col.push(sum * h);
    for _ in 1..levels {
        n *= 2;
        let h_new = (b - a) / n as f64;
        let mut add = 0.0;
        let mut x = a + h_new;
        for _ in 0..(n / 2) {
            add += f(x);
            x += h_new * 2.0;
        }
        sum += add;
        col.push(sum * h_new);
    }
    // Richardson extrapolation: R[i][j] = (4^j R[i][j-1] - R[i-1][j-1]) / (4^j - 1).
    let mut table: Vec<Vec<f64>> = vec![col];
    for j in 1..levels {
        let prev = &table[j - 1];
        let mut next = Vec::with_capacity(levels - j);
        let four_j = (4f64).powi(j as i32);
        for i in 0..(levels - j) {
            let v = (four_j * prev[i + 1] - prev[i]) / (four_j - 1.0);
            next.push(v);
        }
        table.push(next);
    }
    Ok(table[levels - 1][0])
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn derivative_of_quadratic() {
        // f(x) = x^2, f'(x) = 2x
        let f = |x: f64| x * x;
        let d = derivative(f, 3.0, 1e-4);
        assert!((d - 6.0).abs() < 1e-6);
    }

    #[test]
    fn derivative_of_sin() {
        let f = |x: f64| x.sin();
        let d = derivative(f, PI / 2.0, 1e-4);
        assert!((d - 0.0).abs() < 1e-6);
    }

    #[test]
    fn second_derivative_works() {
        // f(x) = x^3, f''(x) = 6x
        let f = |x: f64| x.powi(3);
        let d = second_derivative(f, 2.0, 1e-2);
        assert!((d - 12.0).abs() < 1e-3);
    }

    #[test]
    fn trapezoidal_known_integral() {
        // ∫_0^π sin(x) dx = 2
        let f = |x: f64| x.sin();
        let v = integrate_trap(f, 0.0, PI, 1000);
        assert!((v - 2.0).abs() < 1e-3);
    }

    #[test]
    fn simpson_known_integral() {
        let f = |x: f64| x.sin();
        let v = integrate_simpson(f, 0.0, PI, 1000).unwrap();
        assert!((v - 2.0).abs() < 1e-6);
    }

    #[test]
    fn adaptive_integral_gauss() {
        // ∫_{-∞}^{∞} exp(-x^2) dx = sqrt(pi).  But our adaptive integrator is
        // for finite intervals; instead test ∫_0^1 exp(-x^2) dx ≈ 0.746824.
        let f = |x: f64| (-x * x).exp();
        let v = integrate_adaptive(f, 0.0, 1.0, 1e-9, 30).unwrap();
        assert!((v - 0.7468241328).abs() < 1e-7);
    }

    #[test]
    fn gradient_simple() {
        // f(x, y) = x^2 + y^2 -> grad = (2x, 2y)
        let f = |p: &[f64]| p[0] * p[0] + p[1] * p[1];
        let g = gradient(f, &[3.0, 4.0], 1e-4).unwrap();
        assert!((g[0] - 6.0).abs() < 1e-5);
        assert!((g[1] - 8.0).abs() < 1e-5);
    }

    #[test]
    fn romberg_sin() {
        // ∫_0^π sin(x) dx = 2
        let f = |x: f64| x.sin();
        let v = integrate_romberg(f, 0.0, PI, 8).unwrap();
        assert!((v - 2.0).abs() < 1e-12, "got {}", v);
    }

    #[test]
    fn romberg_polynomial_exact() {
        // Romberg of level k is exact for degree ≤ 2^k − 1.
        // f(x) = x^5 + 3x^3 − 2x + 7 over [0, 1] should be exact at level 4.
        let f = |x: f64| x.powi(5) + 3.0 * x.powi(3) - 2.0 * x + 7.0;
        let analytical = 1.0 / 6.0 + 3.0 / 4.0 - 1.0 + 7.0;
        let v = integrate_romberg(f, 0.0, 1.0, 4).unwrap();
        assert!((v - analytical).abs() < 1e-12, "got {} want {}", v, analytical);
    }

    #[test]
    fn romberg_gaussian() {
        // ∫_0^1 exp(-x²) dx ≈ 0.7468241328124270
        let f = |x: f64| (-x * x).exp();
        let v = integrate_romberg(f, 0.0, 1.0, 10).unwrap();
        assert!((v - 0.7468241328124270).abs() < 1e-10, "got {}", v);
    }
}