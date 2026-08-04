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

/// Fourier series coefficients for a function on `[-L, L]`.
///
/// The series approximation is:
/// `f(x) ≈ a₀/2 + Σ_{n=1}^{N} [aₙ cos(nπx/L) + bₙ sin(nπx/L)]`
#[derive(Debug, Clone)]
pub struct FourierSeries {
    pub a0: f64,
    pub an: Vec<f64>,
    pub bn: Vec<f64>,
    pub l: f64,
}

/// Compute the Fourier series of `f` on `[-L, L]` with `n_terms` harmonics.
///
/// Uses Simpson's rule with a high panel count for accurate coefficient
/// integration.
pub fn fourier_series<F: Fn(f64) -> f64>(
    f: F,
    n_terms: usize,
    l: f64,
) -> Result<FourierSeries> {
    if l <= 0.0 {
        return Err(MathError::InvalidArgument(
            "fourier_series: L must be positive".into(),
        ));
    }
    if n_terms == 0 {
        return Err(MathError::InvalidArgument(
            "fourier_series: n_terms must be >= 1".into(),
        ));
    }

    let panels = 2000.max(n_terms * 40);
    let pi_over_l = std::f64::consts::PI / l;

    let a0 = integrate_simpson(&f, -l, l, panels)? / l;
    let mut an = Vec::with_capacity(n_terms);
    let mut bn = Vec::with_capacity(n_terms);

    for n in 1..=n_terms {
        let n_f = n as f64;
        let omega = n_f * pi_over_l;
        let an_val = integrate_simpson(
            |x| f(x) * (omega * x).cos(),
            -l,
            l,
            panels,
        )? / l;
        let bn_val = integrate_simpson(
            |x| f(x) * (omega * x).sin(),
            -l,
            l,
            panels,
        )? / l;
        an.push(an_val);
        bn.push(bn_val);
    }

    Ok(FourierSeries { a0, an, bn, l })
}

/// Evaluate a Fourier series at point `x`.
pub fn fourier_eval(series: &FourierSeries, x: f64) -> f64 {
    let pi_over_l = std::f64::consts::PI / series.l;
    let mut sum = series.a0 / 2.0;
    for (i, (a, b)) in series.an.iter().zip(series.bn.iter()).enumerate() {
        let n = (i + 1) as f64;
        let omega = n * pi_over_l;
        sum += a * (omega * x).cos() + b * (omega * x).sin();
    }
    sum
}

/// Simple deterministic linear congruential generator for reproducible
/// Monte Carlo sampling without external dependencies.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg { state: seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407) }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    /// Uniform sample in `[0, 1)`.
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Monte Carlo integration of a 1-D function over `[a, b]` using `n_samples`
/// uniform random points and a fixed `seed` for reproducibility.
///
/// Returns `(estimate, standard_error)` where the standard error is the
/// Monte Carlo standard deviation divided by `sqrt(n)`.
pub fn monte_carlo_integrate_1d<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    n_samples: usize,
    seed: u64,
) -> Result<(f64, f64)> {
    if !(a.is_finite() && b.is_finite() && a < b) {
        return Err(MathError::InvalidArgument(format!(
            "monte_carlo_integrate_1d: bad interval [{}, {}]",
            a, b
        )));
    }
    if n_samples == 0 {
        return Err(MathError::InvalidArgument(
            "monte_carlo_integrate_1d: n_samples must be > 0".into(),
        ));
    }
    let mut rng = Lcg::new(seed);
    let width = b - a;
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    for _ in 0..n_samples {
        let x = a + rng.next_f64() * width;
        let fx = f(x);
        sum += fx;
        sum_sq += fx * fx;
    }
    let mean = sum / n_samples as f64;
    let mean_sq = sum_sq / n_samples as f64;
    let variance = (mean_sq - mean * mean).max(0.0);
    let se = (variance / n_samples as f64).sqrt();
    let estimate = mean * width;
    Ok((estimate, se))
}

/// Monte Carlo integration of an N-D function over a hyperrectangle
/// defined by `bounds` (each element is `(lo, hi)` for that dimension).
///
/// Returns `(estimate, standard_error)`.
pub fn monte_carlo_integrate_nd<F: Fn(&[f64]) -> f64>(
    f: F,
    bounds: &[(f64, f64)],
    n_samples: usize,
    seed: u64,
) -> Result<(f64, f64)> {
    if bounds.is_empty() {
        return Err(MathError::InvalidArgument(
            "monte_carlo_integrate_nd: need at least one dimension".into(),
        ));
    }
    if n_samples == 0 {
        return Err(MathError::InvalidArgument(
            "monte_carlo_integrate_nd: n_samples must be > 0".into(),
        ));
    }
    for &(lo, hi) in bounds {
        if !(lo.is_finite() && hi.is_finite() && lo < hi) {
            return Err(MathError::InvalidArgument(format!(
                "monte_carlo_integrate_nd: bad bounds [{}, {}]",
                lo, hi
            )));
        }
    }
    let dim = bounds.len();
    let mut rng = Lcg::new(seed);
    let volume: f64 = bounds.iter().map(|(lo, hi)| hi - lo).product();
    let mut sum = 0.0;
    let mut sum_sq = 0.0;
    let mut point = vec![0.0; dim];
    for _ in 0..n_samples {
        for (i, &(lo, hi)) in bounds.iter().enumerate() {
            point[i] = lo + rng.next_f64() * (hi - lo);
        }
        let fx = f(&point);
        sum += fx;
        sum_sq += fx * fx;
    }
    let mean = sum / n_samples as f64;
    let mean_sq = sum_sq / n_samples as f64;
    let variance = (mean_sq - mean * mean).max(0.0);
    let se = (variance / n_samples as f64).sqrt();
    Ok((mean * volume, se))
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

    #[test]
    fn fourier_square_wave() {
        // Square wave on [-1, 1]: f(x) = 1 for x > 0, -1 for x < 0
        // Known coefficients: a0 = 0, an = 0, bn = (2/(nπ))(1 - cos(nπ))
        // For odd n: bn = 4/(nπ)
        // Note: Simpson's rule converges slowly at the discontinuity, so
        // we use relaxed tolerances for a0 and an.
        let f = |x: f64| if x > 0.0 { 1.0 } else { -1.0 };
        let fs = fourier_series(f, 10, 1.0).unwrap();
        assert!(fs.a0.abs() < 1e-2, "a0 should be ~0, got {}", fs.a0);
        for n in 0..10 {
            let n1 = (n + 1) as f64;
            assert!(fs.an[n].abs() < 1e-2, "a{} should be ~0, got {}", n + 1, fs.an[n]);
            if n1 % 2.0 == 1.0 {
                let expected = 4.0 / (n1 * PI);
                assert!((fs.bn[n] - expected).abs() < 1e-3,
                    "b{} should be ~{}, got {}", n + 1, expected, fs.bn[n]);
            } else {
                assert!(fs.bn[n].abs() < 1e-3, "b{} should be ~0, got {}", n + 1, fs.bn[n]);
            }
        }
    }

    #[test]
    fn fourier_cosine() {
        // f(x) = cos(πx) on [-1, 1] — this is already a Fourier mode.
        // a1 = 1, all other an = 0, all bn = 0, a0 = 0.
        let f = |x: f64| (PI * x).cos();
        let fs = fourier_series(f, 5, 1.0).unwrap();
        assert!(fs.a0.abs() < 1e-6, "a0 should be ~0, got {}", fs.a0);
        assert!((fs.an[0] - 1.0).abs() < 1e-6, "a1 should be 1, got {}", fs.an[0]);
        for n in 1..5 {
            assert!(fs.an[n].abs() < 1e-6, "a{} should be ~0, got {}", n + 1, fs.an[n]);
            assert!(fs.bn[n].abs() < 1e-6, "b{} should be ~0, got {}", n + 1, fs.bn[n]);
        }
        assert!(fs.bn[0].abs() < 1e-6, "b1 should be ~0, got {}", fs.bn[0]);
    }

    #[test]
    fn fourier_eval_accuracy() {
        // f(x) = x on [-1, 1]. Fourier series: bn = 2*(-1)^(n+1)/(nπ)
        // With enough terms, evaluation at x=0.5 should be close to 0.5.
        let f = |x: f64| x;
        let fs = fourier_series(f, 20, 1.0).unwrap();
        let val = fourier_eval(&fs, 0.5);
        assert!((val - 0.5).abs() < 0.05, "fourier_eval at 0.5: got {} want ~0.5", val);
    }

    #[test]
    fn fourier_constant() {
        // f(x) = 3 on [-1, 1]. a0 = 6, all other coefficients = 0.
        let f = |_x: f64| 3.0;
        let fs = fourier_series(f, 5, 1.0).unwrap();
        assert!((fs.a0 - 6.0).abs() < 1e-6, "a0 should be 6, got {}", fs.a0);
        for n in 0..5 {
            assert!(fs.an[n].abs() < 1e-6, "a{} should be ~0, got {}", n + 1, fs.an[n]);
            assert!(fs.bn[n].abs() < 1e-6, "b{} should be ~0, got {}", n + 1, fs.bn[n]);
        }
        let val = fourier_eval(&fs, 0.3);
        assert!((val - 3.0).abs() < 1e-6, "eval should be 3, got {}", val);
    }

    #[test]
    fn fourier_invalid_args() {
        let f = |x: f64| x;
        assert!(fourier_series(f, 5, -1.0).is_err());
        assert!(fourier_series(f, 0, 1.0).is_err());
    }

    #[test]
    fn mc_1d_constant() {
        // ∫_0^1 5 dx = 5, SE should be 0
        let (est, se) = monte_carlo_integrate_1d(|_| 5.0, 0.0, 1.0, 10000, 42).unwrap();
        assert!((est - 5.0).abs() < 1e-10, "got {}", est);
        assert!(se < 1e-10, "se should be ~0, got {}", se);
    }

    #[test]
    fn mc_1d_linear() {
        // ∫_0^1 x dx = 0.5
        let (est, _se) = monte_carlo_integrate_1d(|x| x, 0.0, 1.0, 100000, 42).unwrap();
        assert!((est - 0.5).abs() < 0.01, "got {}", est);
    }

    #[test]
    fn mc_1d_sin() {
        // ∫_0^π sin(x) dx = 2
        let (est, _se) = monte_carlo_integrate_1d(|x| x.sin(), 0.0, PI, 100000, 42).unwrap();
        assert!((est - 2.0).abs() < 0.02, "got {}", est);
    }

    #[test]
    fn mc_1d_reproducible() {
        // Same seed → same result
        let f = |x: f64| x * x;
        let (a, _) = monte_carlo_integrate_1d(f, 0.0, 1.0, 1000, 123).unwrap();
        let (b, _) = monte_carlo_integrate_1d(f, 0.0, 1.0, 1000, 123).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn mc_1d_invalid_args() {
        assert!(monte_carlo_integrate_1d(|x| x, 1.0, 0.0, 100, 42).is_err());
        assert!(monte_carlo_integrate_1d(|x| x, 0.0, 1.0, 0, 42).is_err());
    }

    #[test]
    fn mc_nd_2d_constant() {
        // ∫∫ 3 dA over [0,1]×[0,1] = 3
        let (est, se) =
            monte_carlo_integrate_nd(|_| 3.0, &[(0.0, 1.0), (0.0, 1.0)], 10000, 42).unwrap();
        assert!((est - 3.0).abs() < 1e-10, "got {}", est);
        assert!(se < 1e-10, "se should be ~0, got {}", se);
    }

    #[test]
    fn mc_nd_2d_xy() {
        // ∫∫ xy dA over [0,1]×[0,1] = 1/4
        let (est, _se) =
            monte_carlo_integrate_nd(|p| p[0] * p[1], &[(0.0, 1.0), (0.0, 1.0)], 200000, 42).unwrap();
        assert!((est - 0.25).abs() < 0.01, "got {}", est);
    }

    #[test]
    fn mc_nd_3d_unit_cube() {
        // ∫∫∫ 1 dV over [0,1]³ = 1
        let (est, _) = monte_carlo_integrate_nd(
            |_| 1.0,
            &[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)],
            50000,
            42,
        )
        .unwrap();
        assert!((est - 1.0).abs() < 0.01, "got {}", est);
    }

    #[test]
    fn mc_nd_invalid_args() {
        assert!(monte_carlo_integrate_nd(|_| 1.0, &[], 100, 42).is_err());
        assert!(monte_carlo_integrate_nd(|_| 1.0, &[(0.0, 1.0)], 0, 42).is_err());
        assert!(monte_carlo_integrate_nd(|_| 1.0, &[(1.0, 0.0)], 100, 42).is_err());
    }
}