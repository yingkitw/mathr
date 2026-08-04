//! Descriptive statistics from scratch.
//!
//! Provides mean, median, variance, standard deviation, quartiles,
//! Pearson correlation, and linear regression.

use crate::error::{MathError, Result};

/// Arithmetic mean.
pub fn mean(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(MathError::InvalidArgument("mean: empty data".into()));
    }
    Ok(data.iter().sum::<f64>() / data.len() as f64)
}

/// Median (middle value for odd-length, average of two middle for even).
pub fn median(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(MathError::InvalidArgument("median: empty data".into()));
    }
    let mut sorted: Vec<f64> = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    if n % 2 == 1 {
        Ok(sorted[n / 2])
    } else {
        Ok((sorted[n / 2 - 1] + sorted[n / 2]) / 2.0)
    }
}

/// Population variance (dividing by N).
pub fn variance(data: &[f64]) -> Result<f64> {
    if data.is_empty() {
        return Err(MathError::InvalidArgument("variance: empty data".into()));
    }
    let m = mean(data)?;
    Ok(data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / data.len() as f64)
}

/// Sample variance (dividing by N-1, Bessel's correction).
pub fn variance_sample(data: &[f64]) -> Result<f64> {
    if data.len() < 2 {
        return Err(MathError::InvalidArgument("variance_sample: need at least 2 data points".into()));
    }
    let m = mean(data)?;
    Ok(data.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (data.len() - 1) as f64)
}

/// Population standard deviation.
pub fn stddev(data: &[f64]) -> Result<f64> {
    Ok(variance(data)?.sqrt())
}

/// Sample standard deviation.
pub fn stddev_sample(data: &[f64]) -> Result<f64> {
    Ok(variance_sample(data)?.sqrt())
}

/// Minimum value.
pub fn min(data: &[f64]) -> Result<f64> {
    data.iter()
        .fold(None, |acc, &x| match acc {
            None => Some(x),
            Some(m) if x < m => Some(x),
            other => other,
        })
        .ok_or_else(|| MathError::InvalidArgument("min: empty data".into()))
}

/// Maximum value.
pub fn max(data: &[f64]) -> Result<f64> {
    data.iter()
        .fold(None, |acc, &x| match acc {
            None => Some(x),
            Some(m) if x > m => Some(x),
            other => other,
        })
        .ok_or_else(|| MathError::InvalidArgument("max: empty data".into()))
}

/// Range (max - min).
pub fn range(data: &[f64]) -> Result<f64> {
    Ok(max(data)? - min(data)?)
}

/// Quartiles: returns (Q1, Q2/median, Q3) using linear interpolation.
pub fn quartiles(data: &[f64]) -> Result<(f64, f64, f64)> {
    if data.len() < 4 {
        return Err(MathError::InvalidArgument("quartiles: need at least 4 data points".into()));
    }
    let mut sorted: Vec<f64> = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();

    let percentile = |p: f64| -> f64 {
        let pos = p * (n - 1) as f64;
        let lo = pos.floor() as usize;
        let hi = pos.ceil() as usize;
        if lo == hi {
            sorted[lo]
        } else {
            let frac = pos - lo as f64;
            sorted[lo] * (1.0 - frac) + sorted[hi] * frac
        }
    };

    Ok((percentile(0.25), percentile(0.5), percentile(0.75)))
}

/// Interquartile range (Q3 - Q1).
pub fn iqr(data: &[f64]) -> Result<f64> {
    let (q1, _, q3) = quartiles(data)?;
    Ok(q3 - q1)
}

/// Pearson correlation coefficient between two datasets.
pub fn correlation(x: &[f64], y: &[f64]) -> Result<f64> {
    if x.len() != y.len() {
        return Err(MathError::InvalidArgument("correlation: length mismatch".into()));
    }
    if x.len() < 2 {
        return Err(MathError::InvalidArgument("correlation: need at least 2 points".into()));
    }
    let mx = mean(x)?;
    let my = mean(y)?;
    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        let dy = y[i] - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }
    let denom = (dx2 * dy2).sqrt();
    if denom < 1e-30 {
        return Err(MathError::InvalidArgument("correlation: zero variance".into()));
    }
    Ok(num / denom)
}

/// Simple linear regression: returns (slope, intercept) for y = slope*x + intercept.
pub fn linear_regression(x: &[f64], y: &[f64]) -> Result<(f64, f64)> {
    if x.len() != y.len() {
        return Err(MathError::InvalidArgument("linear_regression: length mismatch".into()));
    }
    if x.len() < 2 {
        return Err(MathError::InvalidArgument("linear_regression: need at least 2 points".into()));
    }
    let mx = mean(x)?;
    let my = mean(y)?;
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..x.len() {
        let dx = x[i] - mx;
        num += dx * (y[i] - my);
        den += dx * dx;
    }
    if den < 1e-30 {
        return Err(MathError::InvalidArgument("linear_regression: zero variance in x".into()));
    }
    let slope = num / den;
    let intercept = my - slope * mx;
    Ok((slope, intercept))
}

/// Summary statistics in one call.
pub struct Summary {
    pub count: usize,
    pub mean: f64,
    pub median: f64,
    pub stddev: f64,
    pub min: f64,
    pub max: f64,
    pub range: f64,
}

pub fn summary(data: &[f64]) -> Result<Summary> {
    if data.is_empty() {
        return Err(MathError::InvalidArgument("summary: empty data".into()));
    }
    Ok(Summary {
        count: data.len(),
        mean: mean(data)?,
        median: median(data)?,
        stddev: stddev(data)?,
        min: min(data)?,
        max: max(data)?,
        range: range(data)?,
    })
}

// ---------------------------------------------------------------------------
// Stochastic primitives: reproducible RNG, sampling, distributions
// ---------------------------------------------------------------------------

/// Reproducible pseudo-random number generator (LCG with PCG-style constants).
pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng {
            state: seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407),
        }
    }

    /// Uniform sample in `[0, 1)`.
    pub fn next_f64(&mut self) -> f64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.state >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform sample in `[lo, hi)`.
    pub fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        lo + self.next_f64() * (hi - lo)
    }

    /// Standard normal sample via Box–Muller transform.
    pub fn standard_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-300);
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        r * theta.cos()
    }

    /// Normal sample with given `mean` and `sigma`.
    pub fn normal(&mut self, mean: f64, sigma: f64) -> f64 {
        mean + sigma * self.standard_normal()
    }

    /// Exponential sample with rate `lambda` via inverse-CDF method.
    pub fn exponential(&mut self, lambda: f64) -> f64 {
        let u = self.next_f64().max(1e-300);
        -u.ln() / lambda
    }
}

// --- Distribution PDFs / CDFs ---

/// Normal probability density function.
pub fn normal_pdf(x: f64, mean: f64, sigma: f64) -> f64 {
    let z = (x - mean) / sigma;
    (-0.5 * z * z).exp() / (sigma * (2.0 * std::f64::consts::PI).sqrt())
}

/// Normal cumulative distribution function via the error function.
pub fn normal_cdf(x: f64, mean: f64, sigma: f64) -> f64 {
    let z = (x - mean) / (sigma * std::f64::consts::SQRT_2);
    0.5 * (1.0 + crate::special::erf(z))
}

/// Exponential probability density function with rate `lambda`.
pub fn exp_pdf(x: f64, lambda: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else {
        lambda * (-lambda * x).exp()
    }
}

/// Exponential cumulative distribution function with rate `lambda`.
pub fn exp_cdf(x: f64, lambda: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else {
        1.0 - (-lambda * x).exp()
    }
}

// --- Moments and cumulants ---

/// Raw moments up to order 4: returns (mean, variance, skewness, excess kurtosis).
///
/// - **Skewness**: `E[(X-μ)³] / σ³` — measures asymmetry
/// - **Excess kurtosis**: `E[(X-μ)⁴] / σ⁴ - 3` — measures tail heaviness (0 for normal)
pub fn moments(data: &[f64]) -> Result<(f64, f64, f64, f64)> {
    let n = data.len();
    if n < 2 {
        return Err(MathError::InvalidArgument(
            "moments: need at least 2 data points".into(),
        ));
    }
    let m1 = mean(data)?;
    let mut m2 = 0.0;
    let mut m3 = 0.0;
    let mut m4 = 0.0;
    for &x in data {
        let d = x - m1;
        let d2 = d * d;
        m2 += d2;
        m3 += d2 * d;
        m4 += d2 * d2;
    }
    let nf = n as f64;
    let var = m2 / nf;
    let sigma = var.sqrt();
    if sigma < 1e-30 {
        return Err(MathError::InvalidArgument("moments: zero variance".into()));
    }
    let skew = (m3 / nf) / (sigma * sigma * sigma);
    let kurt = (m4 / nf) / (var * var) - 3.0;
    Ok((m1, var, skew, kurt))
}

/// Cumulants up to order 4 from central moments.
///
/// Returns (κ1=mean, κ2=variance, κ3, κ4) where:
/// - `κ3 = E[(X-μ)³]`
/// - `κ4 = E[(X-μ)⁴] - 3·Var²`
pub fn cumulants(data: &[f64]) -> Result<(f64, f64, f64, f64)> {
    let n = data.len();
    if n < 2 {
        return Err(MathError::InvalidArgument(
            "cumulants: need at least 2 data points".into(),
        ));
    }
    let m1 = mean(data)?;
    let mut m2 = 0.0;
    let mut m3 = 0.0;
    let mut m4 = 0.0;
    for &x in data {
        let d = x - m1;
        let d2 = d * d;
        m2 += d2;
        m3 += d2 * d;
        m4 += d2 * d2;
    }
    let nf = n as f64;
    let var = m2 / nf;
    let k3 = m3 / nf;
    let k4 = m4 / nf - 3.0 * var * var;
    Ok((m1, var, k3, k4))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn mean_basic() {
        assert!(close(mean(&[1.0, 2.0, 3.0, 4.0, 5.0]).unwrap(), 3.0));
    }

    #[test]
    fn median_odd() {
        assert!(close(median(&[3.0, 1.0, 2.0]).unwrap(), 2.0));
    }

    #[test]
    fn median_even() {
        assert!(close(median(&[1.0, 2.0, 3.0, 4.0]).unwrap(), 2.5));
    }

    #[test]
    fn variance_and_stddev() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let v = variance(&data).unwrap();
        assert!(close(v, 4.0));
        assert!(close(stddev(&data).unwrap(), 2.0));
    }

    #[test]
    fn sample_variance() {
        let data = vec![2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let v = variance_sample(&data).unwrap();
        // population var = 4, sample var = 4 * 8 / 7
        assert!(close(v, 32.0 / 7.0));
    }

    #[test]
    fn min_max_range() {
        let data = vec![3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0];
        assert!(close(min(&data).unwrap(), 1.0));
        assert!(close(max(&data).unwrap(), 9.0));
        assert!(close(range(&data).unwrap(), 8.0));
    }

    #[test]
    fn correlation_perfect() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0];
        assert!(close(correlation(&x, &y).unwrap(), 1.0));
    }

    #[test]
    fn correlation_negative() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 8.0, 6.0, 4.0, 2.0];
        assert!(close(correlation(&x, &y).unwrap(), -1.0));
    }

    #[test]
    fn linear_regression_basic() {
        let x = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let y = vec![1.0, 3.0, 5.0, 7.0, 9.0];
        let (slope, intercept) = linear_regression(&x, &y).unwrap();
        assert!(close(slope, 2.0));
        assert!(close(intercept, 1.0));
    }

    #[test]
    fn quartiles_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let (q1, q2, q3) = quartiles(&data).unwrap();
        assert!(close(q1, 2.75));
        assert!(close(q2, 4.5));
        assert!(close(q3, 6.25));
    }

    #[test]
    fn summary_works() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let s = summary(&data).unwrap();
        assert_eq!(s.count, 5);
        assert!(close(s.mean, 3.0));
        assert!(close(s.median, 3.0));
        assert!(close(s.min, 1.0));
        assert!(close(s.max, 5.0));
    }

    // --- Stochastic primitive tests ---

    #[test]
    fn rng_uniform_range() {
        let mut rng = Rng::new(42);
        for _ in 0..1000 {
            let x = rng.next_f64();
            assert!(x >= 0.0 && x < 1.0, "out of range: {}", x);
        }
    }

    #[test]
    fn rng_uniform_bounds() {
        let mut rng = Rng::new(42);
        for _ in 0..1000 {
            let x = rng.uniform(-5.0, 5.0);
            assert!(x >= -5.0 && x < 5.0, "out of range: {}", x);
        }
    }

    #[test]
    fn rng_reproducible() {
        let mut a = Rng::new(123);
        let mut b = Rng::new(123);
        for _ in 0..100 {
            assert_eq!(a.next_f64(), b.next_f64());
        }
    }

    #[test]
    fn rng_normal_mean() {
        let mut rng = Rng::new(42);
        let n = 100000;
        let samples: Vec<f64> = (0..n).map(|_| rng.standard_normal()).collect();
        let m = mean(&samples).unwrap();
        assert!(m.abs() < 0.02, "mean should be ~0, got {}", m);
        let v = variance(&samples).unwrap();
        assert!((v - 1.0).abs() < 0.03, "variance should be ~1, got {}", v);
    }

    #[test]
    fn rng_normal_shifted() {
        let mut rng = Rng::new(42);
        let n = 100000;
        let samples: Vec<f64> = (0..n).map(|_| rng.normal(5.0, 2.0)).collect();
        let m = mean(&samples).unwrap();
        assert!((m - 5.0).abs() < 0.03, "mean should be ~5, got {}", m);
        let v = variance(&samples).unwrap();
        assert!((v - 4.0).abs() < 0.1, "variance should be ~4, got {}", v);
    }

    #[test]
    fn rng_exponential_mean() {
        let mut rng = Rng::new(42);
        let n = 100000;
        let lambda = 2.0;
        let samples: Vec<f64> = (0..n).map(|_| rng.exponential(lambda)).collect();
        let m = mean(&samples).unwrap();
        // E[Exp(λ)] = 1/λ = 0.5
        assert!((m - 0.5).abs() < 0.01, "mean should be ~0.5, got {}", m);
    }

    #[test]
    fn normal_pdf_peak() {
        // PDF at mean = 1/(sigma * sqrt(2*pi))
        let p = normal_pdf(0.0, 0.0, 1.0);
        let expected = 1.0 / (2.0 * std::f64::consts::PI).sqrt();
        assert!(close(p, expected));
    }

    #[test]
    fn normal_pdf_symmetry() {
        let p1 = normal_pdf(1.0, 0.0, 1.0);
        let p2 = normal_pdf(-1.0, 0.0, 1.0);
        assert!(close(p1, p2));
    }

    #[test]
    fn normal_cdf_bounds() {
        assert!(close(normal_cdf(-100.0, 0.0, 1.0), 0.0));
        assert!(close(normal_cdf(100.0, 0.0, 1.0), 1.0));
        assert!(close(normal_cdf(0.0, 0.0, 1.0), 0.5));
    }

    #[test]
    fn exp_pdf_cdf() {
        let lambda = 2.0;
        assert!(close(exp_pdf(0.0, lambda), lambda));
        assert!(close(exp_pdf(-1.0, lambda), 0.0));
        assert!(close(exp_cdf(0.0, lambda), 0.0));
        assert!(close(exp_cdf(-1.0, lambda), 0.0));
        // CDF at 1/lambda should be 1 - e^{-1} ≈ 0.6321
        let c = exp_cdf(1.0 / lambda, lambda);
        assert!((c - (1.0 - std::f64::consts::E.powi(-1))).abs() < 1e-9);
    }

    #[test]
    fn moments_symmetric() {
        // Symmetric data → skewness ≈ 0
        let data = vec![-2.0, -1.0, 0.0, 1.0, 2.0];
        let (m, v, skew, kurt) = moments(&data).unwrap();
        assert!(close(m, 0.0));
        assert!(close(v, 2.0));
        assert!(skew.abs() < 1e-10, "skew should be 0, got {}", skew);
        // Excess kurtosis for uniform-like = -1.3
        assert!((kurt - (-1.3)).abs() < 0.01, "kurt should be ~-1.3, got {}", kurt);
    }

    #[test]
    fn moments_skewed() {
        // Right-skewed data → positive skewness
        let data = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100.0];
        let (_, _, skew, _) = moments(&data).unwrap();
        assert!(skew > 0.0, "skew should be positive, got {}", skew);
    }

    #[test]
    fn cumulants_basic() {
        // For normal-like data: κ4 ≈ 0
        let mut rng = Rng::new(42);
        let data: Vec<f64> = (0..50000).map(|_| rng.standard_normal()).collect();
        let (_, var, _, k4) = cumulants(&data).unwrap();
        assert!((var - 1.0).abs() < 0.05, "var should be ~1, got {}", var);
        assert!(k4.abs() < 0.1, "k4 should be ~0 for normal, got {}", k4);
    }

    #[test]
    fn moments_zero_variance() {
        let data = vec![5.0, 5.0, 5.0];
        assert!(moments(&data).is_err());
    }
}
