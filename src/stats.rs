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
}
