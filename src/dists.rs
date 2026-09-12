//! Probability distributions and classical hypothesis tests.
//!
//! Complements [`crate::stats`] (descriptive statistics, normal/exponential
//! pdf/cdf) with the distributions needed for inference — Student's t,
//! chi-squared, F, binomial, Poisson, uniform — plus p-value oriented
//! hypothesis tests (t-tests, chi-squared goodness-of-fit, one-way ANOVA).
//!
//! The key special function added here is the **regularized incomplete beta**
//! [`beta_inc`], which underlies the CDFs of the t, F, and binomial
//! distributions. The chi-squared and Poisson CDFs reuse
//! [`crate::special::incomplete_gamma_p`].
//!
//! All functions return plain `f64` (matching [`crate::special`] conventions);
//! parameter validation happens in the hypothesis-test helpers which return
//! `Result`.

use crate::error::Result;
use crate::special::{beta, incomplete_gamma_p, log_gamma};

/// Maximum iterations for the continued fraction in [`beta_inc`].
const MAX_IT: usize = 300;
/// Convergence threshold for the continued fraction.
const EPS: f64 = 3.0e-16;
/// Minimum positive value to avoid division by zero in the continued fraction.
const FPMIN: f64 = 1.0e-300;

/// Regularized incomplete beta function `I_x(a, b) = B(a, b; x) / B(a, b)`.
///
/// This is the CDF of the Beta(a, b) distribution and the building block for
/// Student-t, F, and binomial tail probabilities. Uses the continued fraction
/// from Numerical Recipes with the standard symmetry switch
/// (`x < (a+1)/(a+b+2)` evaluates `I_x(a,b)` directly, otherwise via
/// `1 - I_{1-x}(b,a)`).
pub fn beta_inc(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let ln_front = log_gamma(a + b) - log_gamma(a) - log_gamma(b)
        + a * x.ln()
        + b * (-x).ln_1p();
    let front = ln_front.exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * betacf(a, b, x) / a
    } else {
        1.0 - front * betacf(b, a, 1.0 - x) / b
    }
}

/// Continued fraction for the incomplete beta function (Lentz's method).
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAX_IT {
        let m = m as f64;
        let m2 = 2.0 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            return h;
        }
    }
    h
}

/// Quantile (inverse CDF) of the standard normal distribution.
///
/// Acklam's rational approximation with one Halley refinement step using
/// `erfc`, giving ~1e-15 relative accuracy over the full open interval
/// (0, 1). Returns NaN for p outside (0, 1).
#[allow(clippy::excessive_precision)] // canonical published Acklam coefficients
pub fn normal_ppf(p: f64) -> f64 {
    if !(p > 0.0 && p < 1.0) {
        return f64::NAN;
    }
    const A: [f64; 6] = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383577518672690e+02,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    const B: [f64; 5] = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    const C: [f64; 6] = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    const D: [f64; 4] = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];
    let p_low = 0.02425;
    let q = (-2.0 * p.ln()).sqrt();
    let x = if p < p_low {
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - p_low {
        let p_ = p - 0.5;
        let q_ = p_ * p_;
        (((((A[0] * q_ + A[1]) * q_ + A[2]) * q_ + A[3]) * q_ + A[4]) * q_ + A[5]) * p_
            / (((((B[0] * q_ + B[1]) * q_ + B[2]) * q_ + B[3]) * q_ + B[4]) * q_ + 1.0)
    } else {
        let q2 = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q2 + C[1]) * q2 + C[2]) * q2 + C[3]) * q2 + C[4]) * q2 + C[5])
            / ((((D[0] * q2 + D[1]) * q2 + D[2]) * q2 + D[3]) * q2 + 1.0)
    };
    // One Halley refinement step.
    let e = 0.5 * crate::special::erfc(-x / std::f64::consts::SQRT_2) - p;
    let u = e * (2.0 * std::f64::consts::PI).sqrt() * (x * x / 2.0).exp();
    x - u / (1.0 + x * u / 2.0)
}

/// Generic CDF inversion by bracketed bisection to floating-point precision.
///
/// The initial bracket `[lo, hi]` is expanded geometrically (doubling) until
/// it contains the root, then bisected until the interval can no longer be
/// split (i.e. `lo` and `hi` are adjacent doubles). Convergence is guaranteed
/// for a monotone CDF; each of the quantile functions below halves the
/// bracket on the correct side of the CDF value.
fn invert_cdf<F: Fn(f64) -> f64>(cdf: F, p: f64, lo: f64, hi: f64) -> f64 {
    let (mut lo, mut hi) = (lo, hi);
    let mut it = 0;
    while cdf(lo) > p && it < 200 {
        lo *= 2.0;
        it += 1;
    }
    it = 0;
    while cdf(hi) < p && it < 200 {
        hi *= 2.0;
        it += 1;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break; // adjacent f64s: machine precision reached
        }
        if cdf(mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// Quantile (inverse CDF) of Student's t distribution with `nu` degrees of
/// freedom. Returns NaN for p outside the open interval (0, 1) or nu <= 0.
pub fn student_t_ppf(p: f64, nu: f64) -> f64 {
    if !(p > 0.0 && p < 1.0) || nu <= 0.0 {
        return f64::NAN;
    }
    if p == 0.5 {
        return 0.0;
    }
    invert_cdf(|t| student_t_cdf(t, nu), p, -1.0, 1.0)
}

/// Quantile (inverse CDF) of the chi-squared distribution with `k` degrees
/// of freedom. Returns NaN for p outside the open interval (0, 1) or k <= 0.
pub fn chi2_ppf(p: f64, k: f64) -> f64 {
    if !(p > 0.0 && p < 1.0) || k <= 0.0 {
        return f64::NAN;
    }
    invert_cdf(|x| chi2_cdf(x, k), p, 0.0, 1.0)
}

/// Quantile (inverse CDF) of the F distribution with `(d1, d2)` degrees of
/// freedom. Returns NaN for p outside the open interval (0, 1) or
/// non-positive degrees of freedom.
pub fn f_ppf(p: f64, d1: f64, d2: f64) -> f64 {
    if !(p > 0.0 && p < 1.0) || d1 <= 0.0 || d2 <= 0.0 {
        return f64::NAN;
    }
    invert_cdf(|x| f_cdf(x, d1, d2), p, 0.0, 1.0)
}

/// Quantile (inverse CDF) of the continuous uniform distribution on `[a, b]`.
pub fn uniform_ppf(p: f64, a: f64, b: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) || b <= a {
        return f64::NAN;
    }
    a + p * (b - a)
}

/// PDF of Student's t distribution with `nu` degrees of freedom.
pub fn student_t_pdf(t: f64, nu: f64) -> f64 {
    if nu <= 0.0 {
        return f64::NAN;
    }
    let coef = log_gamma((nu + 1.0) / 2.0) - log_gamma(nu / 2.0)
        - 0.5 * (nu * std::f64::consts::PI).ln();
    (coef - (nu + 1.0) / 2.0 * (1.0 + t * t / nu).ln()).exp()
}

/// CDF of Student's t distribution with `nu` degrees of freedom.
pub fn student_t_cdf(t: f64, nu: f64) -> f64 {
    if nu <= 0.0 {
        return f64::NAN;
    }
    let x = nu / (nu + t * t);
    let tail = 0.5 * beta_inc(nu / 2.0, 0.5, x);
    if t > 0.0 {
        1.0 - tail
    } else if t < 0.0 {
        tail
    } else {
        0.5
    }
}

/// PDF of the chi-squared distribution with `k` degrees of freedom.
pub fn chi2_pdf(x: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    let half = k / 2.0;
    (-0.5 * k * std::f64::consts::LN_2 + (half - 1.0) * x.ln() - x / 2.0
        - log_gamma(half))
        .exp()
}

/// CDF of the chi-squared distribution with `k` degrees of freedom.
pub fn chi2_cdf(x: f64, k: f64) -> f64 {
    if k <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    incomplete_gamma_p(k / 2.0, x / 2.0)
}

/// PDF of the F distribution with `(d1, d2)` degrees of freedom.
pub fn f_pdf(x: f64, d1: f64, d2: f64) -> f64 {
    if d1 <= 0.0 || d2 <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    let num = (d1 / 2.0 * (d1 * x / d2).ln() - (d1 + d2) / 2.0 * (1.0 + d1 * x / d2).ln()).exp();
    num / (x * beta(d1 / 2.0, d2 / 2.0))
}

/// CDF of the F distribution with `(d1, d2)` degrees of freedom.
pub fn f_cdf(x: f64, d1: f64, d2: f64) -> f64 {
    if d1 <= 0.0 || d2 <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    beta_inc(d1 / 2.0, d2 / 2.0, d1 * x / (d1 * x + d2))
}

/// PMF of the binomial distribution `P(X = k)` for `k` successes in `n`
/// trials with success probability `p`.
pub fn binomial_pmf(k: u64, n: u64, p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) || k > n {
        return f64::NAN;
    }
    let k = k as f64;
    let n = n as f64;
    if p == 0.0 {
        return if k == 0.0 { 1.0 } else { 0.0 };
    }
    if p == 1.0 {
        return if k == n { 1.0 } else { 0.0 };
    }
    (log_gamma(n + 1.0) - log_gamma(k + 1.0) - log_gamma(n - k + 1.0)
        + k * p.ln()
        + (n - k) * (-p).ln_1p())
    .exp()
}

/// CDF of the binomial distribution `P(X <= k)`, via the identity
/// `P(X <= k) = I_{1-p}(n-k, k+1)`.
pub fn binomial_cdf(k: u64, n: u64, p: f64) -> f64 {
    if !(0.0..=1.0).contains(&p) || k > n {
        return f64::NAN;
    }
    if k == n {
        return 1.0;
    }
    if p == 0.0 {
        return 1.0;
    }
    beta_inc(n as f64 - k as f64, k as f64 + 1.0, 1.0 - p)
}

/// PMF of the Poisson distribution `P(X = k)` with rate `lambda`.
pub fn poisson_pmf(k: u64, lambda: f64) -> f64 {
    if lambda < 0.0 {
        return f64::NAN;
    }
    if lambda == 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    let k = k as f64;
    (-lambda + k * lambda.ln() - log_gamma(k + 1.0)).exp()
}

/// CDF of the Poisson distribution `P(X <= k)` via the upper incomplete
/// gamma: `P(X <= k) = 1 - P(k+1, lambda)`.
pub fn poisson_cdf(k: u64, lambda: f64) -> f64 {
    if lambda < 0.0 {
        return f64::NAN;
    }
    if lambda == 0.0 {
        return 1.0;
    }
    1.0 - incomplete_gamma_p(k as f64 + 1.0, lambda)
}

/// PDF of the continuous uniform distribution on `[a, b]`.
pub fn uniform_pdf(x: f64, a: f64, b: f64) -> f64 {
    if b <= a {
        return f64::NAN;
    }
    if x >= a && x <= b {
        1.0 / (b - a)
    } else {
        0.0
    }
}

/// CDF of the continuous uniform distribution on `[a, b]`.
pub fn uniform_cdf(x: f64, a: f64, b: f64) -> f64 {
    if b <= a {
        return f64::NAN;
    }
    if x < a {
        0.0
    } else if x > b {
        1.0
    } else {
        (x - a) / (b - a)
    }
}

/// Outcome of a hypothesis test: statistic, degrees of freedom, and p-value.
#[derive(Debug, Clone, PartialEq)]
pub struct TestResult {
    /// Human-readable test name, e.g. `"welch t-test"`.
    pub name: &'static str,
    /// Test statistic.
    pub statistic: f64,
    /// First degrees-of-freedom parameter (none for rank-based tests).
    pub df: Option<f64>,
    /// Second degrees-of-freedom parameter (F test only).
    pub df2: Option<f64>,
    /// Two-sided p-value for t-tests; upper-tail for chi-squared and F.
    pub p_value: f64,
}

/// Two-sided one-sample t-test of `H0: mean(sample) == mu`.
pub fn t_test_one(sample: &[f64], mu: f64) -> Result<TestResult> {
    if sample.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "t-test needs at least 2 samples".into(),
        ));
    }
    let n = sample.len() as f64;
    let m = crate::stats::mean(sample)?;
    let s = crate::stats::stddev_sample(sample)?;
    if s == 0.0 {
        return Err(crate::error::MathError::Eval(
            "sample has zero variance".into(),
        ));
    }
    let t = (m - mu) / (s / n.sqrt());
    let df = n - 1.0;
    let p = 2.0 * (1.0 - student_t_cdf(t.abs(), df));
    Ok(TestResult { name: "one-sample t-test", statistic: t, df: Some(df), df2: None, p_value: p })
}

/// Two-sided Welch's t-test comparing the means of two independent samples
/// (does not assume equal variances).
pub fn t_test_two(a: &[f64], b: &[f64]) -> Result<TestResult> {
    if a.len() < 2 || b.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "two-sample t-test needs at least 2 samples per group".into(),
        ));
    }
    let (n1, n2) = (a.len() as f64, b.len() as f64);
    let (m1, m2) = (crate::stats::mean(a)?, crate::stats::mean(b)?);
    let (v1, v2) = (
        crate::stats::variance_sample(a)?,
        crate::stats::variance_sample(b)?,
    );
    let se2 = v1 / n1 + v2 / n2;
    if se2 == 0.0 {
        return Err(crate::error::MathError::Eval(
            "both samples have zero variance".into(),
        ));
    }
    let t = (m1 - m2) / se2.sqrt();
    let df = se2 * se2
        / ((v1 / n1) * (v1 / n1) / (n1 - 1.0) + (v2 / n2) * (v2 / n2) / (n2 - 1.0));
    let p = 2.0 * (1.0 - student_t_cdf(t.abs(), df));
    Ok(TestResult { name: "welch t-test", statistic: t, df: Some(df), df2: None, p_value: p })
}

/// Two-sided paired t-test comparing the means of paired observations
/// (equivalent to a one-sample t-test on the differences against 0).
pub fn t_test_paired(a: &[f64], b: &[f64]) -> Result<TestResult> {
    if a.len() != b.len() {
        return Err(crate::error::MathError::Eval(
            "paired t-test needs equal-length samples".into(),
        ));
    }
    let diffs: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x - y).collect();
    let mut r = t_test_one(&diffs, 0.0)?;
    r.name = "paired t-test";
    Ok(r)
}

/// Pearson's chi-squared goodness-of-fit test: `observed[i]` vs `expected[i]`
/// counts. Degrees of freedom = `k - 1`.
pub fn chi_square_gof(observed: &[f64], expected: &[f64]) -> Result<TestResult> {
    if observed.len() < 2 || observed.len() != expected.len() {
        return Err(crate::error::MathError::Eval(
            "chi-squared GOF needs equal-length observed/expected counts (>= 2)".into(),
        ));
    }
    if expected.iter().any(|&e| e <= 0.0) {
        return Err(crate::error::MathError::Eval(
            "expected counts must be positive".into(),
        ));
    }
    let chi2: f64 = observed
        .iter()
        .zip(expected.iter())
        .map(|(&o, &e)| (o - e) * (o - e) / e)
        .sum();
    let df = observed.len() as f64 - 1.0;
    let p = 1.0 - chi2_cdf(chi2, df);
    Ok(TestResult { name: "chi-squared GOF", statistic: chi2, df: Some(df), df2: None, p_value: p })
}

/// Chi-squared test of equal cell probabilities: `observed` counts are
/// compared against the uniform expectation `total / k`. Degrees of freedom
/// = `k - 1`.
pub fn chi_square_uniform(observed: &[f64]) -> Result<TestResult> {
    if observed.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "chi-squared test needs at least 2 cells".into(),
        ));
    }
    let total: f64 = observed.iter().sum();
    let expected = vec![total / observed.len() as f64; observed.len()];
    let mut r = chi_square_gof(observed, &expected)?;
    r.name = "chi-squared uniform";
    Ok(r)
}

/// One-way analysis of variance (F-test) comparing the means of `groups`
/// independent samples. Degrees of freedom = `(k - 1, N - k)`.
pub fn anova_oneway(groups: &[&[f64]]) -> Result<TestResult> {
    if groups.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "ANOVA needs at least 2 groups".into(),
        ));
    }
    if groups.iter().any(|g| g.len() < 2) {
        return Err(crate::error::MathError::Eval(
            "each ANOVA group needs at least 2 samples".into(),
        ));
    }
    let k = groups.len() as f64;
    let n_total: f64 = groups.iter().map(|g| g.len() as f64).sum();
    let grand =
        crate::stats::mean(&groups.iter().flat_map(|g| g.iter()).copied().collect::<Vec<f64>>())?;
    let ss_between: Result<f64> = groups
        .iter()
        .map(|g| {
            let m = crate::stats::mean(g)?;
            Ok((m - grand).powi(2) * g.len() as f64)
        })
        .sum();
    let ss_between = ss_between?;
    let ss_within: f64 = groups
        .iter()
        .map(|g| {
            let m = crate::stats::mean(g)?;
            Ok(g.iter().map(|&x| (x - m) * (x - m)).sum::<f64>())
        })
        .sum::<Result<f64>>()?;
    let df1 = k - 1.0;
    let df2 = n_total - k;
    let (statistic, p) = if ss_within == 0.0 {
        if ss_between == 0.0 {
            (0.0, 1.0)
        } else {
            (f64::INFINITY, 0.0)
        }
    } else {
        let f = (ss_between / df1) / (ss_within / df2);
        (f, 1.0 - f_cdf(f, df1, df2))
    };
    Ok(TestResult {
        name: "one-way ANOVA",
        statistic,
        df: Some(df1),
        df2: Some(df2),
        p_value: p,
    })
}

/// Average ranks (1-based); tied values receive their mean rank.
fn average_ranks(values: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..values.len()).collect();
    idx.sort_by(|&a, &b| {
        values[a]
            .partial_cmp(&values[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut ranks = vec![0.0; values.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && values[idx[j + 1]] == values[idx[i]] {
            j += 1;
        }
        let avg = (i + j + 2) as f64 / 2.0;
        for &k in &idx[i..=j] {
            ranks[k] = avg;
        }
        i = j + 1;
    }
    ranks
}

/// Sum of `(t³ − t)` over tie groups (t = group size), for variance
/// corrections in rank-based tests. Zero when there are no ties.
fn tie_correction(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut total = 0.0;
    let mut i = 0;
    while i < sorted.len() {
        let mut j = i;
        while j + 1 < sorted.len() && sorted[j + 1] == sorted[i] {
            j += 1;
        }
        let t = (j - i + 1) as f64;
        if t > 1.0 {
            total += t * t * t - t;
        }
        i = j + 1;
    }
    total
}

/// Two-sided p-value from a standard-normal statistic, with the usual
/// continuity-corrected `z` clamped at zero (avoids p > 1 for near-exact
/// null results).
fn normal_two_sided_p(z: f64) -> f64 {
    (2.0 * (1.0 - crate::stats::normal_cdf(z, 0.0, 1.0))).min(1.0)
}

/// Two-sided Mann–Whitney U test (rank-sum): nonparametric comparison of
/// two independent samples, with continuity correction and tie-corrected
/// variance. The normal approximation is reliable for group sizes >= 8.
pub fn mann_whitney_u(a: &[f64], b: &[f64]) -> Result<TestResult> {
    if a.len() < 2 || b.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "Mann-Whitney U needs at least 2 samples per group".into(),
        ));
    }
    let n1 = a.len() as f64;
    let n2 = b.len() as f64;
    let n_total = n1 + n2;
    let combined: Vec<f64> = a.iter().chain(b.iter()).copied().collect();
    let ranks = average_ranks(&combined);
    let r1: f64 = ranks[..a.len()].iter().sum();
    let u1 = r1 - n1 * (n1 + 1.0) / 2.0;
    let mu = n1 * n2 / 2.0;
    let sigma = ((n1 * n2 / 12.0)
        * ((n_total + 1.0) - tie_correction(&combined) / (n_total * (n_total - 1.0))))
    .sqrt();
    let z = (((u1 - mu).abs() - 0.5).max(0.0)) / sigma;
    Ok(TestResult {
        name: "mann-whitney U",
        statistic: u1,
        df: None,
        df2: None,
        p_value: normal_two_sided_p(z),
    })
}

/// Two-sided Wilcoxon signed-rank test on paired samples: ranks the absolute
/// differences (zeros dropped) and tests the positive-rank sum `W+` against
/// its null distribution (normal approximation, continuity + tie corrected).
pub fn wilcoxon_signed_rank(a: &[f64], b: &[f64]) -> Result<TestResult> {
    if a.len() != b.len() || a.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "Wilcoxon signed-rank needs >= 2 equal-length pairs".into(),
        ));
    }
    let diffs: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x - y).collect();
    let abs: Vec<f64> = diffs.iter().map(|d| d.abs()).filter(|d| *d > 0.0).collect();
    let n = abs.len() as f64;
    if n < 2.0 {
        return Err(crate::error::MathError::Eval(
            "Wilcoxon signed-rank needs at least 2 non-zero differences".into(),
        ));
    }
    let ranks = average_ranks(&abs);
    let mut w_plus = 0.0;
    let mut k = 0;
    for &d in &diffs {
        if d != 0.0 {
            if d > 0.0 {
                w_plus += ranks[k];
            }
            k += 1;
        }
    }
    let mu = n * (n + 1.0) / 4.0;
    let sigma = (n * (n + 1.0) * (2.0 * n + 1.0) / 24.0 - tie_correction(&abs) / 48.0).sqrt();
    let z = (((w_plus - mu).abs() - 0.5).max(0.0)) / sigma;
    Ok(TestResult {
        name: "wilcoxon signed-rank",
        statistic: w_plus,
        df: None,
        df2: None,
        p_value: normal_two_sided_p(z),
    })
}

/// Kruskal–Wallis H test: nonparametric one-way ANOVA by ranks over `k`
/// groups, chi-squared with `k − 1` degrees of freedom (tie corrected).
/// The chi-squared approximation is reliable for group sizes >= 5.
pub fn kruskal_wallis(groups: &[&[f64]]) -> Result<TestResult> {
    if groups.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "Kruskal-Wallis needs at least 2 groups".into(),
        ));
    }
    if groups.iter().any(|g| g.is_empty()) {
        return Err(crate::error::MathError::Eval(
            "Kruskal-Wallis groups must be non-empty".into(),
        ));
    }
    let combined: Vec<f64> = groups.iter().flat_map(|g| g.iter()).copied().collect();
    let n_total = combined.len() as f64;
    let ranks = average_ranks(&combined);
    let mut h = 0.0;
    let mut start = 0usize;
    for g in groups {
        let r_j: f64 = ranks[start..start + g.len()].iter().sum();
        h += r_j * r_j / g.len() as f64;
        start += g.len();
    }
    h = 12.0 / (n_total * (n_total + 1.0)) * h - 3.0 * (n_total + 1.0);
    let ties = tie_correction(&combined);
    if ties > 0.0 {
        h /= 1.0 - ties / (n_total * n_total * n_total - n_total);
    }
    let df = groups.len() as f64 - 1.0;
    Ok(TestResult {
        name: "kruskal-wallis H",
        statistic: h,
        df: Some(df),
        df2: None,
        p_value: 1.0 - chi2_cdf(h, df),
    })
}

/// Spearman rank correlation: Pearson correlation of the average ranks of
/// `x` and `y`. Returns a value in `[-1, 1]` (NaN when a sample is constant).
pub fn spearman_corr(x: &[f64], y: &[f64]) -> Result<f64> {
    if x.len() != y.len() || x.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "Spearman correlation needs >= 2 equal-length samples".into(),
        ));
    }
    let rx = average_ranks(x);
    let ry = average_ranks(y);
    crate::stats::correlation(&rx, &ry)
}

/// Percentile bootstrap confidence interval for an arbitrary statistic.
///
/// Resamples `data` with replacement `iters` times (LCG-seeded [`crate::stats::Rng`],
/// reproducible for a fixed `seed`), evaluates `statistic` on each resample,
/// and returns the `(conf)` percentile interval of the bootstrap draws.
pub fn bootstrap_ci<F>(
    data: &[f64],
    statistic: F,
    iters: usize,
    conf: f64,
    seed: u64,
) -> Result<(f64, f64)>
where
    F: Fn(&[f64]) -> f64,
{
    if data.is_empty() {
        return Err(crate::error::MathError::Eval(
            "bootstrap needs at least one data point".into(),
        ));
    }
    if iters == 0 {
        return Err(crate::error::MathError::Eval(
            "bootstrap needs at least one iteration".into(),
        ));
    }
    if !(0.0..=1.0).contains(&conf) {
        return Err(crate::error::MathError::Eval(
            "confidence level must be in [0, 1]".into(),
        ));
    }
    let n = data.len();
    let mut rng = crate::stats::Rng::new(seed);
    let mut draws = Vec::with_capacity(iters);
    let mut sample = vec![0.0; n];
    for _ in 0..iters {
        for slot in sample.iter_mut() {
            let idx = (rng.uniform(0.0, n as f64) as usize).min(n - 1);
            *slot = data[idx];
        }
        draws.push(statistic(&sample));
    }
    draws.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let alpha = (1.0 - conf) / 2.0;
    let pick = |q: f64| draws[(q * (iters - 1) as f64).round() as usize];
    Ok((pick(alpha), pick(1.0 - alpha)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs().max(a.abs()))
    }

    // ---------- beta_inc ----------

    #[test]
    fn beta_inc_edges() {
        assert_eq!(beta_inc(2.0, 3.0, 0.0), 0.0);
        assert_eq!(beta_inc(2.0, 3.0, 1.0), 1.0);
        assert!(close(beta_inc(0.5, 0.5, 0.5), 0.5, 1e-14));
    }

    #[test]
    fn beta_inc_uniform_is_x() {
        // I_x(1, 1) = x (Beta(1,1) is Uniform(0,1))
        for &x in &[0.1, 0.25, 0.5, 0.75, 0.9] {
            assert!(close(beta_inc(1.0, 1.0, x), x, 1e-14), "x={x}");
        }
    }

    #[test]
    fn beta_inc_symmetry() {
        // I_x(a,b) + I_{1-x}(b,a) = 1
        for &(a, b, x) in &[(2.0, 3.0, 0.3), (0.5, 4.0, 0.7), (5.0, 5.0, 0.5)] {
            assert!(
                close(beta_inc(a, b, x) + beta_inc(b, a, 1.0 - x), 1.0, 1e-13),
                "a={a} b={b} x={x}"
            );
        }
    }

    #[test]
    fn beta_inc_matches_quadrature() {
        // I_x(2,3) = 12 ∫0^x t(1-t)^2 dt = 6x^2 - 8x^3 + 3x^4
        let x = 0.4f64;
        let exact = 6.0 * x * x - 8.0 * x * x * x + 3.0 * x.powi(4);
        assert!(close(beta_inc(2.0, 3.0, x), exact, 1e-13));
    }

    // ---------- normal_ppf ----------

    #[test]
    fn normal_ppf_known_quantiles() {
        assert!(close(normal_ppf(0.975), 1.959963984540054, 1e-10));
        assert!(close(normal_ppf(0.5), 0.0, 1e-12));
        assert!(close(normal_ppf(0.8413447460685429), 1.0, 1e-10));
    }

    #[test]
    fn normal_ppf_round_trip() {
        for &p in &[0.001, 0.025, 0.3, 0.7, 0.99, 0.999] {
            let x = normal_ppf(p);
            let cdf = crate::stats::normal_cdf(x, 0.0, 1.0);
            assert!(close(cdf, p, 1e-10), "p={p}");
        }
    }

    #[test]
    fn normal_ppf_symmetry() {
        assert!(close(normal_ppf(0.3), -normal_ppf(0.7), 1e-12));
        assert!(normal_ppf(0.0).is_nan());
        assert!(normal_ppf(1.0).is_nan());
    }

    // ---------- Student t ----------

    #[test]
    fn student_t_cdf_table_values() {
        // t_{0.975}(10) = 2.228 (standard t table)
        assert!(close(student_t_cdf(2.228, 10.0), 0.975, 1e-3));
        // t_{0.95}(5) = 2.015
        assert!(close(student_t_cdf(2.015, 5.0), 0.95, 1e-3));
        // nu = 1 is the Cauchy distribution: F(t) = 1/2 + atan(t)/pi
        for &t in &[-2.0f64, -0.5, 0.5, 2.0] {
            let cauchy = 0.5 + t.atan() / std::f64::consts::PI;
            assert!(close(student_t_cdf(t, 1.0), cauchy, 1e-12));
        }
        // median 0.5, symmetry
        assert!(close(student_t_cdf(0.0, 7.0), 0.5, 1e-14));
        assert!(close(student_t_cdf(-1.5, 7.0), 1.0 - student_t_cdf(1.5, 7.0), 1e-13));
    }

    #[test]
    fn student_t_pdf_integrates_to_one() {
        let area = crate::calculus::integrate_adaptive(
            |t| student_t_pdf(t, 6.0),
            -30.0,
            30.0,
            1e-10,
            60,
        )
        .unwrap();
        assert!(close(area, 1.0, 1e-5), "area={area}");
    }

    #[test]
    fn student_t_pdf_symmetric() {
        assert!(close(student_t_pdf(1.3, 8.0), student_t_pdf(-1.3, 8.0), 1e-14));
        // nu = 1 is Cauchy: pdf(0) = 1/pi
        assert!(close(student_t_pdf(0.0, 1.0), std::f64::consts::FRAC_1_PI, 1e-13));
        // nu = 2: pdf(0) = Gamma(1.5)/(sqrt(2 pi) Gamma(1)) = 1/(2 sqrt(2))
        assert!(close(student_t_pdf(0.0, 2.0), std::f64::consts::FRAC_1_SQRT_2 / 2.0, 1e-13));
    }

    // ---------- chi-squared ----------

    #[test]
    fn chi2_cdf_table_values() {
        // chi2_{0.95}(1) = 3.8415
        assert!(close(chi2_cdf(3.8415, 1.0), 0.95, 1e-4));
        // chi2_{0.95}(2) = 5.9915
        assert!(close(chi2_cdf(5.9915, 2.0), 0.95, 1e-4));
        // k=2 closed form: F(x) = 1 - e^{-x/2}
        assert!(close(chi2_cdf(2.0, 2.0), 1.0 - (-1.0f64).exp(), 1e-13));
    }

    #[test]
    fn chi2_pdf_integrates_to_one() {
        let area = crate::calculus::integrate_adaptive(
            |x| chi2_pdf(x, 4.0),
            0.0,
            80.0,
            1e-10,
            60,
        )
        .unwrap();
        assert!(close(area, 1.0, 1e-8));
    }

    // ---------- F distribution ----------

    #[test]
    fn f_cdf_consistent_with_t() {
        // F(1, nu) CDF at t^2 equals P(|T(nu)| <= t)
        let (t, nu) = (2.5, 12.0);
        let p_f = f_cdf(t * t, 1.0, nu);
        let p_t = 2.0 * student_t_cdf(t, nu) - 1.0;
        assert!(close(p_f, p_t, 1e-12));
    }

    #[test]
    fn f_pdf_integrates_to_one() {
        let area = crate::calculus::integrate_adaptive(
            |x| f_pdf(x, 3.0, 10.0),
            0.0,
            60.0,
            1e-10,
            60,
        )
        .unwrap();
        assert!(close(area, 1.0, 1e-5), "area={area}");
    }

    // ---------- binomial ----------

    #[test]
    fn binomial_pmf_fair_coin() {
        // B(10, 0.5): pmf(5) = 252/1024
        assert!(close(binomial_pmf(5, 10, 0.5), 252.0 / 1024.0, 1e-14));
        // pmf sums to 1
        let total: f64 = (0..=10).map(|k| binomial_pmf(k, 10, 0.5)).sum();
        assert!(close(total, 1.0, 1e-12));
        // degenerate cases
        assert_eq!(binomial_pmf(0, 5, 0.0), 1.0);
        assert_eq!(binomial_pmf(3, 5, 0.0), 0.0);
        assert_eq!(binomial_pmf(5, 5, 1.0), 1.0);
    }

    #[test]
    fn binomial_cdf_consistency() {
        // cdf(5) = sum pmf(0..=5) for B(10, 0.3)
        let sum: f64 = (0..=5).map(|k| binomial_pmf(k, 10, 0.3)).sum();
        assert!(close(binomial_cdf(5, 10, 0.3), sum, 1e-12));
        assert_eq!(binomial_cdf(10, 10, 0.3), 1.0);
    }

    // ---------- poisson ----------

    #[test]
    fn poisson_pmf_values() {
        // P(X=2 | lambda=2) = 2 e^{-2}
        assert!(close(poisson_pmf(2, 2.0), 2.0 * (-2.0f64).exp(), 1e-14));
        // P(X=0 | lambda=5) = e^{-5}
        assert!(close(poisson_pmf(0, 5.0), (-5.0f64).exp(), 1e-14));
    }

    #[test]
    fn poisson_cdf_values() {
        // P(X <= 2 | lambda=2) = 5 e^{-2}
        assert!(close(poisson_cdf(2, 2.0), 5.0 * (-2.0f64).exp(), 1e-13));
        // pmf sums to cdf
        let lambda = 3.7;
        let sum: f64 = (0..=6).map(|k| poisson_pmf(k, lambda)).sum();
        assert!(close(poisson_cdf(6, lambda), sum, 1e-12));
        assert_eq!(poisson_cdf(10, 0.0), 1.0);
    }

    // ---------- uniform ----------

    #[test]
    fn uniform_pdf_cdf() {
        assert!(close(uniform_pdf(1.5, 1.0, 3.0), 0.5, 1e-15));
        assert_eq!(uniform_pdf(0.5, 1.0, 3.0), 0.0);
        assert_eq!(uniform_pdf(4.0, 1.0, 3.0), 0.0);
        assert_eq!(uniform_cdf(1.0, 1.0, 3.0), 0.0);
        assert!(close(uniform_cdf(2.0, 1.0, 3.0), 0.5, 1e-15));
        assert_eq!(uniform_cdf(5.0, 1.0, 3.0), 1.0);
    }

    // ---------- quantiles ----------

    #[test]
    fn student_t_ppf_table_and_closed_forms() {
        // standard t table: t_{0.975}(10) = 2.228138852
        assert!(close(student_t_ppf(0.975, 10.0), 2.2281388519649385, 1e-9));
        assert!(close(student_t_ppf(0.95, 5.0), 2.015048372669157, 1e-9));
        // nu = 1 is Cauchy: ppf(p) = tan(pi (p - 1/2))
        for &p in &[0.1, 0.25, 0.75, 0.9] {
            let cauchy = (std::f64::consts::PI * (p - 0.5)).tan();
            assert!(close(student_t_ppf(p, 1.0), cauchy, 1e-9), "p={p}");
        }
        // median and symmetry
        assert_eq!(student_t_ppf(0.5, 7.0), 0.0);
        assert!(close(student_t_ppf(0.3, 7.0), -student_t_ppf(0.7, 7.0), 1e-12));
        assert!(student_t_ppf(0.0, 7.0).is_nan());
        assert!(student_t_ppf(1.0, 7.0).is_nan());
    }

    #[test]
    fn chi2_ppf_closed_forms() {
        // k = 1: chi2_p = (normal_ppf_{(1+p)/2})^2
        let expect1 = normal_ppf(0.975).powi(2);
        assert!(close(chi2_ppf(0.95, 1.0), expect1, 1e-10));
        // k = 2: chi2_p = -2 ln(1 - p)
        assert!(close(chi2_ppf(0.95, 2.0), -2.0 * (0.05f64).ln(), 1e-10));
        // standard table: chi2_{0.95}(5) = 11.070497693516351
        assert!(close(chi2_ppf(0.95, 5.0), 11.070497693516351, 1e-9));
        assert!(chi2_ppf(0.0, 3.0).is_nan());
    }

    #[test]
    fn f_ppf_consistency() {
        // F_{0.95}(1, nu) = t_{0.975}(nu)^2
        let t = student_t_ppf(0.975, 10.0);
        assert!(close(f_ppf(0.95, 1.0, 10.0), t * t, 1e-8));
        // round-trip against the CDF
        for &p in &[0.1, 0.5, 0.9, 0.99] {
            let x = f_ppf(p, 3.0, 10.0);
            assert!(close(f_cdf(x, 3.0, 10.0), p, 1e-12), "p={p}");
        }
        assert!(f_ppf(0.0, 3.0, 10.0).is_nan());
    }

    #[test]
    fn ppf_cdf_round_trips() {
        for &p in &[0.001, 0.025, 0.25, 0.75, 0.975, 0.999] {
            let t = student_t_ppf(p, 12.0);
            assert!(close(student_t_cdf(t, 12.0), p, 1e-12), "t p={p}");
            let x = chi2_ppf(p, 7.0);
            assert!(close(chi2_cdf(x, 7.0), p, 1e-12), "chi2 p={p}");
        }
    }

    #[test]
    fn uniform_ppf_basic() {
        assert!(close(uniform_ppf(0.25, 1.0, 3.0), 1.5, 1e-15));
        assert_eq!(uniform_ppf(0.0, 1.0, 3.0), 1.0);
        assert_eq!(uniform_ppf(1.0, 1.0, 3.0), 3.0);
        assert!(uniform_ppf(-0.1, 1.0, 3.0).is_nan());
    }

    // ---------- hypothesis tests ----------

    #[test]
    fn t_test_one_known_statistic() {
        // data 1..=5, mu=0: mean=3, s^2=2.5, t = 3/(sqrt(2.5/5)) = 3*sqrt(2)
        let r = t_test_one(&[1.0, 2.0, 3.0, 4.0, 5.0], 0.0).unwrap();
        assert!(close(r.statistic, 3.0 * 2.0f64.sqrt(), 1e-12));
        assert_eq!(r.df, Some(4.0));
        // two-sided p: |t| = 4.243 is between t_0.99(4) = 3.747 and t_0.995(4) = 4.604
        // -> p in (0.01, 0.02)
        assert!(r.p_value > 0.01 && r.p_value < 0.02);
        // mu equal to the sample mean -> t = 0, p = 1
        let r0 = t_test_one(&[1.0, 2.0, 3.0, 4.0, 5.0], 3.0).unwrap();
        assert!(close(r0.statistic, 0.0, 1e-14));
        assert!(close(r0.p_value, 1.0, 1e-12));
    }

    #[test]
    fn t_test_two_equal_variances_match_pooled_direction() {
        // Two symmetric groups: [1,2,3] vs [4,5,6] -> t negative, |t| large
        let r = t_test_two(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap();
        // s1^2 = s2^2 = 1, se = sqrt(1/3+1/3), t = -3/sqrt(2/3)
        assert!(close(r.statistic, -3.0 / (2.0f64 / 3.0).sqrt(), 1e-12));
        // df = 4 for equal variances and equal sizes
        assert_eq!(r.df, Some(4.0));
        // |t| = 3.674 < t_0.99(4) = 3.747 -> two-sided p > 0.02, clearly significant
        assert!(r.p_value > 0.02 && r.p_value < 0.05);
        // identical means -> t = 0, p = 1
        let r0 = t_test_two(&[1.0, 2.0, 3.0], &[3.0, 2.0, 1.0]).unwrap();
        assert!(close(r0.statistic, 0.0, 1e-14));
        assert!(close(r0.p_value, 1.0, 1e-12));
    }

    #[test]
    fn t_test_paired_is_one_sample_on_diffs() {
        let a = [8.0, 7.0, 6.0, 9.0, 10.0];
        let b = [6.0, 7.0, 5.0, 7.0, 8.0];
        let paired = t_test_paired(&a, &b).unwrap();
        let diffs: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x - y).collect();
        let direct = t_test_one(&diffs, 0.0).unwrap();
        assert!(close(paired.statistic, direct.statistic, 1e-14));
        assert!(close(paired.p_value, direct.p_value, 1e-14));
        assert_eq!(paired.name, "paired t-test");
    }

    #[test]
    fn t_test_validation_errors() {
        assert!(t_test_one(&[1.0], 0.0).is_err());
        assert!(t_test_two(&[1.0], &[1.0, 2.0]).is_err());
        assert!(t_test_paired(&[1.0, 2.0], &[1.0]).is_err());
        // zero variance rejected
        assert!(t_test_one(&[2.0, 2.0, 2.0], 0.0).is_err());
    }

    #[test]
    fn chi_square_uniform_known() {
        // 6 cells, total 88, uniform: chi2 = 2.0 exactly
        let obs = [16.0, 18.0, 16.0, 14.0, 12.0, 12.0];
        let r = chi_square_uniform(&obs).unwrap();
        assert!(close(r.statistic, 2.0, 1e-12));
        assert_eq!(r.df, Some(5.0));
        // 2.0 << 11.07 critical value at 0.05 -> large p
        assert!(r.p_value > 0.5);
        // perfect fit -> chi2 = 0, p = 1
        let r0 = chi_square_uniform(&[10.0, 10.0, 10.0]).unwrap();
        assert!(close(r0.statistic, 0.0, 1e-15));
        assert!(close(r0.p_value, 1.0, 1e-12));
    }

    #[test]
    fn chi_square_gof_two_sided_counts() {
        // Fair die, 60 rolls, expected 10 each:
        // obs = [6, 14, 10, 10, 8, 12] -> chi2 = (16+16+0+0+4+4)/10 = 4.0, df=5
        let obs = [6.0, 14.0, 10.0, 10.0, 8.0, 12.0];
        let exp = vec![10.0; 6];
        let r = chi_square_gof(&obs, &exp).unwrap();
        assert!(close(r.statistic, 4.0, 1e-12));
        assert!(close(r.p_value, 1.0 - chi2_cdf(4.0, 5.0), 1e-12));
        assert!(chi_square_gof(&[1.0, 2.0], &[0.0, 1.0]).is_err());
        assert!(chi_square_gof(&[1.0], &[1.0]).is_err());
    }

    #[test]
    fn anova_known_statistic() {
        // groups [1,2,3], [4,5,6], [7,8,9]: SSB = 3*(9+0+9) = 54, SSW = 6
        // F = (54/2)/(6/6) = 27, df = (2, 6)
        let g1 = [1.0, 2.0, 3.0];
        let g2 = [4.0, 5.0, 6.0];
        let g3 = [7.0, 8.0, 9.0];
        let r = anova_oneway(&[&g1, &g2, &g3]).unwrap();
        assert!(close(r.statistic, 27.0, 1e-12));
        assert_eq!(r.df, Some(2.0));
        assert_eq!(r.df2, Some(6.0));
        assert!(r.p_value < 0.01);
    }

    #[test]
    fn anova_identical_groups() {
        let g1 = [1.0, 2.0, 3.0];
        let g2 = [1.0, 2.0, 3.0];
        let r = anova_oneway(&[&g1, &g2]).unwrap();
        assert!(close(r.statistic, 0.0, 1e-15));
        assert!(close(r.p_value, 1.0, 1e-12));
    }

    #[test]
    fn anova_validation_errors() {
        let g1 = [1.0, 2.0];
        assert!(anova_oneway(&[&g1]).is_err());
        let short = [1.0];
        assert!(anova_oneway(&[&g1, &short]).is_err());
    }

    // ---------- nonparametric tests ----------

    #[test]
    fn mann_whitney_perfect_separation() {
        // disjoint groups of 8 -> U1 = 0 exactly, tiny p (normal approximation)
        let a: Vec<f64> = (1..=8).map(|i| i as f64).collect();
        let b: Vec<f64> = (9..=16).map(|i| i as f64).collect();
        let r = mann_whitney_u(&a, &b).unwrap();
        assert!(close(r.statistic, 0.0, 1e-14));
        assert!(r.df.is_none());
        assert!(r.p_value < 0.001, "p={}", r.p_value);
        // identical distributions -> U = mu, p = 1 (z clamped at 0)
        let r0 = mann_whitney_u(&[1.0, 2.0, 3.0, 4.0], &[1.0, 2.0, 3.0, 4.0]).unwrap();
        assert!(close(r0.statistic, 8.0, 1e-14));
        assert!(close(r0.p_value, 1.0, 1e-12));
    }

    #[test]
    fn mann_whitney_hand_computed_with_ties() {
        // a = [1, 2, 5], b = [2, 3, 4]:
        // pooled sorted: 1 2 2 3 4 5 -> ranks: 1, 2.5, 2.5, 4, 5, 6
        // R1 = 1 + 2.5 + 6 = 9.5 -> U1 = 9.5 - 6 = 3.5
        // N = 6, mu = 4.5, ties (t=2): sigma^2 = (9/12)[7 - 6/30] = 0.75 * 6.8
        let r = mann_whitney_u(&[1.0, 2.0, 5.0], &[2.0, 3.0, 4.0]).unwrap();
        assert!(close(r.statistic, 3.5, 1e-14));
        let sigma = (0.75f64 * (7.0 - 6.0 / 30.0)).sqrt();
        let z = ((3.5f64 - 4.5).abs() - 0.5) / sigma;
        let expect = 2.0 * (1.0 - crate::stats::normal_cdf(z, 0.0, 1.0));
        assert!(close(r.p_value, expect, 1e-14));
        assert!(r.p_value > 0.7 && r.p_value < 0.95, "p={}", r.p_value);
    }

    #[test]
    fn mann_whitney_shifted_groups() {
        // clearly shifted groups, no ties: p small
        let a: Vec<f64> = (1..=10).map(|i| i as f64).collect();
        let b: Vec<f64> = (11..=20).map(|i| i as f64).collect();
        let r = mann_whitney_u(&a, &b).unwrap();
        assert!(close(r.statistic, 0.0, 1e-14));
        assert!(r.p_value < 0.001, "p={}", r.p_value);
    }

    #[test]
    fn mann_whitney_symmetric_p() {
        // swapping groups mirrors U around mu = n1 n2 / 2, same two-sided p
        let a = [1.0, 3.0, 4.0, 7.0];
        let b = [2.0, 5.0, 6.0, 8.0];
        let r1 = mann_whitney_u(&a, &b).unwrap();
        let r2 = mann_whitney_u(&b, &a).unwrap();
        assert!(close(r1.statistic + r2.statistic, 16.0, 1e-12));
        assert!(close(r1.p_value, r2.p_value, 1e-14));
    }

    #[test]
    fn wilcoxon_hand_computed() {
        // diffs a-b = [2,0,1,2,2] -> drop zero; |d| = [2,1,2,2]
        // ranks: 1 -> 1, 2 -> (2+3+4)/3 = 3; all diffs positive
        // W+ = 1 + 3 + 3 + 3 = 10; n=4: mu = 5,
        // tie group is three 2s: t^3-t = 24 -> sigma^2 = 7.5 - 24/48 = 7.0
        let a = [8.0, 7.0, 6.0, 9.0, 10.0];
        let b = [6.0, 7.0, 5.0, 7.0, 8.0];
        let r = wilcoxon_signed_rank(&a, &b).unwrap();
        assert!(close(r.statistic, 10.0, 1e-14));
        let sigma = (7.0f64).sqrt();
        let z = ((10.0f64 - 5.0).abs() - 0.5) / sigma;
        let expect = 2.0 * (1.0 - crate::stats::normal_cdf(z, 0.0, 1.0));
        assert!(close(r.p_value, expect, 1e-14));
        assert!(r.p_value > 0.05 && r.p_value < 0.2, "p={}", r.p_value);
    }

    #[test]
    fn wilcoxon_all_zero_differences_error() {
        assert!(wilcoxon_signed_rank(&[1.0, 2.0], &[1.0, 2.0]).is_err());
        assert!(wilcoxon_signed_rank(&[1.0], &[2.0]).is_err());
        assert!(wilcoxon_signed_rank(&[1.0, 2.0], &[1.0, 2.0, 3.0]).is_err());
    }

    #[test]
    fn kruskal_wallis_hand_computed() {
        // textbook triple: R = [6, 15, 24] -> H = (12/90)(279) - 30 = 7.2, df = 2
        let g1 = [2.0, 4.0, 3.0];
        let g2 = [5.0, 6.0, 7.0];
        let g3 = [8.0, 10.0, 9.0];
        let r = kruskal_wallis(&[&g1, &g2, &g3]).unwrap();
        assert!(close(r.statistic, 7.2, 1e-12));
        assert_eq!(r.df, Some(2.0));
        assert!(close(r.p_value, 1.0 - chi2_cdf(7.2, 2.0), 1e-14));
        assert!(r.p_value < 0.05, "p={}", r.p_value);
        // identical groups -> H = 0, p = 1
        let r0 = kruskal_wallis(&[&g1, &g1]).unwrap();
        assert!(close(r0.statistic, 0.0, 1e-14));
        assert!(close(r0.p_value, 1.0, 1e-12));
    }

    #[test]
    fn kruskal_wallis_validation() {
        let g1 = [1.0, 2.0];
        let empty: [f64; 0] = [];
        assert!(kruskal_wallis(&[&g1]).is_err());
        assert!(kruskal_wallis(&[&g1, &empty]).is_err());
    }

    #[test]
    fn spearman_correlation() {
        // monotone increasing -> 1, decreasing -> -1
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(close(spearman_corr(&x, &x).unwrap(), 1.0, 1e-14));
        let rev = [5.0, 4.0, 3.0, 2.0, 1.0];
        assert!(close(spearman_corr(&x, &rev).unwrap(), -1.0, 1e-14));
        // handmade with rank inversion: x=[1,2,3,4], y=[1,3,2,5] -> rho = 0.8
        let x4 = [1.0, 2.0, 3.0, 4.0];
        let y = [1.0, 3.0, 2.0, 5.0];
        assert!(close(spearman_corr(&x4, &y).unwrap(), 0.8, 1e-12));
        // ties handled via average ranks: y ties -> |rho| < 1
        let y2 = [1.0, 2.0, 2.0, 4.0];
        let r = spearman_corr(&x4, &y2).unwrap();
        assert!(r > 0.85 && r < 1.0, "rho={r}");
        assert!(spearman_corr(&[1.0], &[2.0]).is_err());
    }

    // ---------- bootstrap ----------

    #[test]
    fn bootstrap_ci_mean_deterministic() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let (lo, hi) =
            bootstrap_ci(&data, |s| crate::stats::mean(s).unwrap(), 2000, 0.95, 42).unwrap();
        assert!(lo <= 3.5 && hi >= 3.5, "ci=({lo}, {hi})");
        assert!(lo < hi);
        // reproducible: same seed -> identical interval
        let (lo2, hi2) =
            bootstrap_ci(&data, |s| crate::stats::mean(s).unwrap(), 2000, 0.95, 42).unwrap();
        assert_eq!((lo, hi), (lo2, hi2));
    }

    #[test]
    fn bootstrap_ci_median_and_conf() {
        let data: Vec<f64> = (1..=101).map(|i| i as f64).collect(); // median 51
        let (lo, hi) =
            bootstrap_ci(&data, |s| crate::stats::median(s).unwrap(), 1000, 0.90, 7).unwrap();
        assert!(lo <= 51.0 && hi >= 51.0, "ci=({lo}, {hi})");
        // 90% interval is contained in the 99% interval for the same seed
        let (lo99, hi99) =
            bootstrap_ci(&data, |s| crate::stats::median(s).unwrap(), 1000, 0.99, 7).unwrap();
        assert!(lo99 <= lo && hi <= hi99);
    }

    #[test]
    fn bootstrap_ci_validation() {
        assert!(bootstrap_ci(&[], |s| s.len() as f64, 10, 0.95, 1).is_err());
        assert!(bootstrap_ci(&[1.0], |s| s[0], 0, 0.95, 1).is_err());
        assert!(bootstrap_ci(&[1.0], |s| s[0], 10, 1.5, 1).is_err());
    }
}
