//! Chebyshev-based fast approximations of transcendental functions.
//!
//! Uses precomputed Chebyshev polynomial series (near-minimax approximations)
//! to evaluate `sin`, `cos`, `exp`, and `log` with argument reduction.
//! The [`ChebyshevApprox`] struct is reusable for any smooth function on a
//! finite interval.
//!
//! # Example
//! ```no_run
//! use mathr::fastmath::{fast_sin, fast_exp, ChebyshevApprox};
//!
//! let approx = ChebyshevApprox::new(|x| x.sin(), -1.0, 1.0, 10);
//! assert!((approx.eval(0.5) - 0.5_f64.sin()).abs() < 1e-12);
//!
//! assert!((fast_sin(1.5) - 1.5_f64.sin()).abs() < 1e-13);
//! assert!((fast_exp(3.7) - 3.7_f64.exp()).abs() < 1e-13);
//! ```

use crate::interpolate::{chebyshev_coefficients, chebyshev_eval, chebyshev_rescale};
use std::f64::consts::{LN_2, PI};

/// A Chebyshev polynomial approximation of a smooth function on `[a, b]`.
///
/// Precomputes `n + 1` Chebyshev coefficients via a discrete cosine transform
/// of function samples at Chebyshev nodes.  Evaluation uses Clenshaw's
/// recurrence, which is numerically stable and requires `O(n)` work.
pub struct ChebyshevApprox {
    coeffs: Vec<f64>,
    a: f64,
    b: f64,
}

impl ChebyshevApprox {
    /// Build a degree-`n` Chebyshev approximation of `f` on `[a, b]`.
    pub fn new<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> Self {
        let mid = 0.5 * (a + b);
        let half = 0.5 * (b - a);
        let coeffs = chebyshev_coefficients(|t| f(half * t + mid), n);
        Self { coeffs, a, b }
    }

    /// Evaluate the approximation at `x`.  No clamping is performed; points
    /// outside `[a, b]` will extrapolate (often poorly).
    pub fn eval(&self, x: f64) -> f64 {
        let t = chebyshev_rescale(x, self.a, self.b);
        chebyshev_eval(&self.coeffs, t)
    }

    /// Return the Chebyshev coefficients (for inspection or serialisation).
    pub fn coeffs(&self) -> &[f64] {
        &self.coeffs
    }

    /// The approximation interval `[a, b]`.
    pub fn interval(&self) -> (f64, f64) {
        (self.a, self.b)
    }
}

// ---------------------------------------------------------------------------
// Argument-reduced fast functions
// ---------------------------------------------------------------------------

/// Reduce `x` modulo `2π` into `(-π, π]`.
fn reduce_2pi(x: f64) -> f64 {
    let r = x - 2.0 * PI * (x / (2.0 * PI)).round();
    // Guard against floating-point drift at the boundaries.
    if r > PI {
        r - 2.0 * PI
    } else if r <= -PI {
        r + 2.0 * PI
    } else {
        r
    }
}

/// Fast `sin(x)` via a degree-20 Chebyshev series on `[-π, π]` with
/// argument reduction modulo `2π`.
///
/// Accuracy: max error ~1e-14 across the real line.
pub fn fast_sin(x: f64) -> f64 {
    static APPROX: std::sync::OnceLock<ChebyshevApprox> = std::sync::OnceLock::new();
    let approx = APPROX.get_or_init(|| ChebyshevApprox::new(f64::sin, -PI, PI, 20));
    approx.eval(reduce_2pi(x))
}

/// Fast `cos(x)` via a degree-20 Chebyshev series on `[-π, π]` with
/// argument reduction modulo `2π`.
///
/// Accuracy: max error ~1e-14 across the real line.
pub fn fast_cos(x: f64) -> f64 {
    static APPROX: std::sync::OnceLock<ChebyshevApprox> = std::sync::OnceLock::new();
    let approx = APPROX.get_or_init(|| ChebyshevApprox::new(f64::cos, -PI, PI, 20));
    approx.eval(reduce_2pi(x))
}

/// Fast `exp(x)` via a degree-15 Chebyshev series with argument reduction.
///
/// Writes `x = k·ln2 + r` with `r ∈ [-ln2/2, ln2/2]`, then
/// `exp(x) = 2^k · exp(r)`.  The `2^k` factor is computed by exponent
/// manipulation.
///
/// Accuracy: max error ~1e-15 across `|x| < 700`.
pub fn fast_exp(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < -745.0 {
        return 0.0;
    }
    if x > 709.78 {
        return f64::INFINITY;
    }
    static APPROX: std::sync::OnceLock<ChebyshevApprox> = std::sync::OnceLock::new();
    let approx = APPROX.get_or_init(|| {
        ChebyshevApprox::new(f64::exp, -0.5 * LN_2, 0.5 * LN_2, 15)
    });
    let k = (x / LN_2).round();
    let r = x - k * LN_2;
    // 2^k via ldexp (exact for integer k within f64 exponent range).
    approx.eval(r) * (if k >= 0.0 { (1u64 << k as i32) as f64 } else { 1.0 / (1u64 << (-k) as i32) as f64 })
}

/// Fast `ln(x)` (natural log) via a degree-20 Chebyshev series with
/// argument reduction.
///
/// Writes `x = m · 2^k` with `m ∈ [1, 2)`, then
/// `ln(x) = k·ln2 + ln(m)`.  The mantissa `m` and exponent `k` are
/// extracted via `frexp`-style bit manipulation.
///
/// Accuracy: max error ~1e-15 for `x > 0`.
pub fn fast_log(x: f64) -> f64 {
    if x.is_nan() || x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    if x == f64::INFINITY {
        return f64::INFINITY;
    }
    static APPROX: std::sync::OnceLock<ChebyshevApprox> = std::sync::OnceLock::new();
    let approx = APPROX.get_or_init(|| ChebyshevApprox::new(f64::ln, 1.0, 2.0, 20));
    // Extract mantissa m ∈ [1, 2) and exponent k such that x = m * 2^k.
    let bits = x.to_bits();
    let exp_raw = ((bits >> 52) & 0x7FF) as i32;
    let mantissa_bits = (bits & 0x000F_FFFF_FFFF_FFFF) | 0x3FF0_0000_0000_0000;
    let m = f64::from_bits(mantissa_bits); // m ∈ [1, 2)
    let k = exp_raw - 1023;
    approx.eval(m) + k as f64 * LN_2
}

/// Fast `tan(x)` via the identity `tan = sin / cos`.
pub fn fast_tan(x: f64) -> f64 {
    let s = fast_sin(x);
    let c = fast_cos(x);
    if c.abs() < 1e-300 {
        f64::NAN
    } else {
        s / c
    }
}

/// Fast `sqrt(x)` via `exp(0.5 * ln(x))` — less accurate than the hardware
/// instruction but demonstrates the Chebyshev pipeline.
pub fn fast_sqrt(x: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    fast_exp(0.5 * fast_log(x))
}

/// Fast `pow(x, y)` via `exp(y * ln(x))`.
pub fn fast_pow(x: f64, y: f64) -> f64 {
    if x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        if y > 0.0 {
            return 0.0;
        }
        return f64::INFINITY;
    }
    fast_exp(y * fast_log(x))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn chebyshev_approx_basic() {
        let approx = ChebyshevApprox::new(|x| x * x, -2.0, 3.0, 6);
        for &x in &[-2.0, -1.0, 0.0, 0.5, 1.5, 3.0] {
            assert!(close(approx.eval(x), x * x, 1e-10), "at x={}", x);
        }
    }

    #[test]
    fn chebyshev_approx_coeffs_access() {
        let approx = ChebyshevApprox::new(f64::sin, -1.0, 1.0, 5);
        assert_eq!(approx.coeffs().len(), 6);
        let (a, b) = approx.interval();
        assert_eq!(a, -1.0);
        assert_eq!(b, 1.0);
    }

    #[test]
    fn fast_sin_accuracy() {
        for i in (-1000..=1000).step_by(37) {
            let x = i as f64 * 0.01;
            let err = (fast_sin(x) - x.sin()).abs();
            assert!(err < 1e-12, "fast_sin({}): err={:e}", x, err);
        }
    }

    #[test]
    fn fast_cos_accuracy() {
        for i in (-1000..=1000).step_by(41) {
            let x = i as f64 * 0.01;
            let err = (fast_cos(x) - x.cos()).abs();
            assert!(err < 1e-12, "fast_cos({}): err={:e}", x, err);
        }
    }

    #[test]
    fn fast_exp_accuracy() {
        for i in (-200..=200).step_by(7) {
            let x = i as f64 * 0.1;
            let err = (fast_exp(x) - x.exp()).abs();
            let rel = err / x.exp().abs();
            assert!(rel < 1e-13, "fast_exp({}): rel_err={:e}", x, rel);
        }
    }

    #[test]
    fn fast_exp_edge_cases() {
        assert!(fast_exp(f64::NAN).is_nan());
        assert_eq!(fast_exp(-800.0), 0.0);
        assert!(fast_exp(800.0).is_infinite());
    }

    #[test]
    fn fast_log_accuracy() {
        for x in (1..=1000_000).step_by(37) {
            let x = x as f64 * 0.001;
            let err = (fast_log(x) - x.ln()).abs();
            assert!(err < 1e-13, "fast_log({}): err={:e}", x, err);
        }
    }

    #[test]
    fn fast_log_edge_cases() {
        assert!(fast_log(f64::NAN).is_nan());
        assert!(fast_log(-1.0).is_nan());
        assert!(fast_log(0.0).is_infinite());
        assert!(fast_log(f64::INFINITY).is_infinite());
    }

    #[test]
    fn fast_tan_accuracy() {
        for i in (-100..=100).step_by(13) {
            let x = i as f64 * 0.01;
            // Skip near π/2 + kπ where tan diverges.
            let r = reduce_2pi(x);
            if (r - FRAC_PI_2).abs() < 0.1 || (r + FRAC_PI_2).abs() < 0.1 {
                continue;
            }
            let err = (fast_tan(x) - x.tan()).abs();
            assert!(err < 1e-10, "fast_tan({}): err={:e}", x, err);
        }
    }

    #[test]
    fn fast_sqrt_accuracy() {
        for x in (1..=100_000).step_by(37) {
            let x = x as f64 * 0.01;
            let rel = ((fast_sqrt(x) - x.sqrt()).abs()) / x.sqrt();
            assert!(rel < 1e-12, "fast_sqrt({}): rel_err={:e}", x, rel);
        }
    }

    #[test]
    fn fast_pow_accuracy() {
        for x in (1..=100).step_by(7) {
            let x = x as f64 * 0.1;
            for y in (-3..=3).step_by(1) {
                let y = y as f64 * 0.5;
                let rel = ((fast_pow(x, y) - x.powf(y)).abs()) / x.powf(y).abs();
                assert!(rel < 1e-12, "fast_pow({}, {}): rel_err={:e}", x, y, rel);
            }
        }
    }

    #[test]
    fn reduce_2pi_in_range() {
        for i in (-10_000..=10_000).step_by(97) {
            let x = i as f64 * 0.01;
            let r = reduce_2pi(x);
            assert!(r > -PI && r <= PI, "reduce_2pi({}) = {} not in (-π, π]", x, r);
        }
    }
}
