//! Special functions from scratch.
//!
//! Provides the Gamma function (Lanczos approximation), Beta function,
//! error function (erf) and complementary error function (erfc),
//! and the sinc function.

use std::f64::consts::PI;

/// Lanczos approximation coefficients (g=7, n=9).
const LANCZOS_G: f64 = 7.0;
const LANCZOS_C: [f64; 9] = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
];

/// Gamma function Γ(z) via the Lanczos approximation.
///
/// Accurate to ~15 significant digits for positive real z.
/// For negative non-integer z, uses the reflection formula.
pub fn gamma(z: f64) -> f64 {
    if z < 0.5 {
        // Reflection formula: Γ(z)Γ(1-z) = π / sin(πz)
        let sin_pi_z = (PI * z).sin();
        if sin_pi_z.abs() < 1e-15 {
            return f64::INFINITY; // pole at non-positive integer
        }
        PI / (sin_pi_z * gamma(1.0 - z))
    } else {
        let z = z - 1.0;
        let mut x = LANCZOS_C[0];
        for i in 1..LANCZOS_C.len() {
            x += LANCZOS_C[i] / (z + i as f64);
        }
        let t = z + LANCZOS_G + 0.5;
        (2.0 * PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * x
    }
}

/// Log-gamma function ln(Γ(z)) for positive z.
pub fn log_gamma(z: f64) -> f64 {
    if z <= 0.0 {
        return f64::NAN;
    }
    let g = gamma(z);
    if g.is_infinite() || g <= 0.0 {
        return f64::INFINITY;
    }
    g.ln()
}

/// Beta function B(a, b) = Γ(a)Γ(b) / Γ(a+b).
pub fn beta(a: f64, b: f64) -> f64 {
    gamma(a) * gamma(b) / gamma(a + b)
}

/// Error function erf(x).
///
/// Uses the identity erf(x) = sign(x) * P(0.5, x²), where P is the
/// regularized lower incomplete gamma function.
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x.is_infinite() {
        return if x > 0.0 { 1.0 } else { -1.0 };
    }
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    sign * incomplete_gamma_p(0.5, x * x)
}

/// Complementary error function erfc(x) = 1 - erf(x).
pub fn erfc(x: f64) -> f64 {
    1.0 - erf(x)
}

/// Normalized sinc function: sinc(x) = sin(πx) / (πx), with sinc(0) = 1.
pub fn sinc(x: f64) -> f64 {
    if x.abs() < 1e-15 {
        1.0
    } else {
        (PI * x).sin() / (PI * x)
    }
}

/// Unnormalized sinc: sin(x) / x, with sinc_u(0) = 1.
pub fn sinc_unnorm(x: f64) -> f64 {
    if x.abs() < 1e-15 {
        1.0
    } else {
        x.sin() / x
    }
}

/// Bessel function of the first kind, order 0: J₀(x).
///
/// Uses the Maclaurin series
/// `J_0(x) = Σ (-1)^k (x/2)^(2k) / (k! k!)`
/// for small `|x|` and the asymptotic expansion for larger `|x|`.
pub fn bessel_j0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        // J_0(x) = Σ (-1)^k (x²/4)^k / (k!)²
        let xsq = x * x / 4.0;
        let mut term = 1.0_f64;
        let mut sum = term;
        for k in 1..60 {
            let kf = k as f64;
            term *= -xsq / (kf * kf);
            sum += term;
            if term.abs() < 1e-17 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    } else {
        // Asymptotic: J_0(x) ~ sqrt(2/(πx)) * [P cos θ - Q sin θ] with θ = x - π/4.
        let z = 8.0 / ax;
        let y = z * z;
        let p = 1.0
            + y * (-0.1098628627e-2
                + y * (0.7464519654e-3
                    + y * (-0.4724987825e-4
                        + y * (0.2181196076e-5
                            + y * (-0.6397653302e-7 + y * 0.9538904063e-9)))));
        let q = -0.1562499995e-1
            + y * (0.1430484407e-3
                + y * (-0.4253339102e-4
                    + y * (0.2493458662e-5
                        + y * (-0.1248279047e-6 + y * 0.2860702546e-8))));
        let xx = ax - PI / 4.0;
        let result = (p * xx.cos() - z * q * xx.sin()) / ax.sqrt();
        // The x < 0 case: J_0 is even.
        result.abs() * result.signum()
    }
}

/// Bessel function of the first kind, order 1: J₁(x).
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    if ax < 8.0 {
        // J_1(x) = (x/2) Σ (-1)^k (x²/4)^k / (k!(k+1)!)
        let xsq = x * x / 4.0;
        let mut term = 0.5 * x; // k=0: (x/2)^1 / (0! 1!)
        let mut sum = term;
        for k in 1..60 {
            let kf = k as f64;
            term *= -xsq / (kf * (kf + 1.0));
            sum += term;
            if term.abs() < 1e-17 * sum.abs().max(1e-300) {
                break;
            }
        }
        sum
    } else {
        // Asymptotic.
        let z = 8.0 / ax;
        let y = z * z;
        let p = 1.0
            + y * (0.183105e-2
                + y * (-0.3516396496e-3
                    + y * (0.2457529642e-4
                        + y * (-0.2403370194e-5
                            + y * 0.1058465960e-7))));
        let q = 0.4687499995e-1
            + y * (-0.2002690873e-3
                + y * (0.4717512717e-4
                    + y * (-0.9414049147e-6
                        + y * (0.1344888788e-7 + y * -0.2199534093e-9))));
        let xx = ax - 3.0 * PI / 4.0;
        let result = (p * xx.cos() - z * q * xx.sin()) / ax.sqrt();
        sign * result
    }
}

/// Bessel function of the first kind, integer order n: Jₙ(x).
///
/// Computes J_0, J_1 via polynomial approximation, then uses the
/// forward recurrence `J_{n+1}(x) = (2n/x) J_n(x) - J_{n-1}(x)`
/// (or the series formula for very small `|x|` or for `n > |x|`,
/// where forward recurrence is unstable).
pub fn bessel_jn(n: i32, x: f64) -> f64 {
    if n == 0 {
        return bessel_j0(x);
    }
    if n == 1 {
        return bessel_j1(x);
    }
    if n < 0 {
        // J_{-n}(x) = (-1)^n J_n(x) for integer n.
        let jn = bessel_jn(-n, x);
        return if (-n) % 2 == 1 { -jn } else { jn };
    }
    if x.abs() < 1e-15 {
        return 0.0;
    }
    let n_u = n as u32;

    // Series form: J_n(x) = (x/2)^n * Σ (-1)^k (x/2)^{2k} / (k! (n+k)!).
    // Stable for n ≳ x (small x relative to n).
    let half_x = x.abs() / 2.0;
    if (n as f64) > x.abs() {
        let mut term = 1.0 / factorial_u64(n_u);
        let mut sum = term;
        let xx = half_x * half_x;
        for k in 1..200u32 {
            term *= -xx / (k as f64 * (n_u + k) as f64);
            let next = sum + term;
            if (next - sum).abs() < 1e-18 * sum.abs().max(1e-300) {
                break;
            }
            sum = next;
        }
        let mag = sum * half_x.powi(n);
        if x < 0.0 && n_u % 2 == 1 { -mag } else { mag }
    } else {
        // Forward recurrence: stable for n < x.
        let tox = 2.0 / x;
        let mut prev = bessel_j0(x);
        let mut curr = bessel_j1(x);
        for k in 1..n_u {
            let next = (k as f64) * tox * curr - prev;
            prev = curr;
            curr = next;
            if curr.abs() > 1e150 {
                return 0.0;
            }
        }
        if x < 0.0 && n_u % 2 == 1 { -curr } else { curr }
    }
}

fn factorial_u64(n: u32) -> f64 {
    let mut f = 1.0_f64;
    for i in 2..=n {
        f *= i as f64;
    }
    f
}

/// Incomplete gamma function P(a, x) = γ(a, x) / Γ(a)
/// via series expansion (for x < a+1) or continued fraction (for x >= a+1).
pub fn incomplete_gamma_p(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }

    let gln = log_gamma(a);

    if x < a + 1.0 {
        // Series: P(a, x) = (x^a * e^{-x} / Γ(a)) * Σ x^n / (a(a+1)...(a+n))
        let mut term = 1.0 / a;
        let mut sum = term;
        for n in 1..200 {
            term *= x / (a + n as f64);
            sum += term;
            if term.abs() < 1e-18 * sum.abs() {
                break;
            }
        }
        sum * x.powf(a) * (-x).exp() / gln.exp()
    } else {
        // Continued fraction for Q(a, x) = 1 - P(a, x)
        // Q(a, x) = (e^{-x} x^a / Γ(a)) * CF
        // CF = 1/(x+1-a - 1*(1-a)/(x+3-a - 2*(2-a)/(x+5-a - ...)))
        let tiny = 1e-30;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut f = d;
        for n in 1..300 {
            let an = -(n as f64) * (n as f64 - a);
            b = b + 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let delta = c * d;
            f *= delta;
            if (delta - 1.0).abs() < 1e-16 {
                break;
            }
        }
        let q = f * (-x + a * x.ln() - gln).exp();
        1.0 - q
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    fn rel_close(a: f64, b: f64, tol: f64) -> bool {
        if b.abs() < 1e-15 {
            (a - b).abs() < tol
        } else {
            ((a - b) / b).abs() < tol
        }
    }

    #[test]
    fn gamma_half() {
        // Γ(1/2) = √π
        assert!(close(gamma(0.5), PI.sqrt(), 1e-10));
    }

    #[test]
    fn gamma_integers() {
        // Γ(n) = (n-1)!
        assert!(close(gamma(1.0), 1.0, 1e-10));
        assert!(close(gamma(2.0), 1.0, 1e-10));
        assert!(close(gamma(3.0), 2.0, 1e-10));
        assert!(close(gamma(4.0), 6.0, 1e-10));
        assert!(close(gamma(5.0), 24.0, 1e-10));
        assert!(close(gamma(6.0), 120.0, 1e-10));
    }

    #[test]
    fn gamma_reflection() {
        // Γ(z)Γ(1-z) = π / sin(πz)
        let z = 0.3;
        let product = gamma(z) * gamma(1.0 - z);
        assert!(close(product, PI / (PI * z).sin(), 1e-9));
    }

    #[test]
    fn beta_function() {
        // B(1,1) = 1
        assert!(close(beta(1.0, 1.0), 1.0, 1e-10));
        // B(2,2) = 1/6
        assert!(close(beta(2.0, 2.0), 1.0 / 6.0, 1e-10));
        // B(0.5, 0.5) = π
        assert!(close(beta(0.5, 0.5), PI, 1e-9));
    }

    #[test]
    fn erf_basic() {
        assert!(close(erf(0.0), 0.0, 1e-15));
        // erf(0.5) ≈ 0.5204998778
        assert!(close(erf(0.5), 0.5204998778, 1e-8));
        // erf(1.0) ≈ 0.8427007929
        assert!(close(erf(1.0), 0.8427007929, 1e-8));
        // erf(2.0) ≈ 0.9953222650
        assert!(close(erf(2.0), 0.9953222650, 1e-8));
    }

    #[test]
    fn erf_negative() {
        // erf is odd: erf(-x) = -erf(x)
        assert!(close(erf(-1.0), -erf(1.0), 1e-14));
        assert!(close(erf(-0.5), -erf(0.5), 1e-14));
    }

    #[test]
    fn erfc_basic() {
        assert!(close(erfc(0.0), 1.0, 1e-15));
        // erfc(1.0) ≈ 0.1572992071
        assert!(close(erfc(1.0), 0.1572992071, 1e-8));
        // erfc(x) = 1 - erf(x)
        for &x in &[0.1, 0.5, 1.0, 2.0, 3.0] {
            assert!(close(erfc(x), 1.0 - erf(x), 1e-12));
        }
    }

    #[test]
    fn erf_large() {
        // erf(inf) = 1, erf(-inf) = -1
        assert!(close(erf(f64::INFINITY), 1.0, 1e-15));
        assert!(close(erf(f64::NEG_INFINITY), -1.0, 1e-15));
    }

    #[test]
    fn sinc_basic() {
        assert!(close(sinc(0.0), 1.0, 1e-15));
        // sinc(1) = sin(π)/π = 0
        assert!(close(sinc(1.0), 0.0, 1e-15));
        // sinc(0.5) = sin(π/2)/(π/2) = 2/π
        assert!(close(sinc(0.5), 2.0 / PI, 1e-14));
    }

    #[test]
    fn sinc_unnorm_basic() {
        assert!(close(sinc_unnorm(0.0), 1.0, 1e-15));
        // sin(x)/x at x=π = 0
        assert!(close(sinc_unnorm(PI), 0.0, 1e-14));
    }

    #[test]
    fn log_gamma_positive() {
        // ln(Γ(5)) = ln(24)
        assert!(close(log_gamma(5.0), 24.0f64.ln(), 1e-10));
    }

    #[test]
    fn incomplete_gamma_chi2() {
        // P(a, x) for chi-squared distribution check:
        // P(1, 0) = 0, P(1, ∞) = 1
        assert!(close(incomplete_gamma_p(1.0, 0.0), 0.0, 1e-15));
        // P(2, 2) ≈ 0.2642 (chi-squared with 2 df at x=2)
        let p = incomplete_gamma_p(1.0, 2.0);
        // P(1, x) = 1 - e^{-x}
        assert!(close(p, 1.0 - (-2.0f64).exp(), 1e-8));
    }

    #[test]
    fn bessel_j0_basic() {
        // J_0(0) = 1
        assert!(close(bessel_j0(0.0), 1.0, 1e-12));
        // J_0 has its first zero near x = 2.4048
        assert!(close(bessel_j0(2.4048255576957727), 0.0, 1e-6));
        // J_0(5) ≈ -0.17759677131434
        assert!(close(bessel_j0(5.0), -0.17759677131434, 1e-8));
    }

    #[test]
    fn bessel_j1_basic() {
        // J_1(0) = 0
        assert!(close(bessel_j1(0.0), 0.0, 1e-12));
        // J_1(2) ≈ 0.5767248077568736
        assert!(close(bessel_j1(2.0), 0.5767248077568736, 1e-8));
        // First zero of J_1 near x = 3.83171
        assert!(close(bessel_j1(3.8317059702075125), 0.0, 1e-6));
    }

    #[test]
    fn bessel_jn_positive() {
        // J_n matches J_0 and J_1 for n=0,1
        for &x in &[0.5, 1.0, 2.0, 5.0, 10.0] {
            assert!(close(bessel_jn(0, x), bessel_j0(x), 1e-10));
            assert!(close(bessel_jn(1, x), bessel_j1(x), 1e-10));
        }
        // J_2(5) ≈ 0.0465651
        assert!(close(bessel_jn(2, 5.0), 0.0465651, 1e-6));
        // J_3(2) ≈ 0.128943249
        assert!(close(bessel_jn(3, 2.0), 0.128943249, 1e-6));
        // J_10(5) — nonzero value, sanity test
        let v = bessel_jn(10, 5.0);
        assert!(v.is_finite());
    }

    #[test]
    fn bessel_jn_negative() {
        // J_{-n}(x) = (-1)^n J_n(x)
        for &x in &[1.0, 3.0, 5.0] {
            for n in [1, 2, 3, 4] {
                let jn = bessel_jn(n, x);
                let jneg = bessel_jn(-n, x);
                let expected = if n % 2 == 1 { -jn } else { jn };
                assert!(close(jneg, expected, 1e-10), "n={} x={}: {} vs {}", n, x, jneg, expected);
            }
        }
    }

    #[test]
    fn bessel_jn_at_zero() {
        // J_n(0) = 0 for n >= 1
        assert!(close(bessel_jn(5, 0.0), 0.0, 1e-12));
    }

    #[test]
    fn bessel_recurrence() {
        // Recurrence J_{n+1}(x) = (2n/x) J_n(x) - J_{n-1}(x)
        for x in [1.0, 3.0, 5.0, 8.0] {
            let jnm1 = bessel_jn(1, x);
            let jn = bessel_jn(2, x);
            let jnp1 = bessel_jn(3, x);
            let lhs = jnp1;
            let rhs = 2.0 * 2.0 / x * jn - jnm1;
            assert!(close(lhs, rhs, 1e-8), "x={}: {} vs {}", x, lhs, rhs);
        }
    }
}

// =========================================================================
// Digamma / Trigamma / Polygamma
// =========================================================================

/// Digamma function ψ(x) = d/dx ln Γ(x).
///
/// Reflection for x < 0.5, recurrence ψ(x+1) = ψ(x) + 1/x to push the
/// argument above 6, then the asymptotic series
/// `ψ(x) ≈ ln x − 1/(2x) − 1/(12x²) + 1/(120x⁴) − 1/(252x⁶) + …`.
/// Poles at non-positive integers return NaN.
pub fn digamma(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }
    let mut x = x;
    let mut result = 0.0f64;
    if x < 0.0 {
        // ψ(x) = ψ(1−x) − π/tan(πx)
        result -= PI / (PI * x).tan();
        x = 1.0 - x;
    }
    while x < 12.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    result += x.ln() - 0.5 * inv
        - inv2 * (1.0 / 12.0
            - inv2 * (1.0 / 120.0
                - inv2 * (1.0 / 252.0
                    - inv2 * (1.0 / 240.0 - inv2 * (1.0 / 132.0)))));
    result
}

/// Trigamma function ψ₁(x) = d²/dx² ln Γ(x).
///
/// Same strategy as [`digamma`]: reflection
/// `ψ₁(x) + ψ₁(1−x) = π²/sin²(πx)`, recurrence
/// `ψ₁(x+1) = ψ₁(x) − 1/x²`, asymptotic series for large arguments.
pub fn trigamma(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 && x == x.floor() {
        return f64::NAN;
    }
    if x < 0.0 {
        let sin_term = (PI * x).sin();
        return PI * PI / (sin_term * sin_term) - trigamma(1.0 - x);
    }
    let mut x = x;
    let mut result = 0.0f64;
    while x < 14.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    result += inv + 0.5 * inv2
        + inv2 * inv * (1.0 / 6.0
            - inv2 * (1.0 / 30.0 - inv2 * (1.0 / 42.0 - inv2 * (1.0 / 30.0))));
    result
}

/// Polygamma function of order m: ψ⁽ᵐ⁾(x), the m-th derivative of ψ(x).
///
/// `m = 0` → [`digamma`], `m = 1` → [`trigamma`]; for `m ≥ 2` uses
/// `ψ⁽ᵐ⁾(x) = (−1)^(m+1) · m! · ζ(m+1, x)` with the Hurwitz zeta.
pub fn polygamma(m: u32, x: f64) -> f64 {
    match m {
        0 => digamma(x),
        1 => trigamma(x),
        _ => {
            let sign: f64 = if (m + 1) % 2 == 0 { 1.0 } else { -1.0 };
            let mfact: f64 = (1..=m as u64).product::<u64>() as f64;
            sign * mfact * hurwitz_zeta(m as f64 + 1.0, x)
        }
    }
}

/// Harmonic number H_n = 1 + 1/2 + … + 1/n.
///
/// Direct summation for small `n`; `H_n = ψ(n+1) + γ` for large `n`.
/// `H_0 = 0`.
pub fn harmonic(n: u64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    if n <= 100 {
        return (1..=n).map(|k| 1.0 / k as f64).sum();
    }
    digamma(n as f64 + 1.0) + std::f64::consts::EULER_GAMMA
}

// =========================================================================
// Hurwitz / Riemann zeta
// =========================================================================

/// Hurwitz zeta ζ(s, a) = Σ_{k=0}^∞ 1/(k+a)^s for `a > 0`, `s ≠ 1`,
/// via Euler–Maclaurin summation (12 direct terms + 6 Bernoulli
/// correction terms). Valid for any real `s ≠ 1`; the pole at `s = 1`
/// returns ±inf.
pub fn hurwitz_zeta(s: f64, a: f64) -> f64 {
    if s.is_nan() || a.is_nan() || a <= 0.0 {
        return f64::NAN;
    }
    if s == 1.0 {
        return f64::INFINITY;
    }
    const N: u64 = 12;
    let mut sum: f64 = (0..N).map(|k| (a + k as f64).powf(-s)).sum();
    let t = a + N as f64;
    // Tail: t^{1−s}/(s−1) + ½·t^{−s} + Σ B₂ⱼ/(2ⱼ)! · s^{↑(2ⱼ−2)} · t^{−s−2ⱼ+1}
    sum += t.powf(1.0 - s) / (s - 1.0);
    sum += 0.5 * t.powf(-s);
    // B₂ⱼ/(2ⱼ)! for j = 1..=6
    const BERN: [f64; 6] = [
        1.0 / 12.0,
        -1.0 / 720.0,
        1.0 / 30240.0,
        -1.0 / 1209600.0,
        1.0 / 47900160.0,
        -691.0 / 1307674368000.0,
    ];
    // Term j needs rising = s^{↑(2j−1)} and power = t^{−(s+2j−1)}: the
    // first step is one power of t, subsequent steps add t^{−2}.
    let mut rising = s;
    let mut power = t.powf(-s - 1.0);
    for j in 1..=6u64 {
        sum += BERN[(j - 1) as usize] * rising * power;
        rising *= (s + (2 * j - 1) as f64) * (s + 2.0 * j as f64);
        power /= t * t;
    }
    sum
}

/// Riemann zeta function ζ(s).
///
/// Computed via [`hurwitz_zeta`] at `a = 1` (Euler–Maclaurin, valid for
/// any real `s ≠ 1`). ζ(0) = −½ and negative even integers (trivial
/// zeros) are special-cased to exact values; the pole at `s = 1` returns
/// +inf.
pub fn zeta(s: f64) -> f64 {
    if s.is_nan() {
        return f64::NAN;
    }
    if s == 1.0 {
        return f64::INFINITY;
    }
    if s == 0.0 {
        return -0.5;
    }
    if s < 0.0 && s == s.floor() && ((s as i64) % 2 == 0) {
        return 0.0;
    }
    hurwitz_zeta(s, 1.0)
}

// =========================================================================
// Elliptic integrals (Carlson forms)
// =========================================================================

/// Carlson's symmetric elliptic integral RF(x, y, z) for non-negative
/// arguments (at most one of which is zero), via the duplication theorem.
pub fn carlson_rf(mut x: f64, mut y: f64, mut z: f64) -> f64 {
    if x.min(y).min(z) < 0.0 {
        return f64::NAN;
    }
    const THIRD: f64 = 1.0 / 3.0;
    const C1: f64 = 1.0 / 24.0;
    const C2: f64 = 0.1;
    const C3: f64 = 3.0 / 44.0;
    const C4: f64 = 1.0 / 14.0;
    for _ in 0..100 {
        let sx = x.sqrt();
        let sy = y.sqrt();
        let sz = z.sqrt();
        let alamb = sx * (sy + sz) + sy * sz;
        x = 0.25 * (x + alamb);
        y = 0.25 * (y + alamb);
        z = 0.25 * (z + alamb);
        let ave = (x + y + z) * THIRD;
        let delx = (ave - x) / ave;
        let dely = (ave - y) / ave;
        let delz = (ave - z) / ave;
        if delx.abs().max(dely.abs()).max(delz.abs()) < 1e-9 {
            let e2 = delx * dely - delx * delz - dely * delz;
            let e3 = delx * dely * delz;
            return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / ave.sqrt();
        }
    }
    f64::NAN
}

/// Complete elliptic integral of the first kind, K(k) = ∫₀^{π/2} dθ /
/// √(1 − k²sin²θ), with modulus `k`. K(0) = π/2; K(±1) diverges.
///
/// Uses Carlson's RF duplication algorithm: K(k) = RF(0, 1−k², 1).
pub fn elliptic_k(k: f64) -> f64 {
    let m = k * k;
    if m >= 1.0 {
        return f64::INFINITY;
    }
    carlson_rf(0.0, 1.0 - m, 1.0)
}

/// K(k) and E(k) together via the arithmetic-geometric mean:
/// K = π/(2·a_N), E = K·(1 − Σ 2^{n−1}·(a_n² − b_n²)) with the AGM
/// sequence a₀ = 1, b₀ = k' = √(1−k²).
fn agm_ke(m: f64) -> (f64, f64) {
    let kp = (1.0 - m).sqrt();
    let mut a = 1.0f64;
    let mut b = kp;
    // n = 0 term: 2^{−1}·(a₀² − b₀²) = ½·(1 − k'²) = ½·m = ½·k²
    let mut sum = 0.5 * m;
    let mut pow2 = 1.0f64;
    for _ in 0..64 {
        let an = 0.5 * (a + b);
        let bn = (a * b).sqrt();
        // Cancellation-free difference of squares: a² − b² = (a−b)(a+b).
        let d2 = (an - bn) * (an + bn);
        a = an;
        b = bn;
        // When d2 is at rounding level relative to a², the sequence has
        // converged — adding further noise terms (amplified by pow2)
        // destroys the result.
        if d2 <= 1e-15 * a * a {
            break;
        }
        sum += pow2 * d2;
        pow2 *= 2.0;
    }
    let kk = std::f64::consts::FRAC_PI_2 / a;
    (kk, kk * (1.0 - sum))
}

/// Complete elliptic integral of the second kind, E(k) = ∫₀^{π/2}
/// √(1 − k²sin²θ) dθ, with modulus `k`. E(0) = π/2; E(±1) = 1.
///
/// Uses the AGM identity E(k) = K(k)·(1 − Σ 2^{n−1}(a_n² − b_n²)).
pub fn elliptic_e(k: f64) -> f64 {
    let m = k * k;
    if m > 1.0 {
        return f64::NAN;
    }
    if m >= 1.0 {
        return 1.0;
    }
    agm_ke(m).1
}

/// Incomplete elliptic integral of the first kind, F(φ, k) =
/// ∫₀^φ dθ / √(1 − k²sin²θ). F(π/2, k) = K(k).
///
/// Via Carlson's RF: F(φ, k) = sin φ · RF(cos²φ, 1 − k²sin²φ, 1).
pub fn elliptic_f(phi: f64, k: f64) -> f64 {
    let s = phi.sin();
    let c = phi.cos();
    carlson_rf(c * c, 1.0 - k * k * s * s, 1.0) * s
}

/// Incomplete elliptic integral of the second kind, E(φ, k) = ∫₀^φ
/// √(1 − k²sin²θ) dθ. E(π/2, k) = E(k).
///
/// Computed by composite Simpson quadrature of the smooth integrand.
pub fn elliptic_e_inc(phi: f64, k: f64) -> f64 {
    let m = k * k;
    crate::calculus::integrate_simpson(
        |th| (1.0 - m * th.sin().powi(2)).sqrt(),
        0.0,
        phi,
        2048,
    )
    .unwrap_or(f64::NAN)
}

#[cfg(test)]
mod zeta_gamma_tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn digamma_known_values() {
        // ψ(1) = −γ
        assert!(close(digamma(1.0), -0.5772156649015329, 1e-13));
        // ψ(0.5) = −γ − 2·ln 2
        assert!(close(digamma(0.5), -1.9635100260214235, 1e-13));
        // ψ(2) = 1 − γ
        assert!(close(digamma(2.0), 0.4227843350984671, 1e-13));
        // ψ(10) = 2.251752589066721…
        assert!(close(digamma(10.0), 2.2517525890667211, 1e-11));
        // reflection: ψ(−0.5) = ψ(1.5) + π/tan(π·0.5)... ψ(−0.5) = 2·ln 2 − γ
        assert!(close(digamma(-0.5), 0.03648997397857652, 1e-13));
        assert!(digamma(0.0).is_nan());
        assert!(digamma(-3.0).is_nan());
    }

    #[test]
    fn trigamma_known_values() {
        // ψ₁(1) = π²/6
        assert!(close(trigamma(1.0), std::f64::consts::PI.powi(2) / 6.0, 1e-12));
        // ψ₁(0.5) = π²/2
        assert!(close(trigamma(0.5), std::f64::consts::PI.powi(2) / 2.0, 1e-12));
        // ψ₁(2) = π²/6 − 1
        assert!(close(trigamma(2.0), std::f64::consts::PI.powi(2) / 6.0 - 1.0, 1e-12));
        assert!(trigamma(0.0).is_nan());
    }

    #[test]
    fn polygamma_consistency() {
        assert!(close(polygamma(0, 1.0), digamma(1.0), 1e-14));
        assert!(close(polygamma(1, 1.0), trigamma(1.0), 1e-14));
        // polygamma(2, 1) = −2·ζ(3)
        assert!(close(polygamma(2, 1.0), -2.0 * zeta(3.0), 1e-12));
        // polygamma(3, 1) = 6·ζ(4) = π⁴/15
        assert!(close(polygamma(3, 1.0), std::f64::consts::PI.powi(4) / 15.0, 1e-11));
        // derivative relation: d/dx ψ(x) ≈ ψ₁(x) — numeric check
        let x = 1.3;
        let h = 1e-6;
        let d = (digamma(x + h) - digamma(x - h)) / (2.0 * h);
        assert!(close(d, trigamma(x), 1e-6));
    }

    #[test]
    fn harmonic_known_values() {
        assert_eq!(harmonic(0), 0.0);
        assert_eq!(harmonic(1), 1.0);
        assert!(close(harmonic(2), 1.5, 1e-15));
        assert!(close(harmonic(10), 2.9289682539682538, 1e-14));
        // large n uses digamma: H_100 ≈ 5.187377517639621
        assert!(close(harmonic(100), 5.187377517639621, 1e-12));
        // asymptotic H_n ≈ ln n + γ + 1/(2n)
        let n = 50_000u64;
        let approx = (n as f64).ln() + std::f64::consts::EULER_GAMMA + 1.0 / (2.0 * n as f64);
        assert!(close(harmonic(n), approx, 1e-9));
    }

    #[test]
    fn hurwitz_zeta_reduces_to_riemann() {
        assert!(close(hurwitz_zeta(2.0, 1.0), zeta(2.0), 1e-13));
        assert!(close(hurwitz_zeta(3.0, 1.0), zeta(3.0), 1e-13));
        // ζ(s, 2) = ζ(s) − 1
        assert!(close(hurwitz_zeta(2.0, 2.0), zeta(2.0) - 1.0, 1e-12));
        assert!(close(hurwitz_zeta(2.0, 0.5), zeta(2.0) * 3.0, 1e-12)); // ζ(2, ½) = (2²−1)ζ(2)? = 3ζ(2)? π²/2 ✓ = 3·(π²/6)
    }

    #[test]
    fn zeta_known_values() {
        // ζ(2) = π²/6
        assert!(close(zeta(2.0), std::f64::consts::PI.powi(2) / 6.0, 1e-13));
        // ζ(4) = π⁴/90
        assert!(close(zeta(4.0), std::f64::consts::PI.powi(4) / 90.0, 1e-13));
        // ζ(3) = Apéry ≈ 1.2020569
        assert!(close(zeta(3.0), 1.2020569031595943, 1e-12));
        // ζ(0) = −1/2
        assert!(close(zeta(0.0), -0.5, 1e-15));
        // ζ(−1) = −1/12
        assert!(close(zeta(-1.0), -1.0 / 12.0, 1e-13));
        // ζ(−3) = 1/120
        assert!(close(zeta(-3.0), 1.0 / 120.0, 1e-13));
        // trivial zeros
        assert_eq!(zeta(-2.0), 0.0);
        assert_eq!(zeta(-4.0), 0.0);
        // ζ(½) ≈ −1.4603545088095868
        assert!(close(zeta(0.5), -1.4603545088095868, 1e-12));
        // pole
        assert_eq!(zeta(1.0), f64::INFINITY);
    }

    #[test]
    fn carlson_rf_known_values() {
        // RF(1, 1, 1) = 1
        assert!(close(carlson_rf(1.0, 1.0, 1.0), 1.0, 1e-15));
        // RF(0, 1, 1) = K(0) = π/2
        assert!(close(carlson_rf(0.0, 1.0, 1.0), std::f64::consts::FRAC_PI_2, 1e-14));
        // RF(0, ½, 1) = K(1/√2) = 1.8540746773013719
        assert!(close(carlson_rf(0.0, 0.5, 1.0), 1.8540746773013719, 1e-13));
    }

    #[test]
    fn elliptic_integrals_known_values() {
        use std::f64::consts::FRAC_PI_2;
        // K(0) = E(0) = π/2
        assert!(close(elliptic_k(0.0), FRAC_PI_2, 1e-14));
        assert!(close(elliptic_e(0.0), FRAC_PI_2, 1e-14));
        // K(½) = 1.6857503548125960 (modulus k = 0.5, parameter m = 0.25)
        assert!(close(elliptic_k(0.5), 1.6857503548125960, 1e-13));
        // E(½) = 1.4674622093394272
        assert!(close(elliptic_e(0.5), 1.4674622093394272, 1e-13));
        // K(0.9) = 2.2805491384227702
        assert!(close(elliptic_k(0.9), 2.2805491384227702, 1e-12));
        // E(±1) = 1; K(±1) diverges
        assert!(close(elliptic_e(1.0), 1.0, 1e-15));
        assert_eq!(elliptic_k(1.0), f64::INFINITY);
        // Legendre relation sanity: K(k)·E(k) ≥ π/2 for k ∈ (0,1)
        assert!(elliptic_k(0.7) * elliptic_e(0.7) > FRAC_PI_2);
    }

    #[test]
    fn incomplete_elliptic_reduces_to_complete() {
        let k = 0.6;
        let phi = std::f64::consts::FRAC_PI_2;
        assert!(close(elliptic_f(phi, k), elliptic_k(k), 1e-12));
        assert!(close(elliptic_e_inc(phi, k), elliptic_e(k), 1e-12));
        // φ = 0 → 0
        assert!(close(elliptic_f(0.0, k), 0.0, 1e-15));
        assert!(close(elliptic_e_inc(0.0, k), 0.0, 1e-15));
        // k = 0 → F = E = φ
        let phi = 0.7;
        assert!(close(elliptic_f(phi, 0.0), phi, 1e-14));
        assert!(close(elliptic_e_inc(phi, 0.0), phi, 1e-14));
    }

    #[test]
    fn special_functions_via_eval() {
        let ctx = crate::eval::Context::standard();
        let eval_str = |src: &str, ctx: &crate::eval::Context| -> f64 {
            crate::eval::eval(&crate::parser::Parser::parse(src).unwrap(), ctx).unwrap()
        };
        assert!(close(eval_str("zeta(2)", &ctx), std::f64::consts::PI.powi(2) / 6.0, 1e-12));
        assert!(close(eval_str("harmonic(10)", &ctx), 2.9289682539682538, 1e-14));
        assert!(close(eval_str("digamma(1) + 0.5772156649015329", &ctx), 0.0, 1e-13));
        assert!(close(eval_str("elliptic_k(0.5)", &ctx), 1.6857503548125960, 1e-12));
        assert!(close(eval_str("polygamma(2, 1)", &ctx), -2.0 * zeta(3.0), 1e-12));
    }

    // ===== Additional coverage =====

    #[test]
    fn digamma_large_value() {
        // ψ(100) ≈ ln(100) - 1/(2*100) (asymptotic)
        let v = digamma(100.0);
        let expected = 100.0_f64.ln() - 1.0 / 200.0;
        assert!(close(v, expected, 1e-4), "digamma(100) ≈ {}: got {}", expected, v);
    }

    #[test]
    fn digamma_fractional_negative() {
        // ψ(-1.5) via reflection
        let v = digamma(-1.5);
        assert!(v.is_finite(), "digamma(-1.5) should be finite: {}", v);
    }

    #[test]
    fn trigamma_negative_reflection() {
        // ψ₁(-0.5) via reflection formula: ψ₁(1-x) + ψ₁(x) = π²/sin²(πx)
        // With x = 1.5: ψ₁(-0.5) + ψ₁(1.5) = π²/sin²(3π/2) = π²
        // ψ₁(1.5) = ψ₁(0.5) - 1/0.5² = π²/2 - 4
        // So ψ₁(-0.5) = π² - (π²/2 - 4) = π²/2 + 4
        let v = trigamma(-0.5);
        assert!(v.is_finite(), "trigamma(-0.5) should be finite: {}", v);
        let expected = std::f64::consts::PI.powi(2) / 2.0 + 4.0;
        assert!(close(v, expected, 1e-6), "trigamma(-0.5) ≈ {}: got {}", expected, v);
    }

    #[test]
    fn trigamma_large_value() {
        // ψ₁(100) ≈ 1/100 + 1/(2*100²) (asymptotic)
        let v = trigamma(100.0);
        let expected = 1.0 / 100.0 + 1.0 / (2.0 * 100.0 * 100.0);
        assert!(close(v, expected, 1e-4), "trigamma(100) ≈ {}: got {}", expected, v);
    }

    #[test]
    fn polygamma_recurrence_identity() {
        // ψ^(m)(x+1) = ψ^(m)(x) + (-1)^m * m! / x^(m+1)
        let m: u32 = 1;
        let x = 2.0;
        let lhs = polygamma(m, x + 1.0);
        let mi = m as i32;
        let rhs = polygamma(m, x) + (-1.0_f64).powi(mi) * (1.0_f64).powi(mi) / x.powi(mi + 1);
        // For m=1: trigamma(3) = trigamma(2) - 1/4
        assert!(close(lhs, rhs, 1e-10), "polygamma recurrence: {} vs {}", lhs, rhs);
    }

    #[test]
    fn harmonic_consistency_with_digamma() {
        // H_n = ψ(n+1) + γ for all n
        for n in [1, 5, 10, 50, 100, 500] {
            let h = harmonic(n);
            let psi = digamma((n + 1) as f64) + 0.5772156649015329;
            assert!(close(h, psi, 1e-6), "H_{} = {} but ψ({})+γ = {}", n, h, n + 1, psi);
        }
    }

    #[test]
    fn zeta_negative_odd() {
        // ζ(-5) = -1/252
        let v = zeta(-5.0);
        assert!(close(v, -1.0 / 252.0, 1e-6), "zeta(-5) = {}: got {}", -1.0 / 252.0, v);

        // ζ(-7) = 1/240
        let v2 = zeta(-7.0);
        assert!(close(v2, 1.0 / 240.0, 1e-6), "zeta(-7) = {}: got {}", 1.0 / 240.0, v2);
    }

    #[test]
    fn zeta_even_positive() {
        // ζ(6) = π^6/945
        let v = zeta(6.0);
        let expected = std::f64::consts::PI.powi(6) / 945.0;
        assert!(close(v, expected, 1e-10), "zeta(6) = {}: got {}", expected, v);

        // ζ(8) = π^8/9450
        let v2 = zeta(8.0);
        let expected2 = std::f64::consts::PI.powi(8) / 9450.0;
        assert!(close(v2, expected2, 1e-10), "zeta(8) = {}: got {}", expected2, v2);
    }

    #[test]
    fn hurwitz_zeta_non_integer_a() {
        // ζ(2, 0.5) = π²/2 (known value)
        let v = hurwitz_zeta(2.0, 0.5);
        let expected = std::f64::consts::PI.powi(2) / 2.0;
        assert!(close(v, expected, 1e-6), "hurwitz_zeta(2, 0.5) = {}: got {}", expected, v);
    }

    #[test]
    fn hurwitz_zeta_s_equals_one_errors() {
        // s=1 is the pole of zeta
        let v = hurwitz_zeta(1.0, 1.0);
        assert!(v.is_nan() || v.is_infinite(), "hurwitz_zeta(1, 1) should be NaN/inf: {}", v);
    }

    #[test]
    fn carlson_rf_symmetry() {
        // RF(x, y, z) is symmetric in its arguments
        let a = carlson_rf(1.0, 2.0, 3.0);
        let b = carlson_rf(3.0, 2.0, 1.0);
        let c = carlson_rf(2.0, 3.0, 1.0);
        assert!(close(a, b, 1e-12), "RF symmetry: {} vs {}", a, b);
        assert!(close(a, c, 1e-12), "RF symmetry: {} vs {}", a, c);
    }

    #[test]
    fn carlson_rf_with_one_zero() {
        // RF(0, 1, 1) = π/2 (by the integral definition)
        let v = carlson_rf(0.0, 1.0, 1.0);
        let expected = std::f64::consts::FRAC_PI_2;
        assert!(close(v, expected, 1e-10), "RF(0,1,1) = {}: got {}", expected, v);
    }

    #[test]
    fn elliptic_k_at_zero() {
        // K(0) = π/2
        let v = elliptic_k(0.0);
        assert!(close(v, std::f64::consts::FRAC_PI_2, 1e-12), "K(0) = π/2: got {}", v);
    }

    #[test]
    fn elliptic_e_at_zero() {
        // E(0) = π/2
        let v = elliptic_e(0.0);
        assert!(close(v, std::f64::consts::FRAC_PI_2, 1e-12), "E(0) = π/2: got {}", v);
    }

    #[test]
    fn elliptic_e_at_one() {
        // E(1) = 1
        let v = elliptic_e(1.0);
        assert!(close(v, 1.0, 1e-6), "E(1) = 1: got {}", v);
    }

    #[test]
    fn elliptic_f_at_zero_modulus() {
        // F(φ, 0) = φ
        let v = elliptic_f(std::f64::consts::FRAC_PI_4, 0.0);
        assert!(close(v, std::f64::consts::FRAC_PI_4, 1e-12), "F(π/4, 0) = π/4: got {}", v);
    }

    #[test]
    fn elliptic_e_inc_at_zero_modulus() {
        // E_inc(φ, 0) = φ
        let v = elliptic_e_inc(std::f64::consts::FRAC_PI_4, 0.0);
        assert!(close(v, std::f64::consts::FRAC_PI_4, 1e-6), "E_inc(π/4, 0) = π/4: got {}", v);
    }

    #[test]
    fn elliptic_e_inc_at_quarter_pi() {
        // E_inc(π/2, k) = E(k) (complete elliptic integral)
        let v = elliptic_e_inc(std::f64::consts::FRAC_PI_2, 0.5);
        let expected = elliptic_e(0.5);
        assert!(close(v, expected, 1e-4), "E_inc(π/2, 0.5) = E(0.5) = {}: got {}", expected, v);
    }

    #[test]
    fn special_functions_via_eval_extended() {
        let ctx = crate::eval::Context::standard();
        let eval_str = |src: &str, ctx: &crate::eval::Context| -> f64 {
            crate::eval::eval(&crate::parser::Parser::parse(src).unwrap(), ctx).unwrap()
        };
        // trigamma
        assert!(close(eval_str("trigamma(1)", &ctx), std::f64::consts::PI.powi(2) / 6.0, 1e-10));
        // elliptic_e
        assert!(close(eval_str("elliptic_e(0)", &ctx), std::f64::consts::FRAC_PI_2, 1e-10));
        // elliptic_f
        assert!(close(eval_str("elliptic_f(0, 0.5)", &ctx), 0.0, 1e-10));
    }

    #[test]
    fn bessel_jn_recurrence_large_n() {
        // J_n(0) = 0 for n > 0, J_0(0) = 1
        assert!(close(bessel_jn(10, 0.0), 0.0, 1e-15));
        assert!(close(bessel_jn(0, 0.0), 1.0, 1e-15));
    }

    #[test]
    fn gamma_reflection_negative() {
        // Γ(1-x)Γ(x) = π/sin(πx) for x = 0.5: Γ(0.5)² = π
        let g = gamma(0.5);
        assert!(close(g * g, std::f64::consts::PI, 1e-10), "Γ(0.5)² = π: got {}", g * g);
    }

    #[test]
    fn beta_symmetry() {
        // B(x, y) = B(y, x)
        let a = beta(2.0, 3.0);
        let b = beta(3.0, 2.0);
        assert!(close(a, b, 1e-12), "B(2,3) = B(3,2): {} vs {}", a, b);
    }

    #[test]
    fn erf_at_inf() {
        // erf(∞) = 1, erf(-∞) = -1
        assert!(close(erf(50.0), 1.0, 1e-15), "erf(50) ≈ 1: got {}", erf(50.0));
        assert!(close(erf(-50.0), -1.0, 1e-15), "erf(-50) ≈ -1: got {}", erf(-50.0));
    }

    #[test]
    fn erfc_at_inf() {
        // erfc(∞) = 0, erfc(-∞) = 2
        assert!(close(erfc(50.0), 0.0, 1e-15), "erfc(50) ≈ 0: got {}", erfc(50.0));
        assert!(close(erfc(-50.0), 2.0, 1e-15), "erfc(-50) ≈ 2: got {}", erfc(-50.0));
    }

    #[test]
    fn incomplete_gamma_chi2_round_trip() {
        // For a chi-squared distribution with k degrees of freedom,
        // P(k/2, x/2) is the CDF. Check P(1, x) = 1 - e^(-x)
        for &x in &[0.5, 1.0, 2.0, 5.0] {
            let p = incomplete_gamma_p(1.0, x);
            let expected = 1.0 - (-x).exp();
            assert!(close(p, expected, 1e-10), "P(1, {}) = {}: got {}", x, expected, p);
        }
    }
}
