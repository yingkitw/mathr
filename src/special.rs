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
