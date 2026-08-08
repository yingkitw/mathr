//! Arbitrary-precision integer arithmetic for number theory.
//!
//! Uses [`num_bigint::BigInt`] for unbounded integer operations.
//! Provides big-integer versions of primality testing, factorization,
//! GCD, factorial, Fibonacci, and modular exponentiation.

use crate::error::{MathError, Result};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// Parse a decimal string into a `BigInt`.
pub fn parse(s: &str) -> Result<BigInt> {
    s.trim()
        .parse::<BigInt>()
        .map_err(|_| MathError::Parse(format!("invalid big integer: {}", s)))
}

/// Greatest common divisor of two big integers.
pub fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    if a.is_zero() && b.is_zero() {
        return BigInt::zero();
    }
    let mut x = a.abs();
    let mut y = b.abs();
    while !y.is_zero() {
        let t = y.clone();
        y = x % &t;
        x = t;
    }
    x
}

/// Least common multiple of two big integers.
pub fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    if a.is_zero() || b.is_zero() {
        return BigInt::zero();
    }
    (a / gcd(a, b)) * b
}

/// Modular exponentiation: `base^exp mod m`.
pub fn mod_pow(base: &BigInt, exp: &BigInt, m: &BigInt) -> Result<BigInt> {
    if m.is_zero() {
        return Err(MathError::Eval("modular exponentiation: modulus is zero".into()));
    }
    if m.is_one() {
        return Ok(BigInt::zero());
    }
    if exp.is_negative() {
        return Err(MathError::Eval(
            "modular exponentiation: negative exponent not supported for big integers".into(),
        ));
    }
    let m_abs = m.abs();
    let mut result = BigInt::one();
    let mut base = base.mod_floor(&m_abs);
    let mut exp = exp.clone();

    while !exp.is_zero() {
        if exp.is_odd() {
            result = (result * &base).mod_floor(&m_abs);
        }
        exp >>= 1;
        base = (&base * &base).mod_floor(&m_abs);
    }
    Ok(result)
}

/// Miller–Rabin primality test for big integers.
/// Uses deterministic witnesses for n < 3.3e24, then probabilistic with `k` rounds.
pub fn is_prime(n: &BigInt, k: usize) -> bool {
    let two = BigInt::from(2);
    let three = BigInt::from(3);

    if n <= &BigInt::one() {
        return false;
    }
    if n == &two || n == &three {
        return true;
    }
    if n.is_even() {
        return false;
    }
    let _ = &two; // referenced above

    // Write n-1 = d * 2^r with d odd
    let n_minus_1 = n - BigInt::one();
    let mut d = n_minus_1.clone();
    let mut r = 0u32;
    while d.is_even() {
        d >>= 1;
        r += 1;
    }

    // Deterministic witnesses for n < 3,317,044,064,679,887,385,961,981
    let deterministic: &[u64] = &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    let deterministic_bound = BigInt::from(3317044064679887385u64)
        * BigInt::from(1000000000u64)
        * BigInt::from(1000000u64)
        * BigInt::from(1000u64);

    let witnesses: Vec<BigInt> = if n < &deterministic_bound {
        deterministic.iter().map(|&w| BigInt::from(w)).collect()
    } else {
        // Probabilistic: use first k primes as witnesses
        deterministic
            .iter()
            .take(k)
            .map(|&w| BigInt::from(w))
            .collect()
    };

    for a in witnesses {
        if &a >= n {
            continue;
        }
        let mut x = mod_pow(&a, &d, n).unwrap_or(BigInt::zero());
        if x == BigInt::one() || x == n_minus_1 {
            continue;
        }
        let mut composite = true;
        for _ in 0..r - 1 {
            x = (&x * &x).mod_floor(n);
            if x == n_minus_1 {
                composite = false;
                break;
            }
        }
        if composite {
            return false;
        }
    }
    true
}

/// Factorize a big integer using trial division up to a limit, then Pollard's rho.
/// Returns a sorted list of (prime, exponent) pairs.
pub fn factorize(n: &BigInt) -> Vec<(BigInt, u32)> {
    if n <= &BigInt::one() {
        return vec![];
    }

    let mut n = n.abs();
    let mut factors: Vec<(BigInt, u32)> = Vec::new();

    // Trial division by small primes
    let small_primes: &[u64] = &[
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
        89, 97,
    ];
    for &p in small_primes {
        let p_big = BigInt::from(p);
        if &p_big * &p_big > n {
            break;
        }
        let mut count = 0u32;
        while (&n % &p_big).is_zero() {
            n /= &p_big;
            count += 1;
        }
        if count > 0 {
            factors.push((p_big, count));
        }
    }

    // Pollard's rho for remaining factors
    while n > BigInt::one() {
        if is_prime(&n, 20) {
            factors.push((n.clone(), 1));
            break;
        }
        if let Some(factor) = pollard_rho(&n) {
            // Recursively factor the factor and the cofactor
            let cofactor = &n / &factor;
            let mut sub_factors = factorize(&factor);
            // Merge sub-factors
            for (p, e) in sub_factors.drain(..) {
                if let Some(existing) = factors.iter_mut().find(|(ep, _)| ep == &p) {
                    existing.1 += e;
                } else {
                    factors.push((p, e));
                }
            }
            n = cofactor;
        } else {
            // Fallback: n is prime (pollard_rho failed)
            factors.push((n.clone(), 1));
            break;
        }
    }

    factors.sort_by(|a, b| a.0.cmp(&b.0));
    factors
}

/// Pollard's rho algorithm for finding a non-trivial factor.
fn pollard_rho(n: &BigInt) -> Option<BigInt> {
    if n.is_even() {
        return Some(BigInt::from(2));
    }
    if is_prime(n, 10) {
        return None;
    }

    let mut c = BigInt::from(1);
    loop {
        let mut x = BigInt::from(2);
        let mut y = BigInt::from(2);
        let mut d = BigInt::one();

        while d.is_one() {
            x = pollard_f(&x, &c, n);
            y = pollard_f(&pollard_f(&y, &c, n), &c, n);
            let diff = if x > y { &x - &y } else { &y - &x };
            if diff.is_zero() {
                break; // cycle detected, try new c
            }
            d = gcd(&diff, n);
        }

        if d < *n && d > BigInt::one() {
            return Some(d);
        }
        if d == *n {
            c += BigInt::one();
            if c > BigInt::from(100) {
                return None;
            }
        }
    }
}

/// f(x) = x^2 + c (mod n)
fn pollard_f(x: &BigInt, c: &BigInt, n: &BigInt) -> BigInt {
    ((x * x + c) % n).into()
}

/// Factorial of n: n! = 1 * 2 * ... * n.
pub fn factorial(n: u64) -> BigInt {
    let mut result = BigInt::one();
    for i in 2..=n {
        result *= BigInt::from(i);
    }
    result
}

/// Fibonacci number F(n). Uses fast doubling.
pub fn fibonacci(n: u64) -> BigInt {
    if n == 0 {
        return BigInt::zero();
    }
    fib_pair(n).0
}

/// Returns (F(n), F(n+1)) using fast doubling.
fn fib_pair(n: u64) -> (BigInt, BigInt) {
    if n == 0 {
        return (BigInt::zero(), BigInt::one());
    }
    let (a, b) = fib_pair(n / 2);
    // F(2k) = F(k) * (2*F(k+1) - F(k))
    let f2k = &a * ((&b + &b) - &a);
    // F(2k+1) = F(k)^2 + F(k+1)^2
    let f2k1 = &a * &a + &b * &b;
    if n % 2 == 0 {
        (f2k, f2k1)
    } else {
        let next = &f2k + &f2k1;
        (f2k1, next)
    }
}

/// Binomial coefficient C(n, k) = n! / (k! * (n-k)!).
pub fn binomial(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let k = k.min(n - k);
    let mut result = BigInt::one();
    for i in 0..k {
        result *= BigInt::from(n - i);
        result /= BigInt::from(i + 1);
    }
    result
}

/// Euler's totient function φ(n) for big integers.
/// (Named `totient` to match SymPy's `sympy.ntheory.totient`.)
pub fn totient(n: &BigInt) -> BigInt {
    if n <= &BigInt::zero() {
        return BigInt::zero();
    }
    let factors = factorize(n);
    let mut result = n.abs();
    for (p, _) in &factors {
        result = &result / p * (p - BigInt::one());
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn big_gcd() {
        assert_eq!(gcd(&BigInt::from(12), &BigInt::from(8)), BigInt::from(4));
        assert_eq!(gcd(&BigInt::from(0), &BigInt::from(5)), BigInt::from(5));
        assert_eq!(
            gcd(&BigInt::from(-12), &BigInt::from(8)),
            BigInt::from(4)
        );
    }

    #[test]
    fn big_lcm() {
        assert_eq!(lcm(&BigInt::from(4), &BigInt::from(6)), BigInt::from(12));
        assert_eq!(lcm(&BigInt::from(0), &BigInt::from(5)), BigInt::from(0));
    }

    #[test]
    fn big_mod_pow() {
        assert_eq!(
            mod_pow(&BigInt::from(2), &BigInt::from(10), &BigInt::from(1000)).unwrap(),
            BigInt::from(24)
        );
        assert_eq!(
            mod_pow(&BigInt::from(3), &BigInt::from(7), &BigInt::from(13)).unwrap(),
            BigInt::from(3)
        );
        // Large: 2^100 mod 10^20
        let result = mod_pow(
            &BigInt::from(2),
            &BigInt::from(100),
            &BigInt::from(100000000000000000000u128),
        )
        .unwrap();
        assert!(result < BigInt::from(100000000000000000000u128));
    }

    #[test]
    fn big_mod_pow_errors() {
        assert!(mod_pow(&BigInt::from(2), &BigInt::from(10), &BigInt::zero()).is_err());
        assert!(mod_pow(&BigInt::from(2), &BigInt::from(-1), &BigInt::from(7)).is_err());
    }

    #[test]
    fn big_is_prime_small() {
        assert!(is_prime(&BigInt::from(2), 10));
        assert!(is_prime(&BigInt::from(3), 10));
        assert!(is_prime(&BigInt::from(7), 10));
        assert!(is_prime(&BigInt::from(97), 10));
        assert!(!is_prime(&BigInt::from(1), 10));
        assert!(!is_prime(&BigInt::from(4), 10));
        assert!(!is_prime(&BigInt::from(100), 10));
    }

    #[test]
    fn big_is_prime_large() {
        // Known large primes
        assert!(is_prime(&BigInt::from(1000000007u64), 20));
        assert!(is_prime(&BigInt::from(1000000009u64), 20));
        assert!(is_prime(&BigInt::from(999999999989u64), 20));
        // Known composites
        assert!(!is_prime(&BigInt::from(1000000008u64), 20));
        assert!(!is_prime(&BigInt::from(999999999990u64), 20));
    }

    #[test]
    fn big_is_prime_very_large() {
        // 2^61 - 1 is a Mersenne prime
        let mersenne = BigInt::from(1u64 << 61) - BigInt::one();
        assert!(is_prime(&mersenne, 20));
    }

    #[test]
    fn big_factorize() {
        let factors = factorize(&BigInt::from(12));
        assert_eq!(factors, vec![
            (BigInt::from(2), 2),
            (BigInt::from(3), 1)
        ]);

        let factors = factorize(&BigInt::from(360));
        assert_eq!(factors, vec![
            (BigInt::from(2), 3),
            (BigInt::from(3), 2),
            (BigInt::from(5), 1)
        ]);

        let factors = factorize(&BigInt::from(1));
        assert!(factors.is_empty());

        let factors = factorize(&BigInt::from(17));
        assert_eq!(factors, vec![(BigInt::from(17), 1)]);
    }

    #[test]
    fn big_factorize_large() {
        // 1000000007 is prime
        let factors = factorize(&BigInt::from(1000000007u64));
        assert_eq!(factors, vec![(BigInt::from(1000000007u64), 1)]);

        // Semi-prime: 1000000007 * 1000000009
        let n = BigInt::from(1000000007u64) * BigInt::from(1000000009u64);
        let factors = factorize(&n);
        assert_eq!(factors.len(), 2);
        assert_eq!(factors[0].0, BigInt::from(1000000007u64));
        assert_eq!(factors[1].0, BigInt::from(1000000009u64));
    }

    #[test]
    fn big_factorial() {
        assert_eq!(factorial(0), BigInt::from(1));
        assert_eq!(factorial(1), BigInt::from(1));
        assert_eq!(factorial(5), BigInt::from(120));
        assert_eq!(factorial(10), BigInt::from(3628800));
        // 20! overflows u64
        let f20 = factorial(20);
        assert_eq!(f20.to_string(), "2432902008176640000");
        // 50! — way beyond u64
        let f50 = factorial(50);
        assert!(f50.to_string().len() > 60);
    }

    #[test]
    fn big_fibonacci() {
        assert_eq!(fibonacci(0), BigInt::from(0));
        assert_eq!(fibonacci(1), BigInt::from(1));
        assert_eq!(fibonacci(2), BigInt::from(1));
        assert_eq!(fibonacci(10), BigInt::from(55));
        assert_eq!(fibonacci(20), BigInt::from(6765));
        // F(100) — way beyond u64
        let f100 = fibonacci(100);
        assert_eq!(f100.to_string(), "354224848179261915075");
    }

    #[test]
    fn big_binomial() {
        assert_eq!(binomial(0, 0), BigInt::from(1));
        assert_eq!(binomial(5, 0), BigInt::from(1));
        assert_eq!(binomial(5, 2), BigInt::from(10));
        assert_eq!(binomial(10, 3), BigInt::from(120));
        assert_eq!(binomial(5, 6), BigInt::from(0)); // k > n
        // C(100, 50) — way beyond u64
        let c = binomial(100, 50);
        assert_eq!(c.to_string(), "100891344545564193334812497256");
    }

    #[test]
    fn big_totient() {
        assert_eq!(totient(&BigInt::from(1)), BigInt::from(1));
        assert_eq!(totient(&BigInt::from(7)), BigInt::from(6)); // prime
        assert_eq!(totient(&BigInt::from(12)), BigInt::from(4)); // 12 = 2^2 * 3
        assert_eq!(totient(&BigInt::from(36)), BigInt::from(12)); // 36 = 2^2 * 3^2
    }

    #[test]
    fn big_parse() {
        assert_eq!(parse("42").unwrap(), BigInt::from(42));
        assert_eq!(parse("-17").unwrap(), BigInt::from(-17));
        assert_eq!(
            parse("123456789012345678901234567890").unwrap(),
            BigInt::from(123456789012345678901234567890i128)
        );
        assert!(parse("abc").is_err());
    }
}
