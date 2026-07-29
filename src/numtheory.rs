//! Number theory from scratch.
//!
//! Provides GCD, LCM, primality testing, prime factorization,
//! binomial coefficients, factorial, and Fibonacci numbers.

use crate::error::{MathError, Result};

/// Greatest common divisor (Euclidean algorithm).
pub fn gcd(a: u64, b: u64) -> u64 {
    let mut a = a;
    let mut b = b;
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Least common multiple.
pub fn lcm(a: u64, b: u64) -> u64 {
    if a == 0 || b == 0 {
        return 0;
    }
    a / gcd(a, b) * b
}

/// Extended GCD: returns (g, x, y) such that `a*x + b*y = g = gcd(a, b)`.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        return (a, 1, 0);
    }
    let (g, x1, y1) = extended_gcd(b, a % b);
    (g, y1, x1 - (a / b) * y1)
}

/// Modular inverse: returns `x` such that `a*x ≡ 1 (mod m)`, if it exists.
pub fn mod_inverse(a: i64, m: i64) -> Option<i64> {
    let (g, x, _) = extended_gcd(((a % m) + m) % m, m);
    if g != 1 {
        None
    } else {
        Some(((x % m) + m) % m)
    }
}

/// Trial-division primality test. Good enough for `n < 2^32`.
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    if n < 4 {
        return true;
    }
    if n % 2 == 0 || n % 3 == 0 {
        return false;
    }
    let mut i = 5u64;
    while i * i <= n {
        if n % i == 0 || n % (i + 2) == 0 {
            return false;
        }
        i += 6;
    }
    true
}

/// Prime factorization via trial division. Returns factors in ascending order.
pub fn prime_factors(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    while n % 2 == 0 {
        factors.push(2);
        n /= 2;
    }
    let mut i = 3u64;
    while i * i <= n {
        while n % i == 0 {
            factors.push(i);
            n /= i;
        }
        i += 2;
    }
    if n > 1 {
        factors.push(n);
    }
    factors
}

/// Binomial coefficient C(n, k) = n! / (k! * (n-k)!).
/// Uses iterative computation to avoid overflow for moderate n.
pub fn binomial(n: u64, k: u64) -> Result<u64> {
    if k > n {
        return Ok(0);
    }
    let k = k.min(n - k);
    let mut result: u64 = 1;
    for i in 0..k {
        // result = result * (n - i) / (i + 1)
        // Multiply first, then divide — this stays exact because
        // C(n, i+1) = C(n, i) * (n-i) / (i+1) is always an integer.
        result = result
            .checked_mul(n - i)
            .ok_or_else(|| MathError::InvalidArgument("binomial: overflow".into()))?;
        result /= i + 1;
    }
    Ok(result)
}

/// Factorial n! computed iteratively.
pub fn factorial(n: u64) -> Result<u64> {
    let mut result: u64 = 1;
    for i in 2..=n {
        result = result
            .checked_mul(i)
            .ok_or_else(|| MathError::InvalidArgument(format!("factorial: overflow at {}", i)))?;
    }
    Ok(result)
}

/// nth Fibonacci number (F(0) = 0, F(1) = 1) via fast doubling.
pub fn fibonacci(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    // fast doubling: F(2k) = F(k)*(2*F(k+1) - F(k)), F(2k+1) = F(k+1)^2 + F(k)^2
    fn fib(n: u64) -> (u64, u64) {
        if n == 0 {
            return (0, 1);
        }
        let (a, b) = fib(n / 2);
        let c = a * (2 * b - a);
        let d = a * a + b * b;
        if n % 2 == 0 {
            (c, d)
        } else {
            (d, c + d)
        }
    }
    fib(n).0
}

/// List all primes up to `n` using the Sieve of Eratosthenes.
pub fn sieve_primes(n: u64) -> Vec<u64> {
    if n < 2 {
        return Vec::new();
    }
    let n = n as usize;
    let mut is_composite = vec![false; n + 1];
    let mut primes = Vec::new();
    for i in 2..=n {
        if !is_composite[i] {
            primes.push(i as u64);
            let mut j = i * i;
            while j <= n {
                is_composite[j] = true;
                j += i;
            }
        }
    }
    primes
}

/// Euler's totient function φ(n): count of integers 1..=n coprime to n.
pub fn euler_totient(n: u64) -> u64 {
    if n == 0 {
        return 0;
    }
    let mut result = n;
    let mut m = n;
    let mut p = 2u64;
    while p * p <= m {
        if m % p == 0 {
            while m % p == 0 {
                m /= p;
            }
            result -= result / p;
        }
        p += 1;
    }
    if m > 1 {
        result -= result / m;
    }
    result
}

// --- Modular exponentiation -------------------------------------------------

/// Modular exponentiation: `base^exp mod m` using binary exponentiation.
pub fn mod_pow(base: u64, exp: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut result: u128 = 1;
    let mut base: u128 = (base % m) as u128;
    let m = m as u128;
    let mut exp = exp;
    while exp > 0 {
        if exp & 1 == 1 {
            result = result * base % m;
        }
        exp >>= 1;
        base = base * base % m;
    }
    result as u64
}

/// Miller–Rabin probabilistic primality test.
///
/// `k` is the number of rounds; higher = more accurate.
/// For `k` rounds, the error probability is at most `4^(-k)`.
/// With `k = 20`, the result is deterministic for all `n < 3.3 × 10^24`.
pub fn is_prime_miller_rabin(n: u64, k: usize) -> bool {
    if n < 2 {
        return false;
    }
    if n == 2 || n == 3 {
        return true;
    }
    if n % 2 == 0 {
        return false;
    }

    // Write n-1 as 2^r * d
    let mut d = n - 1;
    let mut r = 0u32;
    while d % 2 == 0 {
        d /= 2;
        r += 1;
    }

    // Deterministic witnesses for n < 3.3e24
    let witnesses: &[u64] = if n < 2047 {
        &[2]
    } else if n < 1_373_653 {
        &[2, 3]
    } else if n < 9_080_191 {
        &[31, 73]
    } else if n < 25_326_001 {
        &[2, 3, 5]
    } else if n < 3_215_031_751 {
        &[2, 3, 5, 7]
    } else if n < 4_759_123_141 {
        &[2, 7, 61]
    } else if n < 1_122_004_669_633 {
        &[2, 13, 23, 1662803]
    } else if n < 2_152_302_898_747 {
        &[2, 3, 5, 7, 11]
    } else if n < 3_474_749_660_383 {
        &[2, 3, 5, 7, 11, 13]
    } else if n < 341_550_071_728_321 {
        &[2, 3, 5, 7, 11, 13, 17]
    } else {
        // For very large n, use first k primes as witnesses
        &[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    };

    let witnesses: Vec<u64> = if (n as u128) < 3_317_044_064_679_887_385_961_981 {
        witnesses.to_vec()
    } else {
        // Random-ish witnesses: use first k primes
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71];
        primes.iter().take(k.max(5)).copied().collect()
    };

    'witness: for &a in &witnesses {
        if a >= n {
            continue;
        }
        let mut x = mod_pow(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..(r - 1) {
            x = mod_pow(x, 2, n);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// Chinese Remainder Theorem: given pairwise coprime moduli and remainders,
/// returns x such that x ≡ r_i (mod m_i) for all i.
pub fn chinese_remainder(remainders: &[u64], moduli: &[u64]) -> Result<u64> {
    if remainders.len() != moduli.len() {
        return Err(MathError::InvalidArgument("chinese_remainder: length mismatch".into()));
    }
    if remainders.is_empty() {
        return Err(MathError::InvalidArgument("chinese_remainder: empty input".into()));
    }
    // Check pairwise coprimality
    for i in 0..moduli.len() {
        for j in (i + 1)..moduli.len() {
            if gcd(moduli[i], moduli[j]) != 1 {
                return Err(MathError::InvalidArgument(format!(
                    "chinese_remainder: moduli {} and {} are not coprime",
                    moduli[i], moduli[j]
                )));
            }
        }
    }
    let m_prod: u64 = moduli.iter().product();
    let mut x: u64 = 0;
    for i in 0..remainders.len() {
        let mi = moduli[i];
        let mi_prod = m_prod / mi;
        let inv = mod_inverse(mi_prod as i64, mi as i64)
            .ok_or_else(|| MathError::InvalidArgument("chinese_remainder: no inverse".into()))?;
        x = (x + (remainders[i] as u128 * mi_prod as u128 % m_prod as u128 * inv as u128 % m_prod as u128) as u64) % m_prod;
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gcd_basic() {
        assert_eq!(gcd(12, 18), 6);
        assert_eq!(gcd(7, 13), 1);
        assert_eq!(gcd(0, 5), 5);
    }

    #[test]
    fn lcm_basic() {
        assert_eq!(lcm(4, 6), 12);
        assert_eq!(lcm(5, 7), 35);
        assert_eq!(lcm(0, 5), 0);
    }

    #[test]
    fn extended_gcd_bezout() {
        let (g, x, y) = extended_gcd(35, 15);
        assert_eq!(g, 5);
        assert_eq!(35 * x + 15 * y, 5);
    }

    #[test]
    fn mod_inverse_works() {
        let inv = mod_inverse(3, 11).unwrap();
        assert_eq!((3 * inv) % 11, 1);
        assert!(mod_inverse(4, 8).is_none());
    }

    #[test]
    fn primality() {
        assert!(!is_prime(0));
        assert!(!is_prime(1));
        assert!(is_prime(2));
        assert!(is_prime(3));
        assert!(!is_prime(4));
        assert!(is_prime(17));
        assert!(is_prime(97));
        assert!(!is_prime(100));
        assert!(is_prime(2147483647));
    }

    #[test]
    fn prime_factors_basic() {
        assert_eq!(prime_factors(12), vec![2, 2, 3]);
        assert_eq!(prime_factors(17), vec![17]);
        assert_eq!(prime_factors(60), vec![2, 2, 3, 5]);
        assert_eq!(prime_factors(1), vec![]);
    }

    #[test]
    fn binomial_basic() {
        assert_eq!(binomial(5, 0).unwrap(), 1);
        assert_eq!(binomial(5, 2).unwrap(), 10);
        assert_eq!(binomial(10, 3).unwrap(), 120);
        assert_eq!(binomial(5, 6).unwrap(), 0);
    }

    #[test]
    fn factorial_basic() {
        assert_eq!(factorial(0).unwrap(), 1);
        assert_eq!(factorial(1).unwrap(), 1);
        assert_eq!(factorial(5).unwrap(), 120);
        assert_eq!(factorial(10).unwrap(), 3628800);
    }

    #[test]
    fn fibonacci_basic() {
        assert_eq!(fibonacci(0), 0);
        assert_eq!(fibonacci(1), 1);
        assert_eq!(fibonacci(2), 1);
        assert_eq!(fibonacci(10), 55);
        assert_eq!(fibonacci(20), 6765);
        assert_eq!(fibonacci(50), 12586269025);
    }

    #[test]
    fn sieve_basic() {
        let primes = sieve_primes(20);
        assert_eq!(primes, vec![2, 3, 5, 7, 11, 13, 17, 19]);
    }

    #[test]
    fn totient_basic() {
        assert_eq!(euler_totient(1), 1);
        assert_eq!(euler_totient(9), 6);
        assert_eq!(euler_totient(10), 4);
        assert_eq!(euler_totient(36), 12);
    }

    #[test]
    fn mod_pow_basic() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 5, 7), 5);
        assert_eq!(mod_pow(7, 0, 11), 1);
        assert_eq!(mod_pow(2, 32, 1), 0);
    }

    #[test]
    fn miller_rabin_small_primes() {
        assert!(is_prime_miller_rabin(2, 10));
        assert!(is_prime_miller_rabin(3, 10));
        assert!(is_prime_miller_rabin(5, 10));
        assert!(is_prime_miller_rabin(7, 10));
        assert!(is_prime_miller_rabin(97, 10));
        assert!(is_prime_miller_rabin(2147483647, 10));
    }

    #[test]
    fn miller_rabin_composites() {
        assert!(!is_prime_miller_rabin(1, 10));
        assert!(!is_prime_miller_rabin(4, 10));
        assert!(!is_prime_miller_rabin(9, 10));
        assert!(!is_prime_miller_rabin(15, 10));
        assert!(!is_prime_miller_rabin(100, 10));
        assert!(!is_prime_miller_rabin(561, 10)); // Carmichael number
        assert!(!is_prime_miller_rabin(1729, 10)); // Carmichael number
    }

    #[test]
    fn miller_rabin_large_prime() {
        // 2^61 - 1 is a Mersenne prime
        assert!(is_prime_miller_rabin(2305843009213693951, 20));
    }

    #[test]
    fn chinese_remainder_basic() {
        // x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7) → x = 23
        let r = vec![2, 3, 2];
        let m = vec![3, 5, 7];
        assert_eq!(chinese_remainder(&r, &m).unwrap(), 23);
    }

    #[test]
    fn chinese_remainder_non_coprime() {
        let r = vec![1, 2];
        let m = vec![4, 6];
        assert!(chinese_remainder(&r, &m).is_err());
    }
}
