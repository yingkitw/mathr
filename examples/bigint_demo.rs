//! Example: arbitrary-precision big integer arithmetic.
//!
//! Run with: `cargo run --example bigint_demo`

use mathr::bigint;
use num_traits::One;

fn main() {
    println!("=== Big Integer Demo ===\n");

    // Factorial — 50! is way beyond u64
    let f50 = bigint::factorial(50);
    println!("50! = {}", f50);
    println!();

    // Fibonacci — F_100 is way beyond u64
    let f100 = bigint::fibonacci(100);
    println!("F_100 = {}", f100);
    println!();

    // Binomial coefficient — C(100, 50)
    let c = bigint::binomial(100, 50);
    println!("C(100, 50) = {}", c);
    println!();

    // Primality testing on large numbers
    let mersenne_61 = num_bigint::BigInt::from(1u64 << 61) - num_bigint::BigInt::one();
    println!("2^61 - 1 = {}", mersenne_61);
    println!("Is prime? {}", bigint::is_prime(&mersenne_61, 20));
    println!();

    // Factorization of a large semiprime
    let semiprime = num_bigint::BigInt::from(1000000007u64)
        * num_bigint::BigInt::from(1000000009u64);
    println!("Factorizing {} ...", semiprime);
    let factors = bigint::factorize(&semiprime);
    let parts: Vec<String> = factors
        .iter()
        .map(|(p, e)| {
            if *e == 1 {
                p.to_string()
            } else {
                format!("{}^{}", p, e)
            }
        })
        .collect();
    println!("  = {}", parts.join(" * "));
    println!();

    // Modular exponentiation with big numbers
    let base = num_bigint::BigInt::from(2);
    let exp = num_bigint::BigInt::from(1000);
    let modulus = num_bigint::BigInt::from(1000000007u64);
    let result = bigint::mod_pow(&base, &exp, &modulus).unwrap();
    println!("2^1000 mod 1000000007 = {}", result);
    println!();

    // Euler's totient
    let n = num_bigint::BigInt::from(360);
    println!("φ({}) = {}", n, bigint::totient(&n));
}
