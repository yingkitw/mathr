use mathr::numtheory;

fn main() {
    // --- Miller–Rabin primality test ---
    println!("Miller–Rabin primality test:");
    let candidates = [
        2u64, 3, 17, 97, 1009, 7919, 104729, 1299709,
        2305843009213693951, // Mersenne prime M_61
        561,                 // Carmichael number (composite)
        1729,                // Carmichael number
    ];
    for n in candidates {
        let is_prime = numtheory::is_prime_miller_rabin(n, 20);
        println!("  {:>20} → {}", n, if is_prime { "prime" } else { "composite" });
    }

    // --- Sieve of Eratosthenes ---
    let primes = numtheory::sieve_primes(100);
    println!("\nPrimes below 100 ({} total):", primes.len());
    println!("  {:?}", primes);

    // --- Prime factorization ---
    println!("\nPrime factorization:");
    for n in [60u64, 360, 1024, 999983, 1234567890] {
        let factors = numtheory::prime_factors(n);
        let s: String = factors.iter().map(|f| f.to_string()).collect::<Vec<_>>().join(" × ");
        println!("  {} = {}", n, s);
    }

    // --- Modular exponentiation ---
    println!("\nModular exponentiation:");
    let bases = [2u64, 3, 7];
    let exps = [10u64, 20, 100];
    let mods = [1000u64, 10000, 100000];
    for i in 0..3 {
        let r = numtheory::mod_pow(bases[i], exps[i], mods[i]);
        println!("  {}^{} mod {} = {}", bases[i], exps[i], mods[i], r);
    }

    // --- Chinese Remainder Theorem ---
    println!("\nChinese Remainder Theorem:");
    let remainders = [2u64, 3, 2];
    let moduli = [3u64, 5, 7];
    match numtheory::chinese_remainder(&remainders, &moduli) {
        Ok(x) => println!("  x ≡ {} (mod 3), x ≡ {} (mod 5), x ≡ {} (mod 7) → x = {}", remainders[0], remainders[1], remainders[2], x),
        Err(e) => println!("  error: {}", e),
    }

    // --- GCD, LCM, binomial, factorial, Fibonacci ---
    println!("\nMiscellaneous:");
    println!("  gcd(48, 36) = {}", numtheory::gcd(48, 36));
    println!("  lcm(4, 6)   = {}", numtheory::lcm(4, 6));
    println!("  C(10, 3)    = {}", numtheory::binomial(10, 3).unwrap());
    println!("  20!         = {}", numtheory::factorial(20).unwrap());
    println!("  fib(50)     = {}", numtheory::fibonacci(50));
    println!("  φ(36)       = {}", numtheory::euler_totient(36));
}
