use mathr::numtheory;

fn main() {
    println!("=== Number Theory Demo ===\n");

    // --- Primality ---
    println!("--- Primality tests ---");
    let primes_to_check = [2, 7, 17, 97, 100, 561, 7919, 104729];
    for n in &primes_to_check {
        let trial = numtheory::is_prime(*n);
        let mr = numtheory::is_prime_miller_rabin(*n, 20);
        println!("  {:6}  trial: {}  Miller-Rabin: {}", n, trial, mr);
    }

    // 561 is a Carmichael number (pseudoprime)
    println!("  Note: 561 is a Carmichael number — trial division catches it");

    // --- Factorization ---
    println!("\n--- Prime factorization ---");
    let nums = [60, 360, 1024, 999983, 1234567890];
    for n in &nums {
        let factors = numtheory::prime_factors(*n);
        let product: u64 = factors.iter().product();
        let factor_str: Vec<String> = factors.iter().map(|f| f.to_string()).collect();
        println!("  {:10} = {}  (product check: {})", n, factor_str.join(" × "), product == *n);
    }

    // --- GCD and LCM ---
    println!("\n--- GCD and LCM ---");
    let pairs = [(48, 36), (17, 5), (100, 75), (1024, 768)];
    for (a, b) in &pairs {
        let g = numtheory::gcd(*a, *b);
        let l = numtheory::lcm(*a, *b);
        println!("  gcd({}, {}) = {}, lcm({}, {}) = {}", a, b, g, a, b, l);
    }

    // --- Fibonacci ---
    println!("\n--- Fibonacci numbers ---");
    for n in [0, 1, 5, 10, 20, 50, 90] {
        println!("  fib({:2}) = {}", n, numtheory::fibonacci(n));
    }

    // --- Binomial coefficients ---
    println!("\n--- Binomial coefficients C(n,k) ---");
    for n in [5, 10, 20] {
        let row: Vec<String> = (0..=n)
            .map(|k| numtheory::binomial(n, k).unwrap().to_string())
            .collect();
        println!("  n={:2}: {}", n, row.join(" "));
    }

    // --- Factorials ---
    println!("\n--- Factorials ---");
    for n in [0, 1, 5, 10, 20] {
        println!("  {}! = {}", n, numtheory::factorial(n).unwrap());
    }

    // --- Sieve of Eratosthenes ---
    println!("\n--- Sieve of Eratosthenes ---");
    let primes = numtheory::sieve_primes(100);
    println!("  Primes below 100 ({} total):", primes.len());
    let prime_str: Vec<String> = primes.iter().map(|p| p.to_string()).collect();
    println!("  {}", prime_str.join(", "));

    // --- Euler's totient ---
    println!("\n--- Euler's totient φ(n) ---");
    for n in [1, 6, 9, 12, 36, 100] {
        let phi = numtheory::euler_totient(n);
        println!("  φ({:3}) = {}", n, phi);
    }

    // --- Chinese Remainder Theorem ---
    println!("\n--- Chinese Remainder Theorem ---");
    let remainders = [2u64, 3, 2];
    let moduli = [3u64, 5, 7];
    let x = numtheory::chinese_remainder(&remainders, &moduli).unwrap();
    println!("  x ≡ {} (mod {}), {} (mod {}), {} (mod {}) → x = {}", 
        remainders[0], moduli[0], remainders[1], moduli[1], remainders[2], moduli[2], x);
    for i in 0..3 {
        println!("  {} mod {} = {} ✓", x, moduli[i], x % moduli[i]);
    }

    // --- Modular exponentiation ---
    println!("\n--- Modular exponentiation ---");
    let result = numtheory::mod_pow(2, 10, 1000);
    println!("  2^10 mod 1000 = {}", result);
    let result = numtheory::mod_pow(3, 100, 7);
    println!("  3^100 mod 7 = {}", result);
}
