use mathr::special;

fn main() {
    // --- Gamma function ---
    println!("Gamma function:");
    for x in [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 0.1] {
        println!("  Γ({:5.1}) = {:.10}", x, special::gamma(x));
    }
    // Γ(1/2) = √π
    let sqrt_pi = std::f64::consts::PI.sqrt();
    println!("  Γ(0.5) = √π = {:.10}  (error: {:.2e})", special::gamma(0.5), (special::gamma(0.5) - sqrt_pi).abs());

    // --- log-Gamma ---
    println!("\nlog-Gamma:");
    for x in [1.0, 2.0, 10.0, 100.0] {
        println!("  ln Γ({:5.1}) = {:.10}", x, special::log_gamma(x));
    }

    // --- Beta function ---
    println!("\nBeta function:");
    // B(1,1) = 1
    println!("  B(1, 1)   = {:.10}  (expected 1.0)", special::beta(1.0, 1.0));
    // B(0.5, 0.5) = π
    println!("  B(0.5, 0.5) = {:.10}  (expected π)", special::beta(0.5, 0.5));
    // B(2, 3) = 1/12
    println!("  B(2, 3)   = {:.10}  (expected 0.0833...)", special::beta(2.0, 3.0));

    // --- Error function ---
    println!("\nError function:");
    for x in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0] {
        println!("  erf({:5.1}) = {:12.10}   erfc({:5.1}) = {:12.10}",
            x, special::erf(x), x, special::erfc(x));
    }

    // --- Sinc function ---
    println!("\nSinc function:");
    for x in [0.0, 0.5, 1.0, 2.0, 3.0] {
        println!("  sinc({:5.1}) = {:.10}", x, special::sinc(x));
    }

    // --- Incomplete gamma P (chi-squared CDF) ---
    // P(2, x) = 1 - e^{-x/2} is the CDF of χ² with 2 degrees of freedom
    println!("\nIncomplete gamma P(2, x) = χ² CDF (2 dof):");
    for x in [1.0, 2.0, 4.0, 6.0, 10.0] {
        let p = special::incomplete_gamma_p(1.0, x / 2.0);
        println!("  P(1, {:4.1}) = {:.10}", x / 2.0, p);
    }
}
