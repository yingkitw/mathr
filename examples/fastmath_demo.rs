//! Demonstrates Chebyshev-based fast approximations of transcendental
//! functions and compares them against the standard library.

use mathr::fastmath::{fast_cos, fast_exp, fast_log, fast_sin, fast_sqrt, fast_tan, ChebyshevApprox};

fn main() {
    println!("=== ChebyshevApprox: custom function ===");
    let approx = ChebyshevApprox::new(|x| x.sin() * x.cos(), -3.0, 3.0, 12);
    for x in [-3.0_f64, -1.5, 0.0, 0.7, 2.0, 3.0] {
        let exact = x.sin() * x.cos();
        let fast = approx.eval(x);
        println!("  x={:6.2}  approx={:.15}  exact={:.15}  err={:.2e}", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_sin vs std sin ===");
    for x in [0.0, 0.5, 1.0, 2.0, 3.14159, 10.0, 100.0, -7.0] {
        let fast = fast_sin(x);
        let exact = x.sin();
        println!("  sin({:8.5}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_cos vs std cos ===");
    for x in [0.0, 0.5, 1.0, 2.0, 3.14159, 10.0, 100.0, -7.0] {
        let fast = fast_cos(x);
        let exact = x.cos();
        println!("  cos({:8.5}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_tan vs std tan ===");
    for x in [0.0, 0.3, 0.7, 1.0, 1.2, -0.5] {
        let fast = fast_tan(x);
        let exact = x.tan();
        println!("  tan({:8.5}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_exp vs std exp ===");
    for x in [-5.0, -1.0, 0.0, 0.5, 2.0, 10.0, 50.0] {
        let fast = fast_exp(x);
        let exact = x.exp();
        println!("  exp({:8.5}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_log vs std ln ===");
    for x in [0.001, 0.5, 1.0, 2.0, 10.0, 100.0, 1e10] {
        let fast = fast_log(x);
        let exact = x.ln();
        println!("  ln({:12.5}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== fast_sqrt vs std sqrt ===");
    for x in [0.5, 1.0, 2.0, 10.0, 100.0, 1e6] {
        let fast = fast_sqrt(x);
        let exact = x.sqrt();
        println!("  sqrt({:10.3}) = {:.15}  (exact: {:.15}, err: {:.2e})", x, fast, exact, (fast - exact).abs());
    }

    println!("\n=== Max error summary ===");
    let sin_max = (0..=100_000)
        .map(|i| (fast_sin(i as f64 * 0.001) - (i as f64 * 0.001).sin()).abs())
        .fold(0.0_f64, f64::max);
    let cos_max = (0..=100_000)
        .map(|i| (fast_cos(i as f64 * 0.001) - (i as f64 * 0.001).cos()).abs())
        .fold(0.0_f64, f64::max);
    let exp_max = (0..=10_000)
        .map(|i| {
            let x = (i as f64 - 5000.0) * 0.01;
            (fast_exp(x) - x.exp()).abs() / x.exp().abs()
        })
        .fold(0.0_f64, f64::max);
    let log_max = (1..=100_000)
        .map(|i| {
            let x = i as f64 * 0.001;
            (fast_log(x) - x.ln()).abs()
        })
        .fold(0.0_f64, f64::max);
    println!("  sin  max abs err: {:.2e}", sin_max);
    println!("  cos  max abs err: {:.2e}", cos_max);
    println!("  exp  max rel err: {:.2e}", exp_max);
    println!("  log  max abs err: {:.2e}", log_max);
}
