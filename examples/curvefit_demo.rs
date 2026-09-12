//! Levenberg–Marquardt nonlinear least-squares curve fitting.
//!
//! Run with: cargo run --example curvefit_demo

use mathr::curvefit::{curve_fit, LmFit, LmOptions};

fn show(label: &str, fit: &LmFit, truth: &[f64]) {
    println!("  {label}");
    for (i, (p, t)) in fit.params.iter().zip(truth.iter()).enumerate() {
        println!(
            "    p{} = {:<14.8}  (true {t})  ± {:.2e}",
            i,
            p,
            fit.std_errors.get(i).copied().unwrap_or(f64::NAN)
        );
    }
    println!(
        "    sse = {:.3e}, iterations = {} ({})",
        fit.sse,
        fit.iterations,
        if fit.converged { "converged" } else { "NOT converged" }
    );
}

fn main() {
    // Exponential decay: y = 3 e^{-0.7 x}, sampled at x = 0..9
    let decay: Vec<(f64, f64)> = (0..10)
        .map(|i| {
            let x = i as f64;
            (x, 3.0 * (-0.7 * x).exp())
        })
        .collect();
    let fit = curve_fit(
        |x, p| Ok(p[0] * (-p[1] * x).exp()),
        &decay,
        &[1.0, 1.0],
        &LmOptions::default(),
    )
    .unwrap();
    show("decay model a·exp(-b·x) from init [1, 1]:", &fit, &[3.0, 0.7]);

    // Michaelis–Menten: y = v·x / (k + x), true v = 4, k = 1
    let mm: Vec<(f64, f64)> = (1..=20)
        .map(|i| {
            let x = 0.25 * i as f64;
            (x, 4.0 * x / (1.0 + x))
        })
        .collect();
    let fit = curve_fit(
        |x, p| Ok(p[0] * x / (p[1] + x)),
        &mm,
        &[1.0, 1.0],
        &LmOptions::default(),
    )
    .unwrap();
    show("michaelis-menten v·x/(k+x) from init [1, 1]:", &fit, &[4.0, 1.0]);

    // Linear models reproduce closed-form regression
    let lin = [(0.0, 1.0), (1.0, 2.1), (2.0, 2.9), (3.0, 4.2)];
    let (slope, intercept) =
        mathr::stats::linear_regression(&[0.0, 1.0, 2.0, 3.0], &[1.0, 2.1, 2.9, 4.2]).unwrap();
    let fit = curve_fit(|x, p| Ok(p[0] * x + p[1]), &lin, &[0.5, 0.5], &LmOptions::default())
        .unwrap();
    show(
        "linear model a·x + b (cross-checked against stats::linear_regression):",
        &fit,
        &[slope, intercept],
    );
}
