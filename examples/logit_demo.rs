//! Logistic regression (IRLS) demo.
//!
//! Run with: cargo run --example logit_demo

use mathr::logit::{logistic_regression, predict_proba, LogitOptions};

fn main() {
    // does study hours predict passing? balanced two-point design:
    // 1/3 pass rate at 1 hour, 2/3 at 2 hours
    let x = [&[1.0, 1.0, 1.0, 2.0, 2.0, 2.0][..]];
    let y = [0.0, 0.0, 1.0, 0.0, 1.0, 1.0];
    let fit = logistic_regression(&x, &y, &LogitOptions::default()).unwrap();

    println!("logistic regression: pass ~ hours");
    println!(
        "  intercept = {:+.6} ± {:.6}",
        fit.coefficients[0], fit.std_errors[0]
    );
    println!(
        "  hours     = {:+.6} ± {:.6}",
        fit.coefficients[1], fit.std_errors[1]
    );
    println!("  log-lik   = {:.6} ({} iterations, {})",
        fit.log_likelihood, fit.iterations,
        if fit.converged { "converged" } else { "NOT converged" });

    for hours in [1.0, 2.0] {
        let p = predict_proba(&fit.coefficients, &[hours]).unwrap();
        println!("  P(pass | {hours} h) = {p:.4}");
    }

    // multivariate: hours + sleep (duplicated design points, provably
    // non-separable; x1 carries no effect after conditioning on x2)
    let x2 = [
        &[0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0][..],
        &[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0][..],
    ];
    let y2 = [0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0];
    let fit2 = logistic_regression(&x2, &y2, &LogitOptions::default()).unwrap();
    println!("\nmultivariate: outcome ~ noise + group");
    for (j, (b, se)) in fit2
        .coefficients
        .iter()
        .zip(fit2.std_errors.iter())
        .enumerate()
    {
        println!("  b{j} = {b:+.6} ± {se:.6}");
    }
}
