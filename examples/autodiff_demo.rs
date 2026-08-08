//! Example: automatic differentiation with dual numbers.
//!
//! Run with: `cargo run --example autodiff_demo`

use mathr::autodiff;
use mathr::eval::Context;
use mathr::parser::Parser;

fn main() {
    println!("=== Automatic Differentiation Demo ===\n");

    // Single-variable derivative
    // f(x) = x^3 + 2x^2 - x + 5, f'(x) = 3x^2 + 4x - 1
    let expr = Parser::parse("x^3 + 2*x^2 - x + 5").unwrap();
    let ctx = Context::standard();
    let d = autodiff::derivative(&expr, "x", 2.0, &ctx).unwrap();
    println!("f(x) = x³ + 2x² - x + 5");
    println!("f(2) = {},  f'(2) = {}", d.val, d.deriv);
    println!();

    // Trig composition: f(x) = sin(x^2), f'(x) = 2x*cos(x^2)
    let expr = Parser::parse("sin(x^2)").unwrap();
    let d = autodiff::derivative(&expr, "x", 1.5, &ctx).unwrap();
    println!("f(x) = sin(x^2)");
    println!("f(1.5) = {}  f'(1.5) = {}", d.val, d.deriv);
    println!();

    // Exponential: f(x) = exp(-(x^2)), f'(x) = -2x*exp(-x^2)
    let expr = Parser::parse("exp(-(x^2))").unwrap();
    let d = autodiff::derivative(&expr, "x", 1.0, &ctx).unwrap();
    println!("f(x) = exp(-x^2)");
    println!("f(1) = {}  f'(1) = {}", d.val, d.deriv);
    println!();

    // Multivariate gradient: f(x,y) = x^2 + y^3
    let expr = Parser::parse("x^2 + y^3").unwrap();
    let mut point = Context::standard();
    point.set("x", 2.0);
    point.set("y", 3.0);
    let grad = autodiff::gradient(&expr, &point).unwrap();
    println!("f(x,y) = x^2 + y^3  at (2, 3)");
    for (name, val) in &grad {
        println!("  df/d{} = {}", name, val);
    }
    println!();

    // Jacobian: f1 = x^2 + y, f2 = x*y^2
    let f1 = Parser::parse("x^2 + y").unwrap();
    let f2 = Parser::parse("x * y^2").unwrap();
    let mut point = Context::standard();
    point.set("x", 2.0);
    point.set("y", 3.0);
    let jac = autodiff::jacobian(&[f1, f2], &point).unwrap();
    println!("Jacobian of [x^2+y, x*y^2] at (2, 3):");
    for (i, row) in jac.iter().enumerate() {
        print!("  row {}: [", i);
        for (j, v) in row.iter().enumerate() {
            if j > 0 {
                print!(", ");
            }
            print!("{:.4}", v);
        }
        println!("]");
    }
}
