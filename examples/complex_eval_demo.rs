//! Demonstrates complex-number evaluation of expression ASTs via `ceval`.

use mathr::ceval::{eval_complex_str, format_complex};

fn main() {
    println!("=== Complex evaluation (i = imaginary unit) ===\n");

    let cases = [
        ("(1 + 2i) * (3 - i)", "(1+2i)(3-i) = 5 + 5i"),
        ("i^2", "i^2 = -1"),
        ("sqrt(-1)", "sqrt(-1) = i (principal branch)"),
        ("exp(i*pi)", "Euler: e^(i*pi) = -1"),
        ("ln(-1)", "ln(-1) = i*pi (principal branch)"),
        ("sin(2i)", "sin(2i) = i*sinh(2)"),
        ("abs(3 + 4i)", "|3+4i| = 5"),
        ("x*i - y with x=2, y=-1", "variables bind as real values"),
    ];

    for (expr, note) in cases {
        let z = eval_complex_str(expr).unwrap();
        println!("  cval {:<28} = {:<22} ({})", expr, format_complex(z), note);
    }
}
