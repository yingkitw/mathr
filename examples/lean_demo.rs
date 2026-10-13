//! Lean 4 export — render mathr expressions as Lean/Mathlib terms.
//!
//! Run with: cargo run --example lean_demo
//!
//! Every output below compiles in a Lean 4 project with Mathlib
//! (`lake new myproj math` then paste after `import Mathlib`).

use mathr::lean::{to_lean, to_lean_doc};
use mathr::parser::Parser;

fn main() {
    println!("=== basic terms ===");
    for src in ["2", "0.5", "x", "2*x + 3", "x^2 - y^2", "x/y + 1"] {
        let e = Parser::parse(src).unwrap();
        println!("{src:12} → {}", to_lean(&e).unwrap());
    }

    println!("\n=== functions ===");
    for src in ["sin(x) + cos(y)", "ln(x)/log(x, 10)", "sqrt(x^2 + 1)", "floor(x/2)"] {
        let e = Parser::parse(src).unwrap();
        println!("{src:18} → {}", to_lean(&e).unwrap());
    }

    println!("\n=== powers: ℕ vs ℝ exponents ===");
    for src in ["x^2", "x^(1/2)", "(-2)^3"] {
        let e = Parser::parse(src).unwrap();
        println!("{src:10} → {}", to_lean(&e).unwrap());
    }

    println!("\n=== pi is exported as its exact decimal (sound) ===");
    let e = Parser::parse("sin(pi/4)").unwrap();
    println!("sin(pi/4) → {}", to_lean(&e).unwrap());

    println!("\n=== self-contained files (paste into Lean + Mathlib) ===");
    let e = Parser::parse("x^2 + sin(y)").unwrap();
    print!("{}", to_lean_doc(&e).unwrap());

    println!("\n=== proof workflow: differentiate, then export ===");
    let e = Parser::parse("x^3 + 2*x").unwrap();
    let d = mathr::symbolic::differentiate(&e, "x").unwrap();
    println!("d/dx (x^3 + 2*x) = {d}");
    print!("{}", to_lean_doc(&d).unwrap());
}
