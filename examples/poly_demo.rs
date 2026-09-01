//! Example: polynomial expansion (distribution) over expressions.
//!
//! Run with: `cargo run --example poly_demo`

use mathr::parser::Parser;
use mathr::poly::expand;

fn main() {
    println!("=== Polynomial Expansion Demo ===\n");

    for src in [
        "(x+1)^3",
        "(x+2)*(x+3)",
        "(x+y)^2",
        "(x+y)*(x-y)",
        "(x - sin(x))/x^3",
    ] {
        let e = Parser::parse(src).unwrap();
        println!("expand({}) = {}", src, expand(&e));
    }
}
