//! Example: arbitrary-precision decimal arithmetic.
//!
//! Run with: `cargo run --example bigdec_demo`

use mathr::bigdec::{self, BigDecimal};
use std::collections::HashMap;

fn main() {
    println!("=== Arbitrary-Precision Decimal Demo ===\n");

    // Constants at arbitrary precision
    println!("π  (100 digits) = {}", bigdec::pi(100).unwrap());
    println!("e  (50 digits)  = {}", bigdec::e(50).unwrap());
    println!();

    // Square roots, exponentials, logarithms
    let two = BigDecimal::from(2);
    println!("√2 (40 digits)  = {}", bigdec::sqrt(&two, 40).unwrap());
    println!("ln 2 (40)       = {}", bigdec::ln(&two, 40).unwrap());
    println!();

    // Exact division: 1/3 to 40 significant digits
    let third = bigdec::div(&BigDecimal::from(1), &BigDecimal::from(3), 40).unwrap();
    println!("1/3 (40 digits) = {}", third);
    println!();

    // Trigonometry with argument reduction across many periods
    let one = BigDecimal::from(1);
    println!("sin(1) (30)     = {}", bigdec::sin(&one, 30).unwrap());
    println!("cos(1) (30)     = {}", bigdec::cos(&one, 30).unwrap());
    println!("tan(1) (30)     = {}", bigdec::tan(&one, 30).unwrap());
    println!();

    // Evaluate a whole expression tree at high precision
    let expr = mathr::parser::Parser::parse("sin(10*pi + pi/6)").unwrap();
    let vars: HashMap<String, BigDecimal> = HashMap::new();
    let result = bigdec::eval_decimal_rounded(&expr, &vars, 40).unwrap();
    println!("sin(10π + π/6)  = {}   (exactly 1/2)", result);
    println!();

    // Power with a non-integer exponent: 2^0.5
    let half = bigdec::parse("0.5").unwrap();
    println!("2^0.5 (40)      = {}", bigdec::pow(&two, &half, 40).unwrap());
}
