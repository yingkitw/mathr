use mathr::limit::{self, LimitValue};
use mathr::parser::Parser;

fn show(src: &str, point: f64) {
    let e = Parser::parse(src).unwrap();
    let v = limit::limit(&e, "x", point).unwrap();
    let pt = limit::fmt_point(point);
    println!("  lim x→{}  {}  =  {}", pt, src, v);
}

fn show_steps(src: &str, point: f64) {
    let e = Parser::parse(src).unwrap();
    let steps = limit::limit_steps(&e, "x", point).unwrap();
    println!("  Steps for lim x→{}  {}:", limit::fmt_point(point), src);
    for s in &steps {
        println!("    • {}", s);
    }
}

fn main() {
    println!("=== Limits Demo ===\n");

    println!("1. Removable singularities (L'Hôpital):");
    show("(x^2 - 1)/(x - 1)", 1.0);
    show("sin(x)/x", 0.0);
    show("(1 - cos(x))/x^2", 0.0);
    show("(exp(x) - 1)/x", 0.0);
    show("ln(x + 1)/x", 0.0);
    println!();

    println!("2. Limits at infinity:");
    show("1/x", f64::INFINITY);
    show("exp(x)/x^2", f64::INFINITY);
    show("ln(x)/x", f64::INFINITY);
    show("x/sqrt(x^2 + 1)", f64::INFINITY);
    show("exp(x)", f64::NEG_INFINITY);
    println!();

    println!("3. Poles and divergences:");
    show("1/x^2", 0.0);
    show("1/(x - 3)", 3.0);
    println!();

    println!("4. Does not exist:");
    show("sin(1/x)", 0.0);
    show("exp(1/x)", 0.0);
    println!();

    println!("5. Step-by-step example:");
    show_steps("(x^2 - 1)/(x - 1)", 1.0);
    println!();

    println!("6. Deep L'Hôpital (3 applications):");
    let e = Parser::parse("(x - sin(x))/x^3").unwrap();
    let v = limit::limit(&e, "x", 0.0).unwrap();
    println!("  lim x→0  (x - sin(x))/x^3  =  {}  (exact: 1/6 ≈ 0.16667)", v);
    println!();

    // Demonstrate the LimitValue enum
    println!("7. Programmatic usage:");
    let e = Parser::parse("1/x").unwrap();
    let v = limit::limit(&e, "x", 0.0).unwrap();
    match v {
        LimitValue::Finite(x) => println!("  1/x at 0 → finite: {}", x),
        LimitValue::PosInfinity => println!("  1/x at 0 → +∞"),
        LimitValue::NegInfinity => println!("  1/x at 0 → -∞"),
        LimitValue::DoesNotExist => println!("  1/x at 0 → does not exist"),
    }
}
