use mathr::apart;
use mathr::parser::Parser;

fn show(src: &str) {
    let e = Parser::parse(src).unwrap();
    let result = apart::apart(&e, "x").unwrap();
    println!("  {}  =  {}", src, result);
}

fn show_steps(src: &str) {
    let e = Parser::parse(src).unwrap();
    let steps = apart::apart_steps(&e, "x").unwrap();
    println!("  Steps for apart({}):", src);
    for s in &steps {
        println!("    • {}", s);
    }
}

fn main() {
    println!("=== Partial Fraction Decomposition Demo ===\n");

    println!("1. Simple linear factors:");
    show("1/(x*(x + 1))");
    show("1/(x^2 - 1)");
    show("5/(x^2 - 9)");
    println!();

    println!("2. Repeated linear factors:");
    show("1/(x*(x - 1)^2)");
    show("1/(x^2*(x - 1)^2)");
    println!();

    println!("3. Irreducible quadratic factors:");
    show("1/(x^2 + 1)");
    show("x/((x - 1)*(x^2 + 1))");
    println!();

    println!("4. Improper fractions (long division first):");
    show("(x^3 + 1)/(x + 1)");
    show("(3*x^2 + 2)/(2*x - 1)");
    println!();

    println!("5. Higher degree denominators:");
    show("1/(x^3 - 1)");
    show("1/(x^4 - 1)");
    println!();

    println!("6. Step-by-step example:");
    show_steps("1/(x^2 - 1)");
    println!();

    // Verify numeric equivalence
    println!("7. Numeric verification:");
    let e = Parser::parse("(x^2 + 2*x + 3)/(x^3 - x)").unwrap();
    let result = apart::apart(&e, "x").unwrap();
    let mut ctx = mathr::eval::Context::standard();
    for &x in &[2.0, 5.0, -3.0, 0.5] {
        ctx.set("x", x);
        let a = mathr::eval::eval(&e, &ctx).unwrap();
        let b = mathr::eval::eval(&result, &ctx).unwrap();
        println!("  x = {:>5}: original = {:.6}, apart = {:.6}, match = {}", x, a, b, (a - b).abs() < 1e-9);
    }
}
