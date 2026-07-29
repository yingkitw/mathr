use mathr::parser::Parser;
use mathr::eval::{eval, Context};

fn main() {
    println!("=== Special Functions Demo ===\n");

    let ctx = Context::standard();

    // --- Gamma function ---
    println!("--- Gamma function Γ(x) ---");
    for x in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0] {
        let e = Parser::parse(&format!("gamma({})", x)).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  Γ({}) = {:.10}", x, v);
    }

    // --- Error function ---
    println!("\n--- Error function erf(x) ---");
    for x in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0] {
        let e = Parser::parse(&format!("erf({})", x)).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  erf({:5.1}) = {:.10}", x, v);
    }

    // --- Complementary error function ---
    println!("\n--- erfc(x) = 1 - erf(x) ---");
    for x in [0.0, 0.5, 1.0, 2.0] {
        let e = Parser::parse(&format!("erfc({})", x)).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  erfc({}) = {:.10}", x, v);
    }

    // --- sinc function ---
    println!("\n--- sinc function ---");
    for x in [-1.0, -0.5, 0.0, 0.5, 1.0, 3.14159] {
        let e = Parser::parse(&format!("sinc({})", x)).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  sinc({:7.5}) = {:.10}", x, v);
    }

    // --- Using TeX syntax for special functions ---
    println!("\n--- TeX syntax for special functions ---");
    let tex_exprs = [
        (r"\Gamma{0.5}", "Γ(0.5) = √π"),
        (r"\Gamma{10}", "Γ(10) = 9! = 362880"),
        (r"\operatorname{erf}(1.0)", "erf(1)"),
        (r"\operatorname{erfc}(0.0)", "erfc(0) = 1"),
        (r"\operatorname{sinc}(0.0)", "sinc(0) = 1"),
    ];
    for (tex, desc) in &tex_exprs {
        let e = Parser::parse(tex).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  {:30} {} = {:.10}", desc, tex, v);
    }

    // --- Combined expressions ---
    println!("\n--- Combined expressions ---");
    let combined = [
        ("gamma(0.5) * gamma(0.5)", "Γ(0.5)² = π"),
        ("erf(1.0) + erfc(1.0)", "erf(1) + erfc(1) = 1"),
        ("sinc(pi)", "sinc(π) = 0"),
        ("gamma(5) / gamma(3)", "Γ(5)/Γ(3) = 24/2 = 12"),
    ];
    for (expr, desc) in &combined {
        let e = Parser::parse(expr).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  {:30} = {:.10}", desc, v);
    }
}
