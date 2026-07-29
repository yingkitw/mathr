use mathr::parser::Parser;
use mathr::symbolic;
use mathr::simplify;
use mathr::taylor;
use mathr::expr::Expr;

fn main() {
    // --- Symbolic differentiation ---
    let exprs = ["x^3 + 2*x^2 - x + 5", "sin(x^2)", "exp(x) * cos(x)", "ln(x) / x"];
    for e in exprs {
        let expr = Parser::parse(e).unwrap();
        let deriv = symbolic::differentiate(&expr, "x").unwrap();
        let simplified = simplify::simplify(&deriv);
        println!("d/dx[{}] = {}", e, simplified);
    }

    // --- Higher-order derivative ---
    let f = Parser::parse("x^4 - 3*x^2 + 2*x - 7").unwrap();
    let d1 = symbolic::differentiate(&f, "x").unwrap();
    let d2 = symbolic::differentiate(&d1, "x").unwrap();
    let d3 = symbolic::differentiate(&d2, "x").unwrap();
    println!("\nHigher-order derivatives of x^4 - 3x^2 + 2x - 7:");
    println!("  f'(x)  = {}", simplify::simplify(&d1));
    println!("  f''(x) = {}", simplify::simplify(&d2));
    println!("  f'''(x)= {}", simplify::simplify(&d3));

    // --- Taylor series ---
    let series = taylor::taylor_series_str("exp(x)", "x", 0.0, 6).unwrap();
    println!("\nTaylor series of exp(x) around 0 (6 terms):");
    println!("  {}", series);

    let series = taylor::taylor_series_str("sin(x)", "x", 0.0, 7).unwrap();
    println!("\nTaylor series of sin(x) around 0 (7 terms):");
    println!("  {}", series);

    // --- Expression equality ---
    let a = Expr::add(Expr::var("x"), Expr::var("y"));
    let b = Expr::add(Expr::var("y"), Expr::var("x"));
    println!("\nx + y == y + x: {}", a.equals(&b));
}
