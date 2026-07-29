use mathr::parser::Parser;
use mathr::eval::{eval, Context};
use mathr::symbolic;
use mathr::simplify;
use mathr::taylor;
use mathr::solver::{newton_central, bisect, SolveOptions};

fn main() {
    println!("=== Calculus Demo: diff, integrate, solve, taylor ===\n");

    // --- Symbolic differentiation ---
    println!("--- Differentiation ---");
    let exprs = [
        "x^5",
        "sin(x^2)",
        "exp(x) * cos(x)",
        "ln(x) / x",
        "x^3 - 2*x^2 + x - 7",
        "(x^2 + 1) / (x - 1)",
    ];
    for e in &exprs {
        let expr = Parser::parse(e).unwrap();
        let d = symbolic::differentiate(&expr, "x").unwrap();
        let s = simplify::simplify(&d);
        println!("  d/dx[{}] = {}", e, s);
    }

    // --- Higher-order derivatives ---
    println!("\n--- Higher-order derivatives of x^6 ---");
    let f = Parser::parse("x^6").unwrap();
    let mut d = f.clone();
    for n in 1..=4 {
        d = symbolic::differentiate(&d, "x").unwrap();
        let s = simplify::simplify(&d);
        println!("  d^{}(x^6)/dx^{} = {}", n, n, s);
    }

    // --- Numerical integration ---
    println!("\n--- Numerical integration ---");
    let integrals = [
        ("sin(x)", 0.0, std::f64::consts::PI, "∫₀^π sin(x) dx = 2"),
        ("exp(x)", 0.0, 1.0, "∫₀¹ exp(x) dx = e-1"),
        ("x^2", 0.0, 3.0, "∫₀³ x² dx = 9"),
        ("1/(1+x^2)", 0.0, 1.0, "∫₀¹ 1/(1+x²) dx = π/4"),
    ];
    let ctx = Context::standard();
    for (expr, a, b, desc) in &integrals {
        let e = Parser::parse(expr).unwrap();
        let ctx_clone = ctx.clone();
        let expr_clone = e.clone();
        let f = move |x: f64| {
            let mut c = ctx_clone.clone();
            c.set("x", x);
            eval(&expr_clone, &c).unwrap_or(f64::NAN)
        };
        let result = mathr::calculus::integrate_adaptive(f, *a, *b, 1e-10, 30).unwrap();
        println!("  {} = {:.10}", desc, result);
    }

    // --- Root finding ---
    println!("\n--- Root finding ---");
    let equations = [
        ("x^2 - 4", 1.0, "x² = 4"),
        ("x^3 - x - 2", 1.0, "x³ - x - 2 = 0"),
        ("exp(x) - 2", 0.0, "eˣ = 2"),
        ("sin(x) - 0.5", 1.0, "sin(x) = 0.5"),
    ];
    for (expr, guess, desc) in &equations {
        let e = Parser::parse(expr).unwrap();
        let ctx = Context::standard();
        let expr_clone = e.clone();
        let ctx_clone = ctx.clone();
        let f = move |x: f64| {
            let mut c = ctx_clone.clone();
            c.set("x", x);
            eval(&expr_clone, &c).unwrap_or(f64::NAN)
        };
        let (root, residual) = newton_central(f, *guess, SolveOptions::default()).unwrap();
        println!("  {:20} → x ≈ {:.10} (f = {:.2e})", desc, root, residual);
    }

    // --- Bisection ---
    println!("\n--- Bisection ---");
    let e = Parser::parse("x^3 - x - 2").unwrap();
    let ctx = Context::standard();
    let expr_clone = e.clone();
    let ctx_clone = ctx.clone();
    let f = move |x: f64| {
        let mut c = ctx_clone.clone();
        c.set("x", x);
        eval(&expr_clone, &c).unwrap_or(f64::NAN)
    };
    let (root, residual) = bisect(f, 1.0, 2.0, SolveOptions::default()).unwrap();
    println!("  x³ - x - 2 on [1,2] → x ≈ {:.10} (f = {:.2e})", root, residual);

    // --- Taylor series ---
    println!("\n--- Taylor series ---");
    let taylor_exprs = [
        ("exp(x)", 0.0, 6),
        ("sin(x)", 0.0, 7),
        ("cos(x)", 0.0, 6),
        ("ln(x+1)", 0.0, 5),
    ];
    for (expr, a, order) in &taylor_exprs {
        let series = taylor::taylor_series_str(expr, "x", *a, *order).unwrap();
        println!("  Taylor[{}, a={}] ({} terms) = {}", expr, a, order, series);
    }

    // --- Partial derivatives ---
    println!("\n--- Partial derivatives ---");
    let f = Parser::parse("x^2 * y + y^3 * x").unwrap();
    let df_dx = symbolic::differentiate(&f, "x").unwrap();
    let df_dy = symbolic::differentiate(&f, "y").unwrap();
    println!("  f(x,y) = x²y + y³x");
    println!("  ∂f/∂x = {}", simplify::simplify(&df_dx));
    println!("  ∂f/∂y = {}", simplify::simplify(&df_dy));
}
