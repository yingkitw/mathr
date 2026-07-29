use mathr::parser::Parser;
use mathr::eval::{eval, Context};
use mathr::symbolic;
use mathr::simplify;
use mathr::solver::{newton_central, SolveOptions};

fn main() {
    println!("=== TeX / Markdown Input Demo ===\n");

    // --- Evaluate TeX expressions ---
    let tex_exprs = [
        (r"\frac{1}{2} + \frac{3}{4}", "fraction addition"),
        (r"\sqrt{144}", "square root"),
        (r"\sin(\pi / 4)", "sin of pi/4"),
        (r"\cos(0) + \tan(\pi / 4)", "trig identities"),
        (r"2 \cdot \pi", "2 times pi"),
        (r"\left( 3 + 4 \right) \cdot 2", "grouped multiplication"),
        (r"\log_2{1024}", "log base 2"),
        (r"\Gamma{5}", "Gamma function"),
        (r"\operatorname{erf}(1.0)", "error function"),
        (r"\frac{\frac{1}{2}}{\frac{3}{4}}", "nested fractions"),
        (r"2^{10}", "power with braces"),
        (r"3\sqrt{16} + \frac{1}{2}", "mixed TeX"),
    ];

    let ctx = Context::standard();
    for (tex, desc) in &tex_exprs {
        let e = Parser::parse(tex).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  {:40} {} = {}", desc, tex, v);
    }

    // --- Markdown delimiters ---
    println!("\n--- Markdown delimiters ---");
    let md_exprs = [
        (r"$\sin(\pi / 4)$", "inline $...$"),
        (r"$$\frac{1}{2} + \frac{3}{4}$$", "display $$...$$"),
        (r"\(\sqrt{16}\)", "inline \\(...\\)"),
        (r"\[\cos(0) + \sin(0)\]", "display \\[...\\]"),
    ];
    for (md, desc) in &md_exprs {
        let e = Parser::parse(md).unwrap();
        let v = eval(&e, &ctx).unwrap();
        println!("  {:30} {} = {}", desc, md, v);
    }

    // --- Differentiate TeX input ---
    println!("\n--- Symbolic differentiation of TeX input ---");
    let diff_exprs = [
        r"\sin(x^2)",
        r"\frac{x^2 + 1}{x - 1}",
        r"\exp(x) \cdot \cos(x)",
        r"\sqrt{x}",
    ];
    for tex in &diff_exprs {
        let e = Parser::parse(tex).unwrap();
        let d = symbolic::differentiate(&e, "x").unwrap();
        let s = simplify::simplify(&d);
        println!("  d/dx[{}] = {}", tex, s);
    }

    // --- Solve TeX input ---
    println!("\n--- Solving TeX equations ---");
    let solve_exprs = [
        (r"\frac{x^2 - 4}{1}", "x^2 - 4 = 0"),
        (r"\sqrt{x} - 2", "sqrt(x) = 2"),
        (r"\sin(x)", "sin(x) = 0"),
    ];
    for (tex, desc) in &solve_exprs {
        let e = Parser::parse(tex).unwrap();
        let ctx = Context::standard();
        let expr_clone = e.clone();
        let ctx_clone = ctx.clone();
        let f = move |x: f64| {
            let mut c = ctx_clone.clone();
            c.set("x", x);
            eval(&expr_clone, &c).unwrap_or(f64::NAN)
        };
        match newton_central(f, 1.0, SolveOptions::default()) {
            Ok((root, residual)) => println!("  {:25} {} → root ≈ {:.6} (f = {:.2e})", desc, tex, root, residual),
            Err(e) => println!("  {:25} {} → {}", desc, tex, e),
        }
    }

    // --- TeX matches plain text ---
    println!("\n--- TeX vs plain text equivalence ---");
    let pairs = [
        (r"\frac{x+1}{x-1}", "(x+1)/(x-1)"),
        (r"\sin(x) + \cos(x)", "sin(x) + cos(x)"),
        (r"\sqrt{x^2 + 1}", "sqrt(x^2 + 1)"),
        (r"2 \cdot \pi + \frac{1}{2}", "2*pi + 1/2"),
    ];
    for (tex, plain) in &pairs {
        let te = Parser::parse(tex).unwrap();
        let pe = Parser::parse(plain).unwrap();
        let eq = te.canonicalize() == pe.canonicalize();
        println!("  {} ≡ {} : {}", tex, plain, if eq { "✓" } else { "✗" });
    }
}
