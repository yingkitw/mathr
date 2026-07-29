use mathr::parser::Parser;
use mathr::eval::{eval, Context};

fn main() {
    let e = Parser::parse("x^3 - 2*x - 5").unwrap();
    eprintln!("parsed: {:#?}", e);
    let mut ctx = Context::standard();
    ctx.set("x", 2.0);
    let f = |x: f64| {
        let mut c = ctx.clone();
        c.set("x", x);
        eval(&e, &c).unwrap_or(f64::NAN)
    };
    println!("f(2) = {}", f(2.0));
    println!("f(2.1) = {}", f(2.1));

    let (r, fv) = mathr::solver::newton_central(f, 2.0, mathr::solver::SolveOptions::default()).unwrap();
    println!("root = {}, f(root) = {}", r, fv);
}
