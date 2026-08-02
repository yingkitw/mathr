//! An interactive REPL that ties together the parser, evaluator, symbolic
//! differentiation, numeric integration, root-finders and plotter.
//!
//! Commands available at the prompt:
//!   `\<expr\>`             evaluate
//!   `let x = 2`            bind a variable
//!   `fn f(x) = x^2`        define a function
//!   `diff \<expr\> \[wrt\]`  symbolic derivative
//!   `int \<expr\> a b`     numerical integral
//!   `solve \<expr\> \[wrt\] \[guess\]`
//!   `simplify \<expr\>`
//!   `plot \<expr\> a b \[out.png\]`
//!   `vars | funcs | clear | help | quit`

use crate::error::Result;
use crate::eval::{eval, Context, Func};
use crate::expr::Expr;
use crate::parser::Parser;
use crate::simplify::simplify;
use crate::symbolic::differentiate;
use rustyline::completion::FilenameCompleter;
use rustyline::error::ReadlineError;
use rustyline::highlight::MatchingBracketHighlighter;
use rustyline::hint::HistoryHinter;
use rustyline::validate::MatchingBracketValidator;
use rustyline::{Completer, Editor, Helper, Highlighter, Hinter, Validator};
use std::borrow::Cow;
use std::collections::HashSet;

#[derive(Helper, Completer, Highlighter, Hinter, Validator)]
struct ReplHelper {
    #[rustyline(Completer)]
    completer: FilenameCompleter,
    _highlighter: MatchingBracketHighlighter,
    #[rustyline(Validator)]
    validator: MatchingBracketValidator,
    #[rustyline(Hinter)]
    hinter: HistoryHinter,
}

impl Default for ReplHelper {
    fn default() -> Self {
        Self {
            completer: FilenameCompleter::new(),
            _highlighter: MatchingBracketHighlighter::new(),
            validator: MatchingBracketValidator::new(),
            hinter: HistoryHinter::new(),
        }
    }
}

pub fn run() -> Result<()> {
    let helper = ReplHelper::default();
    let mut rl = Editor::new().map_err(|e| crate::error::MathError::Other(format!("repl: {}", e)))?;
    rl.set_helper(Some(helper));
    let _ = rl.load_history(".mathr_history");
    println!("mathr {} — type `help` for a list of commands.", env!("CARGO_PKG_VERSION"));

    let mut ctx = Context::standard();
    let mut stdout = std::io::stdout();

    loop {
        let prompt = if ctx.vars.is_empty() {
            "\nmathr> "
        } else {
            "\nmathr* "
        };
        let line = match rl.readline(prompt) {
            Ok(line) => line,
            Err(ReadlineError::Interrupted) => continue,
            Err(ReadlineError::Eof) => break,
            Err(e) => {
                eprintln!("input error: {}", e);
                break;
            }
        };
        if line.trim().is_empty() {
            continue;
        }
        let _ = rl.add_history_entry(line.as_str());

        match dispatch(&line, &mut ctx) {
            Ok(Some(s)) => {
                println!("{}", s);
                use std::io::Write;
                let _ = stdout.flush();
            }
            Ok(None) => {}
            Err(e) => eprintln!("error: {}", e),
        }
    }
    let _ = rl.save_history(".mathr_history");
    Ok(())
}

fn dispatch(line: &str, ctx: &mut Context) -> Result<Option<String>> {
    dispatch_inner(line, ctx)
}

/// Dispatch a single input string against a context.
/// Public so the CLI binary can reuse the same smart-dispatch logic.
pub fn dispatch_str(line: &str, mut ctx: Context) -> Result<Option<String>> {
    dispatch_inner(line, &mut ctx)
}

fn dispatch_inner(line: &str, ctx: &mut Context) -> Result<Option<String>> {
    let line = line.trim();
    if line == "quit" || line == "exit" {
        std::process::exit(0);
    }
    if line == "help" || line == "?" {
        return Ok(Some(HELP.to_string()));
    }
    if line == "clear" {
        *ctx = Context::standard();
        return Ok(Some("context cleared".into()));
    }
    if line == "vars" {
        return Ok(Some(list_vars(ctx)));
    }
    if line == "funcs" {
        return Ok(Some(list_funcs(ctx)));
    }

    if let Some(rest) = line.strip_prefix("let ") {
        return bind_var(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("fn ") {
        return define_fn(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("diff ") {
        return do_diff(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("simplify ") {
        let e = Parser::parse(rest.trim())?;
        return Ok(Some(simplify(&e).to_string()));
    }
    if let Some(rest) = line.strip_prefix("int ") {
        return do_integrate(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("solve ") {
        return do_solve(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("plot ") {
        return do_plot(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("fft ") {
        return do_fft(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("taylor ") {
        return do_taylor(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("gcd ") {
        return do_numtheory(rest.trim(), "gcd");
    }
    if let Some(rest) = line.strip_prefix("lcm ") {
        return do_numtheory(rest.trim(), "lcm");
    }
    if let Some(rest) = line.strip_prefix("is-prime ") {
        return do_numtheory(rest.trim(), "is-prime");
    }
    if let Some(rest) = line.strip_prefix("factor ") {
        return do_numtheory(rest.trim(), "factor");
    }
    if let Some(rest) = line.strip_prefix("fib ") {
        return do_numtheory(rest.trim(), "fib");
    }
    if let Some(rest) = line.strip_prefix("binom ") {
        return do_numtheory(rest.trim(), "binom");
    }
    if let Some(rest) = line.strip_prefix("fact ") {
        return do_numtheory(rest.trim(), "fact");
    }
    if let Some(rest) = line.strip_prefix("mr-prime ") {
        return do_numtheory(rest.trim(), "mr-prime");
    }
    if let Some(rest) = line.strip_prefix("conv ") {
        return do_conv(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("stats ") {
        return do_stats(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("poly-roots ") {
        return do_poly_roots(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("isolate-roots ") {
        return do_isolate_roots(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("lu ") {
        return do_lu(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("rank ") {
        return do_rank(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("spline ") {
        return do_spline(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("jacobi ") {
        return do_numtheory(rest.trim(), "jacobi");
    }
    if let Some(rest) = line.strip_prefix("cf ") {
        return do_numtheory(rest.trim(), "cf");
    }
    if let Some(rest) = line.strip_prefix("diophantine ") {
        return do_numtheory(rest.trim(), "diophantine");
    }
    if let Some(rest) = line.strip_prefix("dlog ") {
        return do_numtheory(rest.trim(), "dlog");
    }
    if let Some(rest) = line.strip_prefix("cholesky ") {
        return do_cholesky(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("eig ") {
        return do_eig(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("symlig ") {
        return do_symlig(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("hessenberg ") {
        return do_hessenberg(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("schur ") {
        return do_schur(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("chebyshev ") {
        return do_chebyshev(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("romberg ") {
        return do_romberg(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("svd ") {
        return do_svd(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("legendre ") {
        return do_legendre(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("integrate ") {
        return do_integrate_sym(rest.trim(), ctx);
    }

    // Default: evaluate the expression and print the value
    let e = Parser::parse(line)?;
    let v = eval(&e, ctx)?;
    Ok(Some(format_value(v)))
}

fn bind_var(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    // expects `name = expr`
    let eq_pos = rest
        .find('=')
        .ok_or_else(|| crate::error::MathError::Eval("`let` expects `name = expr`".into()))?;
    let name = rest[..eq_pos].trim();
    if name.is_empty() || !name.chars().all(|c| c.is_alphanumeric() || c == '_') {
        return Err(crate::error::MathError::Eval(format!("invalid variable name: {}", name)));
    }
    let e = Parser::parse(rest[eq_pos + 1..].trim())?;
    let v = eval(&e, ctx)?;
    ctx.set(name, v);
    Ok(Some(format!("{} = {}", name, format_value(v))))
}

fn define_fn(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    // expects `name(args) = expr`
    let eq_pos = rest
        .find('=')
        .ok_or_else(|| crate::error::MathError::Eval("`fn` expects `name(args) = expr`".into()))?;
    let lhs = rest[..eq_pos].trim();
    let rhs = rest[eq_pos + 1..].trim();
    let open = lhs.find('(').ok_or_else(|| crate::error::MathError::Eval("missing `(` in `fn`".into()))?;
    let close = lhs.rfind(')').ok_or_else(|| crate::error::MathError::Eval("missing `)` in `fn`".into()))?;
    let name = lhs[..open].trim();
    let params: Vec<String> = lhs[open + 1..close]
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    let body = Parser::parse(rhs)?;
    ctx.define(name, body, params.clone());
    Ok(Some(format!("defined {}({})", name, params.join(", "))))
}

fn do_diff(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`diff` needs an expression".into()));
    }
    // Parse from the right: optional var (single alpha), rest is expr
    let mut wrt = "x".to_string();
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        let candidate = tokens[expr_end - 1];
        if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
            wrt = candidate.to_string();
            expr_end -= 1;
        }
    }
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    // Substitute any user-bound variables so we don't print unwieldy
    // intermediate forms. Variables still appear if no binding exists.
    let bound: Vec<(String, Expr)> = ctx
        .vars
        .iter()
        .filter(|(k, _)| **k != wrt && !["pi", "e", "tau", "inf"].contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), Expr::num(*v)))
        .collect();
    let mut e = e;
    for (k, v) in &bound {
        e = e.substitute(k, v);
    }
    let d = differentiate(&e, &wrt)?;
    let s = simplify(&d);
    Ok(Some(s.to_string()))
}

fn do_integrate(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval("`int` needs: <expr> <a> <b>".into()));
    }
    // The last two tokens are bounds a, b (can be expressions like "pi").
    // The rest is the expression to integrate.
    // Try different split points in case the expression has trailing numbers.
    for split in [2usize, 3, 4] {
        if tokens.len() <= split {
            continue;
        }
        let b_src = tokens[tokens.len() - 1];
        let a_src = tokens[tokens.len() - 2];
        let expr_src = tokens[..tokens.len() - split].join(" ");
        let e = match Parser::parse(&expr_src) {
            Ok(e) => e,
            Err(_) => continue,
        };
        // Evaluate bounds (support constants like pi, e)
        let a_e = match Parser::parse(a_src) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let b_e = match Parser::parse(b_src) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let a = crate::eval::eval(&a_e, ctx)?;
        let b = crate::eval::eval(&b_e, ctx)?;

        let bound: Vec<(String, Expr)> = ctx
            .vars
            .iter()
            .filter(|(k, _)| !["pi", "e", "tau", "inf"].contains(&k.as_str()))
            .map(|(k, v)| (k.clone(), Expr::num(*v)))
            .collect();
        let mut e = e;
        for (k, v) in &bound {
            e = e.substitute(k, v);
        }
        let ctx2 = ctx.clone();
        let f = move |x: f64| {
            let mut cx = ctx2.clone();
            cx.set("x", x);
            crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
        };
        let v = crate::calculus::integrate_adaptive(f, a, b, 1e-9, 30)?;
        return Ok(Some(format!("∫ = {}", format_value(v))));
    }
    Err(crate::error::MathError::Eval(format!("could not parse: {}", rest)))
}

fn do_solve(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`solve` needs an expression".into()));
    }
    // Try parsing from the right: optional guess (number), optional var (single alpha).
    // If the expression fails to parse with trailing tokens consumed, retry with fewer consumed.
    for consume in [2usize, 1, 0] {
        if tokens.len() <= consume {
            continue;
        }
        let mut guess = 1.0;
        let mut wrt = "x".to_string();
        let mut expr_end = tokens.len();

        if consume >= 1 {
            if let Ok(g) = tokens[expr_end - 1].parse::<f64>() {
                guess = g;
                expr_end -= 1;
            } else {
                continue;
            }
        }
        if consume >= 2 {
            let candidate = tokens[expr_end - 1];
            if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
                wrt = candidate.to_string();
                expr_end -= 1;
            } else {
                continue;
            }
        }
        let expr_src = tokens[..expr_end].join(" ");
        let e = match Parser::parse(&expr_src) {
            Ok(e) => e,
            Err(_) => continue,
        };
        let ctx2 = ctx.clone();
        let f = move |x: f64| {
            let mut cx = ctx2.clone();
            cx.set(&wrt, x);
            crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
        };
        let (root, fval) = crate::solver::newton_central(f, guess, crate::solver::SolveOptions::default())?;
        return Ok(Some(format!("root ≈ {} (f = {})", format_value(root), format_value(fval))));
    }
    Err(crate::error::MathError::Eval(format!("could not parse: {}", rest)))
}

fn do_plot(rest: &str, _ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`plot` needs an expression".into()));
    }
    // Parse from the right: optional filename, b (number), a (number), rest is expr
    let mut expr_end = tokens.len();
    let file = if expr_end > 3 && tokens[expr_end - 1].ends_with(".png") {
        expr_end -= 1;
        tokens[expr_end].to_string()
    } else {
        "plot.png".to_string()
    };
    if expr_end < 3 {
        return Err(crate::error::MathError::Eval("`plot` needs: <expr> <a> <b> [out.png]".into()));
    }
    let b: f64 = tokens[expr_end - 1].parse()
        .map_err(|_| crate::error::MathError::Eval("could not parse b".into()))?;
    let a: f64 = tokens[expr_end - 2].parse()
        .map_err(|_| crate::error::MathError::Eval("could not parse a".into()))?;
    let expr_src = tokens[..expr_end - 2].join(" ");

    let wrt = guess_var(&expr_src);
    let e = Parser::parse(&expr_src)?;
    crate::plot::plot_function(&file, &e, &wrt, a, b, 800, &format!("y = {}", expr_src))?;
    Ok(Some(format!("wrote {}", file)))
}

fn do_fft(rest: &str) -> Result<Option<String>> {
    let mut samples: Vec<f64> = Vec::new();
    for tok in rest.split(|c: char| c == ',' || c.is_whitespace()) {
        let tok = tok.trim();
        if tok.is_empty() {
            continue;
        }
        match tok.parse::<f64>() {
            Ok(v) => samples.push(v),
            Err(_) => return Err(crate::error::MathError::Eval(format!("not a number: {}", tok))),
        }
    }
    let mags = crate::fft::magnitude_spectrum(&samples)?;
    let mut out = String::new();
    for (k, m) in mags.iter().enumerate() {
        out.push_str(&format!("X[{}] = {:.4}\n", k, m));
    }
    Ok(Some(out))
}

fn do_taylor(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`taylor` needs an expression".into()));
    }
    // Parse from the right: optional order (int), optional around (float), rest is expr
    let mut around = 0.0;
    let mut order = 5usize;
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        if let Ok(o) = tokens[expr_end - 1].parse::<usize>() {
            order = o;
            expr_end -= 1;
        }
    }
    if expr_end > 1 {
        if let Ok(a) = tokens[expr_end - 1].parse::<f64>() {
            around = a;
            expr_end -= 1;
        }
    }
    let expr_src = tokens[..expr_end].join(" ");
    let series = crate::taylor::taylor_series_str(&expr_src, "x", around, order)?;
    Ok(Some(series.to_string()))
}

fn do_numtheory(rest: &str, op: &str) -> Result<Option<String>> {
    use crate::numtheory;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    match op {
        "gcd" => {
            let nums = parse_u64_list(&tokens)?;
            if nums.len() < 2 { return Err(crate::error::MathError::Eval("gcd needs ≥2 numbers".into())); }
            let g = nums.iter().copied().reduce(numtheory::gcd).unwrap();
            Ok(Some(g.to_string()))
        }
        "lcm" => {
            let nums = parse_u64_list(&tokens)?;
            if nums.len() < 2 { return Err(crate::error::MathError::Eval("lcm needs ≥2 numbers".into())); }
            let l = nums.iter().copied().reduce(numtheory::lcm).unwrap();
            Ok(Some(l.to_string()))
        }
        "is-prime" => {
            let n = parse_u64_single(&tokens)?;
            Ok(Some(numtheory::is_prime(n).to_string()))
        }
        "factor" => {
            let n = parse_u64_single(&tokens)?;
            let factors = numtheory::prime_factors(n);
            if factors.is_empty() {
                Ok(Some("1".into()))
            } else {
                let strs: Vec<String> = factors.iter().map(|f| f.to_string()).collect();
                Ok(Some(strs.join(" * ")))
            }
        }
        "fib" => {
            let n = parse_u64_single(&tokens)?;
            Ok(Some(numtheory::fibonacci(n).to_string()))
        }
        "binom" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("binom needs n k".into())); }
            let n: u64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
            let k: u64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad k".into()))?;
            let r = numtheory::binomial(n, k)?;
            Ok(Some(r.to_string()))
        }
        "fact" => {
            let n = parse_u64_single(&tokens)?;
            let r = numtheory::factorial(n)?;
            Ok(Some(r.to_string()))
        }
        "mr-prime" => {
            let n = parse_u64_single(&tokens)?;
            let rounds: usize = tokens.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
            Ok(Some(numtheory::is_prime_miller_rabin(n, rounds).to_string()))
        }
        "jacobi" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("jacobi needs a n".into())); }
            let a: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
            let n: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
            Ok(Some(numtheory::jacobi_symbol(a, n).to_string()))
        }
        "cf" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("cf needs p q".into())); }
            let p: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad p".into()))?;
            let q: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad q".into()))?;
            let cf = numtheory::continued_fraction(p, q)?;
            let strs: Vec<String> = cf.iter().map(|v| v.to_string()).collect();
            Ok(Some(format!("[{}]", strs.join("; "))))
        }
        "diophantine" => {
            if tokens.len() < 3 { return Err(crate::error::MathError::Eval("diophantine needs a b c".into())); }
            let a: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
            let b: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad b".into()))?;
            let c: i64 = tokens[2].parse().map_err(|_| crate::error::MathError::Eval("bad c".into()))?;
            let (x, y) = numtheory::diophantine(a, b, c)?;
            Ok(Some(format!("x = {}, y = {} (a*x + b*y = {})", x, y, a * x + b * y)))
        }
        "dlog" => {
            if tokens.len() < 3 { return Err(crate::error::MathError::Eval("dlog needs g h p".into())); }
            let g: u64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad g".into()))?;
            let h: u64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad h".into()))?;
            let p: u64 = tokens[2].parse().map_err(|_| crate::error::MathError::Eval("bad p".into()))?;
            match numtheory::discrete_log(g, h, p) {
                Some(x) => Ok(Some(format!("x = {} (g^x mod p = {})", x, numtheory::mod_pow(g, x, p)))),
                None => Ok(Some("(no discrete log found)".into())),
            }
        }
        _ => Err(crate::error::MathError::Eval(format!("unknown op: {}", op))),
    }
}

fn do_cholesky(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let c = m.cholesky()?;
    let l = c.l_factor();
    let reconstructed = c.reconstruct();
    let rows = m.rows;
    let cols = m.cols;
    let mut max_diff = 0.0_f64;
    for i in 0..rows {
        for j in 0..cols {
            let d = (m[(i, j)] - reconstructed[(i, j)]).abs();
            if d > max_diff {
                max_diff = d;
            }
        }
    }
    Ok(Some(format!(
        "cholesky ok (max reconstruction error = {:.2e})\nL =\n{}",
        max_diff, l
    )))
}

fn do_eig(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let result = m.power_iteration(crate::matrix::PowerIterOptions::default())?;
    let v_strs: Vec<String> = result.vector.iter().map(|x| format!("{:.6}", x)).collect();
    Ok(Some(format!(
        "λ ≈ {} (dominant eigenvalue)\nv = [{}]",
        format_value(result.value),
        v_strs.join(", ")
    )))
}

fn do_symlig(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let (eigenvalues, eigenvectors) = m.symmetric_eig()?;
    let val_strs: Vec<String> = eigenvalues.iter().map(|x| format_value(*x)).collect();
    let mut rows: Vec<String> = Vec::new();
    for i in 0..eigenvectors.rows {
        let row: Vec<String> = (0..eigenvectors.cols)
            .map(|j| format!("{:.6}", eigenvectors[(i, j)]))
            .collect();
        rows.push(format!("[{}]", row.join(", ")));
    }
    Ok(Some(format!(
        "eigenvalues (ascending):\n  [{}]\neigenvectors (columns):\n{}",
        val_strs.join(", "),
        rows.join("\n")
    )))
}

fn do_hessenberg(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let (h, q) = m.hessenberg()?;
    let mut out = String::from("H (upper Hessenberg):\n");
    for i in 0..h.rows {
        let row: Vec<String> = (0..h.cols).map(|j| format!("{:.6}", h[(i, j)])).collect();
        out.push_str(&format!("[{}]\n", row.join(", ")));
    }
    out.push_str("Q (orthogonal):\n");
    for i in 0..q.rows {
        let row: Vec<String> = (0..q.cols).map(|j| format!("{:.6}", q[(i, j)])).collect();
        out.push_str(&format!("[{}]\n", row.join(", ")));
    }
    Ok(Some(out))
}

fn do_schur(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let (t, q) = m.schur()?;
    let mut out = String::from("T (quasi-upper triangular):\n");
    for i in 0..t.rows {
        let row: Vec<String> = (0..t.cols).map(|j| format!("{:.6}", t[(i, j)])).collect();
        out.push_str(&format!("[{}]\n", row.join(", ")));
    }
    out.push_str("Q (orthogonal):\n");
    for i in 0..q.rows {
        let row: Vec<String> = (0..q.cols).map(|j| format!("{:.6}", q[(i, j)])).collect();
        out.push_str(&format!("[{}]\n", row.join(", ")));
    }
    Ok(Some(out))
}

fn do_chebyshev(rest: &str) -> Result<Option<String>> {
    // Format: "n x" — Chebyshev T_n(x). Or "n" alone for an array of nodes.
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("chebyshev needs n [x]".into()));
    }
    let n: u32 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
    if tokens.len() == 1 {
        let nodes = crate::interpolate::chebyshev_nodes(n as usize);
        let strs: Vec<String> = nodes.iter().map(|x| format!("{:.6}", x)).collect();
        return Ok(Some(format!("T_{} nodes: [{}]", n, strs.join(", "))));
    }
    let x: f64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad x".into()))?;
    Ok(Some(format!(
        "T_{}({}) = {}",
        n, x, format_value(crate::interpolate::chebyshev_t(n, x))
    )))
}

fn do_romberg(rest: &str) -> Result<Option<String>> {
    // Format: "<expr> a b [levels=8]"
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval("romberg needs: <expr> a b [levels]".into()));
    }
    let b: f64 = tokens[tokens.len() - 1].parse()
        .map_err(|_| crate::error::MathError::Eval("bad b".into()))?;
    let a: f64 = tokens[tokens.len() - 2].parse()
        .map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
    let levels: usize = if tokens.len() >= 4 {
        tokens[tokens.len() - 3].parse().unwrap_or(8)
    } else {
        8
    };
    let levels = if tokens.len() == 3 { 8 } else if tokens.len() == 4 { 8 } else { levels };
    let expr_end = if tokens.len() >= 4 { tokens.len() - 3 } else { tokens.len() - 2 };
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    let ctx2 = Context::standard();
    let f = move |x: f64| {
        let mut cx = ctx2.clone();
        cx.set("x", x);
        crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
    };
    let v = crate::calculus::integrate_romberg(f, a, b, levels)?;
    Ok(Some(format!("∫ ({}, {}; levels = {}) = {}", format_value(a), format_value(b), levels, format_value(v))))
}

fn do_svd(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let svd = m.svd()?;
    let s_strs: Vec<String> = svd.singular_values.iter().map(|s| format!("{:.6}", s)).collect();
    Ok(Some(format!(
        "σ = [{}]\nU = {} rows × {} cols\nV = {}",
        s_strs.join(", "),
        svd.u.rows, svd.u.cols, svd.v
    )))
}

fn do_legendre(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("legendre needs: n [x]".into()));
    }
    let n: u32 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
    if tokens.len() == 1 {
        // Show first few Gauss–Legendre nodes for order n.
        let (nodes, weights) = crate::interpolate::gauss_legendre(n as usize);
        let n_strs: Vec<String> = nodes.iter().map(|x| format!("{:.4}", x)).collect();
        let w_strs: Vec<String> = weights.iter().map(|w| format!("{:.4}", w)).collect();
        return Ok(Some(format!(
            "Gauss–Legendre {} nodes:\n  x = [{}]\n  w = [{}]",
            n, n_strs.join(", "), w_strs.join(", ")
        )));
    }
    let x: f64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad x".into()))?;
    Ok(Some(format!(
        "P_{}({}) = {}",
        n, x, format_value(crate::interpolate::legendre_p(n, x))
    )))
}

fn do_integrate_sym(rest: &str, _ctx: &mut Context) -> Result<Option<String>> {
    // Format: "<expr> [var]" — symbolic integration.
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("integrate needs: <expr> [var]".into()));
    }
    let var = if tokens.len() >= 2 {
        tokens[tokens.len() - 1].to_string()
    } else {
        "x".to_string()
    };
    let expr_end = if tokens.len() >= 2
        && tokens[tokens.len() - 1].len() == 1
        && tokens[tokens.len() - 1].chars().all(|c| c.is_ascii_alphabetic())
    {
        tokens.len() - 1
    } else {
        tokens.len()
    };
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    let result = crate::symbolic::integrate(&e, &var)?;
    Ok(Some(format!("∫ {} d{} = {}", expr_src, var, result)))
}

fn parse_u64_single(tokens: &[&str]) -> Result<u64> {
    tokens.first()
        .ok_or_else(|| crate::error::MathError::Eval("missing number".into()))?
        .parse::<u64>()
        .map_err(|_| crate::error::MathError::Eval("could not parse as integer".into()))
}

fn parse_u64_list(tokens: &[&str]) -> Result<Vec<u64>> {
    tokens.iter()
        .map(|s| s.parse::<u64>().map_err(|_| crate::error::MathError::Eval(format!("not an integer: {}", s))))
        .collect()
}

fn do_conv(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    let split_pos = tokens.iter().position(|t| *t == "x" || *t == "X")
        .ok_or_else(|| crate::error::MathError::Eval("conv needs 'x' separator: conv 1 2 3 x 1 1".into()))?;
    let a: Vec<f64> = tokens[..split_pos].iter()
        .map(|s| s.parse::<f64>().map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", s))))
        .collect::<Result<_>>()?;
    let b: Vec<f64> = tokens[split_pos + 1..].iter()
        .map(|s| s.parse::<f64>().map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", s))))
        .collect::<Result<_>>()?;
    let c = crate::fft::convolve(&a, &b)?;
    let strs: Vec<String> = c.iter().map(|v| format!("{:.6}", v)).collect();
    Ok(Some(strs.join("  ")))
}

fn do_stats(rest: &str) -> Result<Option<String>> {
    let nums = parse_f64_list(rest)?;
    if nums.is_empty() { return Err(crate::error::MathError::Eval("stats needs numbers".into())); }
    let s = crate::stats::summary(&nums)?;
    let mut out = format!(
        "count={}\nmean={}\nmedian={}\nstddev={}\nmin={}\nmax={}\nrange={}",
        s.count, format_value(s.mean), format_value(s.median),
        format_value(s.stddev), format_value(s.min), format_value(s.max), format_value(s.range)
    );
    if let Ok(v) = crate::stats::variance_sample(&nums) {
        out.push_str(&format!("\nvar(s)={}", format_value(v)));
    }
    if let Ok((q1, q2, q3)) = crate::stats::quartiles(&nums) {
        out.push_str(&format!("\nQ1={}\nQ2={}\nQ3={}\nIQR={}", format_value(q1), format_value(q2), format_value(q3), format_value(q3 - q1)));
    }
    Ok(Some(out))
}

fn do_poly_roots(rest: &str) -> Result<Option<String>> {
    let coeffs = parse_f64_list(rest)?;
    if coeffs.is_empty() { return Err(crate::error::MathError::Eval("poly-roots needs coefficients".into())); }
    let r = crate::solver::polynomial_roots(&coeffs)?;
    if r.is_empty() {
        Ok(Some("(no real roots found)".into()))
    } else {
        let lines: Vec<String> = r.iter().map(|(x, fx)| format!("x = {}  (f = {})", format_value(*x), format_value(*fx))).collect();
        Ok(Some(lines.join("\n")))
    }
}

fn do_isolate_roots(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("isolate-roots needs integer coefficients".into()));
    }
    let coeffs: Vec<i64> = tokens.iter()
        .map(|t| t.parse::<i64>().map_err(|_| crate::error::MathError::Eval(format!("invalid integer: {}", t))))
        .collect::<crate::error::Result<_>>()?;
    let intervals = crate::solver::isolate_real_roots(&coeffs)?;
    if intervals.is_empty() {
        Ok(Some("(no real roots found)".into()))
    } else {
        let lines: Vec<String> = intervals.iter().map(|(lo, hi)| {
            if (lo - hi).abs() < 1e-10 {
                format!("x = {}", format_value(*lo))
            } else {
                format!("x in ({}, {})", format_value(*lo), format_value(*hi))
            }
        }).collect();
        Ok(Some(lines.join("\n")))
    }
}

/// Parse a matrix from `rest` formatted as "1 2 3 | 4 5 6" (rows separated by |).
fn parse_matrix(rest: &str) -> Result<crate::matrix::Matrix> {
    let mut rows: Vec<Vec<f64>> = Vec::new();
    for row_src in rest.split('|') {
        let row_src = row_src.trim();
        if row_src.is_empty() {
            continue;
        }
        let row: Vec<f64> = row_src
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.trim().is_empty())
            .map(|s| {
                s.trim().parse::<f64>()
                    .map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", s)))
            })
            .collect::<Result<_>>()?;
        rows.push(row);
    }
    crate::matrix::Matrix::from_rows(&rows)
}

fn do_lu(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let fact = m.lu()?;
    let det = fact.determinant();
    let inv = fact.inverse()?;
    Ok(Some(format!(
        "det = {}\nA⁻¹ =\n{}",
        format_value(det),
        inv
    )))
}

fn do_rank(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let r = m.rank(1e-10);
    Ok(Some(format!("rank = {}", r)))
}

fn do_spline(rest: &str) -> Result<Option<String>> {
    // Format: "x1 y1 x2 y2 ...  [eval x_value]"
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 || tokens.len() % 2 == 0 {
        return Err(crate::error::MathError::Eval(
            "spline expects: x1 y1 x2 y2 ... [x_at]".into(),
        ));
    }
    let pts: Vec<(f64, f64)> = (0..(tokens.len() / 2))
        .map(|i| {
            let x: f64 = tokens[2 * i].parse().map_err(|_| {
                crate::error::MathError::Eval(format!("bad x at {}", tokens[2 * i]))
            })?;
            let y: f64 = tokens[2 * i + 1].parse().map_err(|_| {
                crate::error::MathError::Eval(format!("bad y at {}", tokens[2 * i + 1]))
            })?;
            Ok((x, y))
        })
        .collect::<Result<_>>()?;
    let sp = crate::interpolate::CubicSpline::new(&pts)?;
    let last = tokens.last().unwrap();
    let x_at: f64 = last.parse().map_err(|_| {
        crate::error::MathError::Eval(format!("bad x_at: {}", last))
    })?;
    Ok(Some(format!("spline({}) = {}", x_at, format_value(sp.eval(x_at)))))
}

fn parse_f64_list(rest: &str) -> Result<Vec<f64>> {
    rest.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.trim().is_empty())
        .map(|t| t.trim().parse::<f64>().map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", t))))
        .collect()
}

fn guess_var(expr: &str) -> String {
    // Pick the first single-letter identifier as the plotting variable.
    let mut seen = HashSet::new();
    let mut candidates = Vec::new();
    let mut chars = expr.chars().peekable();
    while let Some(c) = chars.next() {
        if c.is_ascii_alphabetic() && c.is_ascii_lowercase() && !seen.contains(&c) {
            let mut name = c.to_string();
            while let Some(&c2) = chars.peek() {
                if c2.is_ascii_alphanumeric() {
                    name.push(c2);
                    chars.next();
                } else {
                    break;
                }
            }
            if name.len() == 1 {
                candidates.push(name.clone());
                seen.insert(c);
            }
        }
    }
    candidates.first().cloned().unwrap_or_else(|| "x".into())
}

fn list_vars(ctx: &Context) -> String {
    if ctx.vars.is_empty() {
        return "(no variables)".into();
    }
    let mut out = String::from("variables:\n");
    let mut names: Vec<&String> = ctx.vars.keys().collect();
    names.sort();
    for n in names {
        out.push_str(&format!("  {:8} = {}\n", n, format_value(ctx.vars[n])));
    }
    out
}

fn list_funcs(ctx: &Context) -> String {
    if ctx.funcs.is_empty() {
        return "(no functions)".into();
    }
    let mut out = String::from("functions:\n");
    let mut names: Vec<&String> = ctx.funcs.keys().collect();
    names.sort();
    for n in names {
        let kind = match &ctx.funcs[n] {
            Func::Builtin(_) => "builtin",
            Func::User(_, p) => {
                let _ = p;
                "user"
            }
        };
        out.push_str(&format!("  {:8} ({})\n", n, kind));
    }
    out
}

fn format_value(v: f64) -> String {
    if v.is_nan() {
        "NaN".into()
    } else if v.is_infinite() {
        if v > 0.0 { "inf".into() } else { "-inf".into() }
    } else if v == v.trunc() && v.abs() < 1e15 {
        format!("{}", v as i64)
    } else {
        format!("{:.10}", v)
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    }
}

const HELP: &str = "\
commands:
  <expr>              evaluate (e.g. sin(pi/4), gamma(0.5), bessel_j0(5))
  let x = <expr>      bind a variable
  fn f(x) = <expr>    define a function
  diff <expr> [var]   symbolic derivative
  simplify <expr>     constant-fold & simplify
  int <expr> a b      numerical integral over [a, b]
  solve <expr> [var] [guess]
  plot <expr> a b [out.png]
  taylor <expr> [a] [order]
                      Taylor series around a (default 0, order 5)
  fft <numbers...>    magnitude spectrum of a real signal
  conv <a...> x <b...>  convolution of two signals
  stats <numbers...>  descriptive statistics
  poly-roots <coeffs...>  polynomial roots (highest degree first)
  isolate-roots <ints...>  real root isolation (VAS, integer coefficients)
  lu <rows>           matrix LU decomposition (rows separated by '|')
                      computes determinant and inverse
  rank <rows>         matrix rank
  spline x1 y1 x2 y2 ... x_at
                      natural cubic spline interpolant at x_at
  cholesky <rows>     Cholesky decomposition of symmetric positive-definite
  eig <rows>          dominant eigenvalue + eigenvector (power iteration)
  svd <rows>          singular value decomposition
  chebyshev n [x]     Chebyshev T_n(x), or T_n nodes on [-1, 1]
  legendre n [x]      Legendre P_n(x), or Gauss–Legendre n-node weights
  integrate <expr> [var]
                      symbolic integration of <expr> with respect to var
  romberg <expr> a b  numerical integral via Romberg (Richardson)
  gcd <n...>          greatest common divisor
  lcm <n...>          least common multiple
  is-prime <n>        primality test
  factor <n>          prime factorization
  fib <n>             Fibonacci number
  binom <n> <k>       binomial coefficient
  fact <n>            factorial
  mr-prime <n> [r]    Miller–Rabin primality test
  jacobi <a> <n>      Jacobi symbol (a/n)
  cf <p> <q>          continued fraction of p/q
  diophantine <a> <b> <c>  solve a*x + b*y = c
  dlog <g> <h> <p>    discrete logarithm x: g^x ≡ h (mod p)
  vars / funcs        show bindings
  clear               reset context
  help                this help
  quit                leave the REPL
constants: pi, e, tau, inf
functions: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh,
           exp, ln, log, log2, log10, sqrt, cbrt, abs, floor,
           ceil, round, sign, min, max, pow, mod, fract,
           gamma, erf, erfc, sinc, bessel_j0, bessel_j1, bessel_j
";

/// Suppresses an unused warning for the `Cow` import in environments where
/// rustc may otherwise complain; the type is used implicitly through Helper.
#[allow(dead_code)]
fn _cow_silence<'a>(_: Cow<'a, str>) {}