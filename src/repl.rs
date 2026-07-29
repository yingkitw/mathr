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
    highlighter: MatchingBracketHighlighter,
    #[rustyline(Validator)]
    validator: MatchingBracketValidator,
    #[rustyline(Hinter)]
    hinter: HistoryHinter,
}

impl Default for ReplHelper {
    fn default() -> Self {
        Self {
            completer: FilenameCompleter::new(),
            highlighter: MatchingBracketHighlighter::new(),
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
        _ => Err(crate::error::MathError::Eval(format!("unknown op: {}", op))),
    }
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
  <expr>              evaluate (e.g. sin(pi/4), gamma(0.5), erf(1.0))
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
  gcd <n...>          greatest common divisor
  lcm <n...>          least common multiple
  is-prime <n>        primality test
  factor <n>          prime factorization
  fib <n>             Fibonacci number
  binom <n> <k>       binomial coefficient
  fact <n>            factorial
  mr-prime <n> [r]    Miller–Rabin primality test
  vars / funcs        show bindings
  clear               reset context
  help                this help
  quit                leave the REPL
constants: pi, e, tau, inf
functions: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh,
           exp, ln, log, log2, log10, sqrt, cbrt, abs, floor,
           ceil, round, sign, min, max, pow, mod, fract,
           gamma, erf, erfc, sinc
";

/// Suppresses an unused warning for the `Cow` import in environments where
/// rustc may otherwise complain; the type is used implicitly through Helper.
#[allow(dead_code)]
fn _cow_silence<'a>(_: Cow<'a, str>) {}