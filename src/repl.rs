//! An interactive REPL that ties together the parser, evaluator, symbolic
//! differentiation, numeric integration, root-finders and plotter.
//!
//! Commands available at the prompt:
//!   <expr>             evaluate
//!   let x = 2          bind a variable
//!   fn f(x) = x^2      define a function
//!   diff <expr> [wrt]  symbolic derivative
//!   int <expr> a b     numerical integral
//!   solve <expr> [wrt] [guess]
//!   simplify <expr>
//!   plot <expr> a b [out.png]
//!   vars | funcs | clear | help | quit

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
    let _ = rl.load_history(".maths_history");
    println!("maths {} — type `help` for a list of commands.", env!("CARGO_PKG_VERSION"));

    let mut ctx = Context::standard();
    let mut stdout = std::io::stdout();

    loop {
        let prompt = if ctx.vars.is_empty() {
            "\nmaths> "
        } else {
            "\nmaths* "
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
    let _ = rl.save_history(".maths_history");
    Ok(())
}

fn dispatch(line: &str, ctx: &mut Context) -> Result<Option<String>> {
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
    let mut parts = rest.split_whitespace();
    let expr_src = parts.next().ok_or_else(|| crate::error::MathError::Eval("`diff` needs an expression".into()))?;
    let wrt = parts.next().unwrap_or("x");
    let e = Parser::parse(expr_src)?;
    // Substitute any user-bound variables so we don't print unwieldy
    // intermediate forms. Variables still appear if no binding exists.
    let bound: Vec<(String, Expr)> = ctx
        .vars
        .iter()
        .filter(|(k, _)| *k != wrt && !["pi", "e", "tau", "inf"].contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), Expr::num(*v)))
        .collect();
    let mut e = e;
    for (k, v) in &bound {
        e = e.substitute(k, v);
    }
    let d = differentiate(&e, wrt)?;
    let s = simplify(&d);
    Ok(Some(s.to_string()))
}

fn do_integrate(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let mut parts = rest.split_whitespace();
    let expr_src = parts.next().ok_or_else(|| crate::error::MathError::Eval("`int` needs an expression".into()))?;
    let a: f64 = parts
        .next()
        .ok_or_else(|| crate::error::MathError::Eval("`int` needs bounds `a b`".into()))?
        .parse()
        .map_err(|_| crate::error::MathError::Eval("could not parse a".into()))?;
    let b: f64 = parts
        .next()
        .ok_or_else(|| crate::error::MathError::Eval("`int` needs bounds `a b`".into()))?
        .parse()
        .map_err(|_| crate::error::MathError::Eval("could not parse b".into()))?;

    let e = Parser::parse(expr_src)?;
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
    Ok(Some(format!("∫ = {}", format_value(v))))
}

fn do_solve(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let mut parts = rest.split_whitespace();
    let expr_src = parts.next().ok_or_else(|| crate::error::MathError::Eval("`solve` needs an expression".into()))?;
    let wrt = parts.next().unwrap_or("x");
    let guess: f64 = parts
        .next()
        .map(|s| s.parse().unwrap_or(0.0))
        .unwrap_or(0.0);
    let e = Parser::parse(expr_src)?;
    let ctx2 = ctx.clone();
    let f = move |x: f64| {
        let mut cx = ctx2.clone();
        cx.set(wrt, x);
        crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
    };
    let (root, fval) = crate::solver::newton_central(f, guess, crate::solver::SolveOptions::default())?;
    Ok(Some(format!("root ≈ {} (f = {})", format_value(root), format_value(fval))))
}

fn do_plot(rest: &str, _ctx: &mut Context) -> Result<Option<String>> {
    let mut parts = rest.split_whitespace();
    let expr_src = parts.next().ok_or_else(|| crate::error::MathError::Eval("`plot` needs an expression".into()))?;
    let a: f64 = parts.next().unwrap_or("-6.28".into()).parse().map_err(|_| crate::error::MathError::Eval("could not parse a".into()))?;
    let b: f64 = parts.next().unwrap_or("6.28".into()).parse().map_err(|_| crate::error::MathError::Eval("could not parse b".into()))?;
    let file = parts.next().unwrap_or("plot.png");

    let wrt = guess_var(expr_src);
    let e = Parser::parse(expr_src)?;
    crate::plot::plot_function(file, &e, &wrt, a, b, 800, &format!("y = {}", expr_src))?;
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
  <expr>              evaluate
  let x = <expr>      bind a variable
  fn f(x) = <expr>    define a function
  diff <expr> [var]   symbolic derivative
  simplify <expr>     constant-fold & simplify
  int <expr> a b      numerical integral over [a, b]
  solve <expr> [var] [guess]
  plot <expr> a b [out.png]
  fft <numbers...>    magnitude spectrum of a real signal
  vars / funcs        show bindings
  clear               reset context
  help                this help
  quit                leave the REPL
constants: pi, e, tau, inf
functions: sin, cos, tan, asin, acos, atan, sinh, cosh, tanh,
           exp, ln, log, log2, log10, sqrt, cbrt, abs, floor,
           ceil, round, sign, min, max, pow, mod, fract
";

/// Suppresses an unused warning for the `Cow` import in environments where
/// rustc may otherwise complain; the type is used implicitly through Helper.
#[allow(dead_code)]
fn _cow_silence<'a>(_: Cow<'a, str>) {}