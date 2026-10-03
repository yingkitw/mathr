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
//!   `dec \<expr\> \[prec \<n\>\] \[with \<var\>=\<val\>,...\]`  arbitrary precision
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

/// Dispatch a single input string against a mutable context.
/// Unlike [`dispatch_str`], this mutates the context in-place so that
/// `let` bindings and `fn` definitions persist across calls.
pub fn dispatch_with_ctx(line: &str, ctx: &mut Context) -> Result<Option<String>> {
    dispatch_inner(line, ctx)
}

/// Dispatch a single input and return step-by-step output.
/// Each element of the returned Vec is one step (rendered as a separate line in the UI).
pub fn dispatch_steps(line: &str, ctx: Context) -> Result<Vec<String>> {
    let line = line.trim();
    if line.is_empty() {
        return Ok(vec![]);
    }

    // diff: show original, derivative (unsimplified), simplified
    if let Some(rest) = line.strip_prefix("diff ") {
        return diff_steps(rest.trim(), ctx);
    }

    // simplify: show original and result
    if let Some(rest) = line.strip_prefix("simplify ") {
        let e = Parser::parse(rest.trim())?;
        let s = simplify(&e);
        return Ok(vec![
            format!("simplify: {}", e),
            format!("= {}", s),
        ]);
    }

    // expand: show original and distributed polynomial
    if let Some(rest) = line.strip_prefix("expand ") {
        let e = Parser::parse(rest.trim())?;
        let s = crate::poly::expand(&e);
        return Ok(vec![
            format!("expand: {}", e),
            format!("= {}", s),
        ]);
    }

    // apart: partial fraction decomposition steps
    if let Some(rest) = line.strip_prefix("apart ") {
        return crate::apart::apart_steps_str(rest.trim());
    }

    // solve: show equation, method, root
    if let Some(rest) = line.strip_prefix("solve ") {
        return solve_steps(rest.trim(), ctx);
    }

    // taylor: show expression, expansion point, order, series
    if let Some(rest) = line.strip_prefix("taylor ") {
        return taylor_steps(rest.trim());
    }

    // integrate (symbolic): show integral and result
    if let Some(rest) = line.strip_prefix("integrate ") {
        return integrate_sym_steps(rest.trim());
    }

    // rat: show operands, operation, result
    if let Some(rest) = line.strip_prefix("rat ") {
        return rat_steps(rest.trim());
    }

    // laurent: show expression, center, pole order, series
    if let Some(rest) = line.strip_prefix("laurent ") {
        return laurent_steps(rest.trim());
    }

    // limit: show substitution / L'Hôpital / probe steps
    if let Some(rest) = line.strip_prefix("limit ") {
        return limit_steps_str(rest.trim());
    }

    // qsolve: show equation, coefficients, discriminant, roots
    if let Some(rest) = line.strip_prefix("qsolve ") {
        return crate::qsolve::qsolve_steps_str(rest.trim());
    }

    // idiff: implicit differentiation steps (F_x, F_y, dy/dx)
    if let Some(rest) = line.strip_prefix("idiff ") {
        return crate::symbolic::idiff_steps_str(rest.trim());
    }

    // For all other REPL commands (int, romberg, fft, plot, stats, etc.),
    // fall back to dispatch_inner and wrap the result as a single step.
    let cmd_keywords = [
        "int ", "romberg ", "fft ", "conv ", "plot ", "stats ",
        "poly-roots ", "isolate-roots ", "lu ", "qr ", "cholesky ", "svd ",
        "eig ", "symlig ", "hessenberg ", "schur ", "rank ", "tikhonov ",
        "spline ", "pchip ", "minimize ", "chebyshev ", "legendre ", "fourier ", "mc ",
        "bspline ", "hermite ",
        "sample ", "dist ", "qtile ", "fit ", "ttest ", "chitest ", "anova ", "mwu ", "wilcoxon ",
        "kw ", "boot ", "logit ", "pdiff ", "gradient ", "let ", "fn ",
        "gcd ", "lcm ", "is-prime ", "factor ", "fib ", "binom ",
        "fact ", "mr-prime ", "jacobi ", "cf ", "diophantine ", "dlog ",
        "det ",
        "fast ",
        "cond ",
        "null ",
        "big ",
        "ad ",
        "mathml ",
        "serialize ",
        "interval ",
        "cval ",
        "qsolve ",
        "sum ",
        "prod ",
        "spcg ",
        "spbicg ",
        "idiff ",
    ];
    if cmd_keywords.iter().any(|kw| line.starts_with(kw)) || line == "vars" || line == "funcs" {
        let result = dispatch_inner(line, &mut ctx.clone())?;
        return Ok(match result {
            Some(s) if !s.is_empty() => vec![s],
            _ => vec![],
        });
    }

    // Default: try exact rational evaluation first, fall back to f64, then simplify
    let e = Parser::parse(line)?;
    if let Some(r) = crate::rational::eval_rational(&e) {
        return Ok(vec![format!("= {}", r)]);
    }
    match eval(&e, &ctx) {
        Ok(v) => Ok(vec![format!("= {}", format_value(v))]),
        Err(_) => {
            // Evaluation failed (e.g. unbound variables) — try simplification
            let s = simplify(&e);
            Ok(vec![format!("= {}", s)])
        }
    }
}

fn diff_steps(rest: &str, ctx: Context) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`diff` needs an expression".into()));
    }
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
    Ok(vec![
        format!("f({}) = {}", wrt, e),
        format!("d/d{} f({})", wrt, wrt),
        format!("= {}", d),
        format!("simplified = {}", s),
    ])
}

fn solve_steps(rest: &str, ctx: Context) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`solve` needs an expression".into()));
    }
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
        let wrt_clone = wrt.clone();
        let f = move |x: f64| {
            let mut cx = ctx2.clone();
            cx.set(&wrt_clone, x);
            crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
        };
        let (root, fval) = crate::solver::newton_central(f, guess, crate::solver::SolveOptions::default())?;
        return Ok(vec![
            format!("solve: {}({}) = 0", expr_src, wrt),
            format!("method: Newton-Raphson, initial guess = {}", format_value(guess)),
            format!("root: {} ≈ {}", wrt, format_value(root)),
            format!("residual: f({}) = {}", format_value(root), format_value(fval)),
        ]);
    }
    Err(crate::error::MathError::Eval(format!("could not parse: {}", rest)))
}

fn taylor_steps(rest: &str) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`taylor` needs an expression".into()));
    }
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
    Ok(vec![
        format!("f(x) = {}", expr_src),
        format!("Taylor expansion around a = {}", format_value(around)),
        format!("order = {}", order),
        format!("f(x) ≈ {}", series),
    ])
}

fn integrate_sym_steps(rest: &str) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("integrate needs: <expr> [var]".into()));
    }
    let var = if tokens.len() >= 2
        && tokens[tokens.len() - 1].len() == 1
        && tokens[tokens.len() - 1].chars().all(|c| c.is_ascii_alphabetic())
    {
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
    Ok(vec![
        format!("integrate: {}", expr_src),
        format!("∫ {} d{}", expr_src, var),
        format!("= {} + C", result),
    ])
}

fn rat_steps(rest: &str) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`rat` needs: <a> <op> <b>  (e.g. rat 1/2 + 1/3)".into(),
        ));
    }
    let a = crate::rational::parse_rational(tokens[0])?;
    let b = crate::rational::parse_rational(tokens[2])?;
    let op = tokens[1];
    let result = match op {
        "+" => a + b,
        "-" => a - b,
        "*" => a * b,
        "/" => {
            if b.num() == 0 {
                return Err(crate::error::MathError::Eval(
                    "`rat`: division by zero".into(),
                ));
            }
            a / b
        }
        _ => {
            return Err(crate::error::MathError::Eval(
                format!("`rat`: unknown operator '{}', use + - * /", op),
            ));
        }
    };
    Ok(vec![
        format!("{} {} {} = {}", a, op, b, result),
    ])
}

fn laurent_steps(rest: &str) -> Result<Vec<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`laurent` needs: <expr> [center] [pole_order] [n_positive]".into(),
        ));
    }
    let mut n_positive = 5usize;
    let mut pole_order = 1usize;
    let mut center = 0.0f64;
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        if let Ok(n) = tokens[expr_end - 1].parse::<usize>() {
            n_positive = n;
            expr_end -= 1;
        }
    }
    if expr_end > 1 {
        if let Ok(k) = tokens[expr_end - 1].parse::<usize>() {
            pole_order = k;
            expr_end -= 1;
        }
    }
    if expr_end > 1 {
        if let Ok(c) = tokens[expr_end - 1].parse::<f64>() {
            center = c;
            expr_end -= 1;
        }
    }
    if expr_end == 0 {
        return Err(crate::error::MathError::Eval(
            "`laurent` needs: <expr> [center] [pole_order] [n_positive]".into(),
        ));
    }
    let expr_src = tokens[..expr_end].join(" ");
    let ls = crate::laurent::laurent_series_str(&expr_src, "x", center, pole_order, n_positive)?;
    Ok(vec![
        format!("f(x) = {}", expr_src),
        format!("Laurent expansion around a = {}", format_value(center)),
        format!("pole order = {}, positive terms = {}", pole_order, n_positive),
        format!("f(x) = {}", ls.to_string()),
    ])
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
    if let Some(rest) = line.strip_prefix("pdiff ") {
        return do_pdiff(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("gradient ") {
        return do_gradient(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("idiff ") {
        return do_idiff(rest.trim());
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
    if let Some(rest) = line.strip_prefix("laurent ") {
        return do_laurent(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("rat ") {
        return do_rat(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("fourier ") {
        return do_fourier(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("mc ") {
        return do_mc(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("sample ") {
        return do_sample(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("dist ") {
        return do_dist(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("qtile ") {
        return do_qtile(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("fit ") {
        return do_fit(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("ttest ") {
        return do_ttest(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("chitest ") {
        return do_chitest(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("anova ") {
        return do_anova(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("mwu ") {
        return do_mwu(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("wilcoxon ") {
        return do_wilcoxon(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("kw ") {
        return do_kw(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("boot ") {
        return do_boot(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("logit ") {
        return do_logit(rest.trim());
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
    if let Some(rest) = line.strip_prefix("qr ") {
        return do_qr(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("tikhonov ") {
        return do_tikhonov(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("spcg ") {
        return do_spcg(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("spbicg ") {
        return do_spbicg(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("rank ") {
        return do_rank(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("cond ") {
        return do_cond(rest.trim());
    }
    if line == "null" || line.starts_with("null ") {
        let rest = line.strip_prefix("null").unwrap_or("").trim();
        return do_null(rest);
    }
    if let Some(rest) = line.strip_prefix("spline ") {
        return do_spline(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("pchip ") {
        return do_pchip(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("bspline ") {
        return do_bspline(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("hermite ") {
        return do_hermite(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("minimize ") {
        return do_minimize(rest.trim(), ctx);
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
    if let Some(rest) = line.strip_prefix("det ") {
        return do_det(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("fast ") {
        return do_fast(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("big ") {
        return do_big(rest.trim());
    }
    if line == "dec" {
        return Err(crate::error::MathError::Eval(
            "dec needs: <expr> [prec <n>] [with <var>=<val>,...]".into(),
        ));
    }
    if let Some(rest) = line.strip_prefix("dec ") {
        return do_dec(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("limit ") {
        return do_limit(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("expand ") {
        let e = Parser::parse(rest.trim())?;
        return Ok(Some(crate::poly::expand(&e).to_string()));
    }
    if let Some(rest) = line.strip_prefix("apart ") {
        return do_apart(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("ad ") {
        return do_ad(rest.trim(), &ctx.clone());
    }
    if let Some(rest) = line.strip_prefix("mathml ") {
        return do_mathml(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("serialize ") {
        return do_serialize(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("interval ") {
        return do_interval(rest.trim());
    }
    if line == "cval" {
        return Err(crate::error::MathError::Eval(
            "cval needs: <expr> [with <var>=<val>,...]".into(),
        ));
    }
    if let Some(rest) = line.strip_prefix("cval ") {
        return do_cval(rest.trim());
    }
    if line == "qsolve" {
        return Err(crate::error::MathError::Eval(
            "qsolve needs: <expr> [= rhs]  (linear/quadratic, one variable)".into(),
        ));
    }
    if let Some(rest) = line.strip_prefix("qsolve ") {
        let mut steps = crate::qsolve::qsolve_steps_str(rest.trim())?;
        return Ok(Some(steps.pop().unwrap_or_default()));
    }
    if line == "sum" || line == "prod" {
        return Err(crate::error::MathError::Eval(
            "sum/prod need: <expr> [<var>] <a> <b>".into(),
        ));
    }
    if let Some(rest) = line.strip_prefix("sum ") {
        return do_sumprod(rest.trim(), true);
    }
    if let Some(rest) = line.strip_prefix("prod ") {
        return do_sumprod(rest.trim(), false);
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

fn do_pdiff(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval("`pdiff` needs: <expr> <var>".into()));
    }
    let wrt = tokens[tokens.len() - 1];
    let expr_src = tokens[..tokens.len() - 1].join(" ");
    let e = Parser::parse(&expr_src)?;
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

fn do_gradient(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    if rest.is_empty() {
        return Err(crate::error::MathError::Eval("`gradient` needs an expression".into()));
    }
    let e = Parser::parse(rest)?;
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
    let grad = crate::symbolic::gradient(&e)?;
    if grad.is_empty() {
        return Ok(Some("(constant — no variables)".into()));
    }
    let lines: Vec<String> = grad.iter().map(|(v, d)| format!("d/d{} = {}", v, simplify(d))).collect();
    Ok(Some(lines.join("\n")))
}

fn do_idiff(rest: &str) -> Result<Option<String>> {
    let steps = crate::symbolic::idiff_steps_str(rest)?;
    // Final step: `dy/dx = -(∂F/∂x)/(∂F/∂y) = <expr>` — return the result expr.
    let last = steps.last().cloned().unwrap_or_default();
    Ok(Some(
        last.rsplit_once(" = ")
            .map(|(_, r)| r.trim().to_string())
            .unwrap_or(last),
    ))
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

fn do_laurent(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`laurent` needs: <expr> [center] [pole_order] [n_positive]".into(),
        ));
    }
    // Parse from the right: optional n_positive, optional pole_order, optional center
    let mut n_positive = 5usize;
    let mut pole_order = 1usize;
    let mut center = 0.0f64;
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        if let Ok(n) = tokens[expr_end - 1].parse::<usize>() {
            n_positive = n;
            expr_end -= 1;
        }
    }
    if expr_end > 1 {
        if let Ok(k) = tokens[expr_end - 1].parse::<usize>() {
            pole_order = k;
            expr_end -= 1;
        }
    }
    if expr_end > 1 {
        if let Ok(c) = tokens[expr_end - 1].parse::<f64>() {
            center = c;
            expr_end -= 1;
        }
    }
    if expr_end == 0 {
        return Err(crate::error::MathError::Eval(
            "`laurent` needs: <expr> [center] [pole_order] [n_positive]".into(),
        ));
    }
    let expr_src = tokens[..expr_end].join(" ");
    let ls = crate::laurent::laurent_series_str(&expr_src, "x", center, pole_order, n_positive)?;
    Ok(Some(ls.to_string()))
}

fn do_rat(rest: &str) -> Result<Option<String>> {
    // Format: rat <a> <op> <b>
    // where a, b are rationals (integers, "n/d", or decimals) and op is +, -, *, /
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`rat` needs: <a> <op> <b>  (e.g. rat 1/2 + 1/3)".into(),
        ));
    }
    let a = crate::rational::parse_rational(tokens[0])?;
    let b = crate::rational::parse_rational(tokens[2])?;
    let result = match tokens[1] {
        "+" => a + b,
        "-" => a - b,
        "*" => a * b,
        "/" => {
            if b.num() == 0 {
                return Err(crate::error::MathError::Eval(
                    "`rat`: division by zero".into(),
                ));
            }
            a / b
        }
        _ => {
            return Err(crate::error::MathError::Eval(
                format!("`rat`: unknown operator '{}', use + - * /", tokens[1]),
            ));
        }
    };
    Ok(Some(result.to_string()))
}

fn do_fourier(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`fourier` needs: <expr> <L> <n_terms> [x_eval]".into(),
        ));
    }
    // Parse from the right: optional x_eval, then n_terms, then L
    let mut x_eval: Option<f64> = None;
    let mut expr_end = tokens.len();
    if let Ok(xv) = tokens[expr_end - 1].parse::<f64>() {
        x_eval = Some(xv);
        expr_end -= 1;
    }
    if expr_end < 3 {
        return Err(crate::error::MathError::Eval(
            "`fourier` needs: <expr> <L> <n_terms> [x_eval]".into(),
        ));
    }
    let n_terms: usize = tokens[expr_end - 1].parse::<usize>().map_err(|_| {
        crate::error::MathError::Eval("n_terms must be a positive integer".into())
    })?;
    expr_end -= 1;
    let l: f64 = tokens[expr_end - 1].parse::<f64>().map_err(|_| {
        crate::error::MathError::Eval("L must be a number".into())
    })?;
    expr_end -= 1;
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    let bound: Vec<(String, Expr)> = ctx
        .vars
        .iter()
        .filter(|(k, _)| k.as_str() != "x" && !["pi", "e", "tau", "inf"].contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), Expr::num(*v)))
        .collect();
    let mut e = e;
    for (k, v) in &bound {
        e = e.substitute(k, v);
    }
    let eval_fn = |x: f64| -> f64 {
        let mut local_ctx = ctx.clone();
        local_ctx.set("x", x);
        crate::eval::eval(&e, &local_ctx).unwrap_or(0.0)
    };
    let fs = crate::calculus::fourier_series(eval_fn, n_terms, l)?;
    let mut lines = Vec::new();
    lines.push(format!("a0 = {}", format_value(fs.a0)));
    for (i, (a, b)) in fs.an.iter().zip(fs.bn.iter()).enumerate() {
        lines.push(format!(
            "a{} = {}  b{} = {}",
            i + 1,
            format_value(*a),
            i + 1,
            format_value(*b)
        ));
    }
    if let Some(xv) = x_eval {
        let val = crate::calculus::fourier_eval(&fs, xv);
        lines.push(format!("f({}) ≈ {}", format_value(xv), format_value(val)));
    }
    Ok(Some(lines.join("\n")))
}

fn do_mc(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`mc` needs: <expr> <a> <b> <n_samples> [seed]".into(),
        ));
    }
    // Parse from the right: optional seed, then n_samples, then b, a, rest is expr
    let mut seed: u64 = 42;
    let mut expr_end = tokens.len();
    if expr_end > 4 {
        if let Ok(s) = tokens[expr_end - 1].parse::<u64>() {
            seed = s;
            expr_end -= 1;
        }
    }
    if expr_end < 4 {
        return Err(crate::error::MathError::Eval(
            "`mc` needs: <expr> <a> <b> <n_samples> [seed]".into(),
        ));
    }
    let n_samples: usize = tokens[expr_end - 1].parse::<usize>().map_err(|_| {
        crate::error::MathError::Eval("n_samples must be a positive integer".into())
    })?;
    let b: f64 = tokens[expr_end - 2].parse::<f64>()
        .map_err(|_| crate::error::MathError::Eval("b must be a number".into()))?;
    let a: f64 = tokens[expr_end - 3].parse::<f64>()
        .map_err(|_| crate::error::MathError::Eval("a must be a number".into()))?;
    expr_end -= 3;
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    let bound: Vec<(String, Expr)> = ctx
        .vars
        .iter()
        .filter(|(k, _)| k.as_str() != "x" && !["pi", "e", "tau", "inf"].contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), Expr::num(*v)))
        .collect();
    let mut e = e;
    for (k, v) in &bound {
        e = e.substitute(k, v);
    }
    let eval_fn = |x: f64| -> f64 {
        let mut local_ctx = ctx.clone();
        local_ctx.set("x", x);
        crate::eval::eval(&e, &local_ctx).unwrap_or(0.0)
    };
    let (est, se) = crate::calculus::monte_carlo_integrate_1d(eval_fn, a, b, n_samples, seed)?;
    Ok(Some(format!(
        "estimate = {}\nstd_error = {}",
        format_value(est),
        format_value(se)
    )))
}

fn do_sample(rest: &str) -> Result<Option<String>> {
    use crate::stats::Rng;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`sample` needs: <dist> <params...> <n> [seed]".into(),
        ));
    }
    let dist = tokens[0];
    // Parse from the right: optional seed, then n
    let mut seed: u64 = 42;
    let mut end = tokens.len();
    if end > 2 {
        if let Ok(s) = tokens[end - 1].parse::<u64>() {
            seed = s;
            end -= 1;
        }
    }
    if end < 2 {
        return Err(crate::error::MathError::Eval(
            "`sample` needs: <dist> <params...> <n> [seed]".into(),
        ));
    }
    let n: usize = tokens[end - 1].parse::<usize>().map_err(|_| {
        crate::error::MathError::Eval("n must be a positive integer".into())
    })?;
    let params = &tokens[1..end - 1];
    let mut rng = Rng::new(seed);
    let samples: Vec<f64> = match dist {
        "uniform" => {
            if params.len() != 2 {
                return Err(crate::error::MathError::Eval("uniform needs: lo hi".into()));
            }
            let lo: f64 = params[0].parse().map_err(|_| crate::error::MathError::Eval("lo must be a number".into()))?;
            let hi: f64 = params[1].parse().map_err(|_| crate::error::MathError::Eval("hi must be a number".into()))?;
            (0..n).map(|_| rng.uniform(lo, hi)).collect()
        }
        "normal" => {
            if params.len() != 2 {
                return Err(crate::error::MathError::Eval("normal needs: mean sigma".into()));
            }
            let mu: f64 = params[0].parse().map_err(|_| crate::error::MathError::Eval("mean must be a number".into()))?;
            let sigma: f64 = params[1].parse().map_err(|_| crate::error::MathError::Eval("sigma must be a number".into()))?;
            (0..n).map(|_| rng.normal(mu, sigma)).collect()
        }
        "exponential" | "exp" => {
            if params.len() != 1 {
                return Err(crate::error::MathError::Eval("exponential needs: lambda".into()));
            }
            let lambda: f64 = params[0].parse().map_err(|_| crate::error::MathError::Eval("lambda must be a number".into()))?;
            (0..n).map(|_| rng.exponential(lambda)).collect()
        }
        _ => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown distribution '{}': use uniform, normal, or exponential",
                dist
            )));
        }
    };
    let s = crate::stats::summary(&samples)?;
    Ok(Some(format!(
        "n={}\nmean={}\nstddev={}\nmin={}\nmax={}",
        s.count, format_value(s.mean), format_value(s.stddev),
        format_value(s.min), format_value(s.max)
    )))
}

fn do_dist(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`dist` needs: <dist> <x> <params...>".into(),
        ));
    }
    let dist = tokens[0];
    let x: f64 = tokens[1].parse().map_err(|_| {
        crate::error::MathError::Eval("x must be a number".into())
    })?;
    let params = &tokens[2..];
    let need = |n: usize| -> Result<Vec<f64>> {
        if params.len() != n {
            return Err(crate::error::MathError::Eval(format!(
                "`{dist}` expects {n} parameter(s)"
            )));
        }
        params
            .iter()
            .map(|p| p.parse::<f64>().map_err(|_| crate::error::MathError::Eval("parameters must be numbers".into())))
            .collect()
    };
    let (density_label, density, cdf) = match dist {
        "normal" => {
            let p = need(2)?;
            ( "pdf", crate::stats::normal_pdf(x, p[0], p[1]), crate::stats::normal_cdf(x, p[0], p[1]))
        }
        "exponential" | "exp" => {
            let p = need(1)?;
            ("pdf", crate::stats::exp_pdf(x, p[0]), crate::stats::exp_cdf(x, p[0]))
        }
        "uniform" => {
            let p = need(2)?;
            ("pdf", crate::dists::uniform_pdf(x, p[0], p[1]), crate::dists::uniform_cdf(x, p[0], p[1]))
        }
        "t" | "student" => {
            let p = need(1)?;
            ("pdf", crate::dists::student_t_pdf(x, p[0]), crate::dists::student_t_cdf(x, p[0]))
        }
        "chi2" | "chisq" => {
            let p = need(1)?;
            ("pdf", crate::dists::chi2_pdf(x, p[0]), crate::dists::chi2_cdf(x, p[0]))
        }
        "f" => {
            let p = need(2)?;
            ("pdf", crate::dists::f_pdf(x, p[0], p[1]), crate::dists::f_cdf(x, p[0], p[1]))
        }
        "binom" | "binomial" => {
            let p = need(2)?;
            let (n, k) = (p[0], x);
            if n.fract() != 0.0 || n < 0.0 || k.fract() != 0.0 || k < 0.0 || k > n {
                return Err(crate::error::MathError::Eval("binomial needs integer 0 <= k <= n".into()));
            }
            (
                "pmf",
                crate::dists::binomial_pmf(k as u64, n as u64, p[1]),
                crate::dists::binomial_cdf(k as u64, n as u64, p[1]),
            )
        }
        "poisson" => {
            let p = need(1)?;
            let k = x;
            if k.fract() != 0.0 || k < 0.0 {
                return Err(crate::error::MathError::Eval("poisson k must be a non-negative integer".into()));
            }
            ("pmf", crate::dists::poisson_pmf(k as u64, p[0]), crate::dists::poisson_cdf(k as u64, p[0]))
        }
        _ => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown distribution '{}': use normal, exp, uniform, t, chi2, f, binom, or poisson",
                dist
            )));
        }
    };
    Ok(Some(format!(
        "{density_label} = {}\ncdf = {}",
        format_value(density),
        format_value(cdf)
    )))
}

fn do_fit(rest: &str) -> Result<Option<String>> {
    // Syntax: fit <expr> [in <var>] [with p=v, ...] x1 y1 x2 y2 ...
    // Data is always the trailing run of numeric tokens (x y pairs).
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    let mut end = tokens.len();
    while end > 0 && tokens[end - 1].parse::<f64>().is_ok() {
        end -= 1;
    }
    if end == tokens.len() || !(tokens.len() - end).is_multiple_of(2) {
        return Err(crate::error::MathError::Eval(
            "`fit` needs x y data pairs at the end: fit <expr> [in <var>] [with p=v,...] x1 y1 x2 y2 ...".into(),
        ));
    }
    let mut data = Vec::with_capacity((tokens.len() - end) / 2);
    let mut k = end;
    while k < tokens.len() {
        let x: f64 = tokens[k].parse().map_err(|_| crate::error::MathError::Eval("data must be numbers".into()))?;
        let y: f64 = tokens[k + 1].parse().map_err(|_| crate::error::MathError::Eval("data must be numbers".into()))?;
        data.push((x, y));
        k += 2;
    }
    let head = tokens[..end].join(" ");
    let (head, assignments) = match head.rfind(" with ") {
        Some(pos) => (head[..pos].to_string(), head[pos + 6..].to_string()),
        None => (head, String::new()),
    };
    let (expr_src, var) = match head.rfind(" in ") {
        Some(pos) => (
            head[..pos].trim().to_string(),
            head[pos + 4..].trim().to_string(),
        ),
        None => (head.trim().to_string(), "x".to_string()),
    };
    let expr = Parser::parse(&expr_src)?;

    // Parameters: explicit `with p=v, ...` guesses, or inferred as every
    // free variable except the independent one (initial value 1).
    let (names, init): (Vec<String>, Vec<f64>) = if assignments.is_empty() {
        let mut names: Vec<String> = expr
            .variables()
            .into_iter()
            .filter(|v| *v != var)
            .collect();
        names.sort();
        if names.is_empty() {
            return Err(crate::error::MathError::Eval(
                "no free parameters to fit; expression must contain at least one variable besides the data variable".into(),
            ));
        }
        let init = vec![1.0; names.len()];
        (names, init)
    } else {
        let mut names = Vec::new();
        let mut init = Vec::new();
        for a in split_assignments(&assignments) {
            let (name, val) = a.split_once('=').ok_or_else(|| {
                crate::error::MathError::Eval(format!("`with` entry '{a}' must be name=value"))
            })?;
            names.push(name.trim().to_string());
            init.push(val.trim().parse().map_err(|_| {
                crate::error::MathError::Eval(format!("initial guess '{val}' is not a number"))
            })?);
        }
        if names.is_empty() {
            return Err(crate::error::MathError::Eval("`with` needs at least one parameter".into()));
        }
        (names, init)
    };

    let mut mctx = Context::standard();
    let var_c = var.clone();
    let names_c = names.clone();
    let expr_c = expr.clone();
    let model = move |x: f64, p: &[f64]| -> Result<f64> {
        mctx.set(var_c.clone(), x);
        for (name, v) in names_c.iter().zip(p.iter()) {
            mctx.set(name.clone(), *v);
        }
        eval(&expr_c, &mctx)
    };
    let fit = crate::curvefit::curve_fit(model, &data, &init, &crate::curvefit::LmOptions::default())?;

    let mut out = String::new();
    for (name, (v, se)) in names
        .iter()
        .zip(fit.params.iter().zip(fit.std_errors.iter()))
    {
        out.push_str(&format!("{name} = {} ± {}\n", format_value(*v), format_value(*se)));
    }
    for (name, v) in names.iter().zip(fit.params.iter()).skip(fit.std_errors.len()) {
        out.push_str(&format!("{name} = {}\n", format_value(*v)));
    }
    out.push_str(&format!(
        "sse = {}\niterations = {} ({})",
        format_value(fit.sse),
        fit.iterations,
        if fit.converged { "converged" } else { "NOT converged" }
    ));
    Ok(Some(out))
}

fn do_mwu(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() != 2 {
        return Err(crate::error::MathError::Eval(
            "`mwu` needs: mwu <group1...> | <group2...>".into(),
        ));
    }
    let (a, b) = (parse_group(groups[0])?, parse_group(groups[1])?);
    let r = crate::dists::mann_whitney_u(&a, &b)?;
    Ok(Some(format_test_result(&r)))
}

fn do_wilcoxon(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() != 2 {
        return Err(crate::error::MathError::Eval(
            "`wilcoxon` needs: wilcoxon <before...> | <after...> (paired)".into(),
        ));
    }
    let (a, b) = (parse_group(groups[0])?, parse_group(groups[1])?);
    let r = crate::dists::wilcoxon_signed_rank(&a, &b)?;
    Ok(Some(format_test_result(&r)))
}

fn do_kw(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`kw` needs: kw <group1...> | <group2...> [| ...]".into(),
        ));
    }
    let parsed: Vec<Vec<f64>> = groups.iter().map(|g| parse_group(g)).collect::<Result<_>>()?;
    let refs: Vec<&[f64]> = parsed.iter().map(|g| g.as_slice()).collect();
    let r = crate::dists::kruskal_wallis(&refs)?;
    Ok(Some(format_test_result(&r)))
}

fn do_boot(rest: &str) -> Result<Option<String>> {
    use crate::stats::{mean, median};
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`boot` needs: boot <mean|median> <values...> [seed <n>]".into(),
        ));
    }
    // optional trailing `seed <n>` keyword (a bare trailing number is data)
    let mut seed = 42u64;
    let mut end = tokens.len();
    if end > 3 && tokens[end - 2] == "seed" {
        seed = tokens[end - 1].parse::<u64>().map_err(|_| {
            crate::error::MathError::Eval("seed must be a positive integer".into())
        })?;
        end -= 2;
    }
    let values = parse_floats(&tokens[1..end])?;
    const ITERS: usize = 10_000;
    const CONF: f64 = 0.95;
    let (stat_name, lo, hi, point) = match tokens[0] {
        "mean" => {
            let (lo, hi) = crate::dists::bootstrap_ci(&values, |s| mean(s).unwrap(), ITERS, CONF, seed)?;
            ("mean", lo, hi, mean(&values)?)
        }
        "median" => {
            let (lo, hi) = crate::dists::bootstrap_ci(&values, |s| median(s).unwrap(), ITERS, CONF, seed)?;
            ("median", lo, hi, median(&values)?)
        }
        other => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown statistic '{}': use mean or median",
                other
            )));
        }
    };
    Ok(Some(format!(
        "{stat_name} = {}\n95% ci = [{}, {}]\niters = {} (seed {})",
        format_value(point),
        format_value(lo),
        format_value(hi),
        ITERS,
        seed
    )))
}

fn do_logit(rest: &str) -> Result<Option<String>> {
    let (resp_str, cols_str) = rest.split_once(" with ").ok_or_else(|| {
        crate::error::MathError::Eval(
            "`logit` needs: logit <y0 1...> with <x1...> [| <x2...>]".into(),
        )
    })?;
    let y_tokens: Vec<&str> = resp_str.split_whitespace().collect();
    let y = parse_floats(&y_tokens)?;
    let cols: Vec<Vec<f64>> = cols_str.split('|').map(parse_group).collect::<Result<_>>()?;
    if cols.iter().any(|c| c.len() != y.len()) {
        return Err(crate::error::MathError::Eval(
            "each predictor column must have the same length as the response".into(),
        ));
    }
    let refs: Vec<&[f64]> = cols.iter().map(|c| c.as_slice()).collect();
    let fit =
        crate::logit::logistic_regression(&refs, &y, &crate::logit::LogitOptions::default())?;
    let mut out = String::new();
    out.push_str(&format!(
        "intercept = {} ± {}\n",
        format_value(fit.coefficients[0]),
        format_value(fit.std_errors[0])
    ));
    for (j, (b, se)) in fit.coefficients[1..]
        .iter()
        .zip(fit.std_errors[1..].iter())
        .enumerate()
    {
        out.push_str(&format!("b{} = {} ± {}\n", j + 1, format_value(*b), format_value(*se)));
    }
    out.push_str(&format!(
        "log-lik = {}\niterations = {} ({})\n",
        format_value(fit.log_likelihood),
        fit.iterations,
        if fit.converged { "converged" } else { "NOT converged" }
    ));
    let probs: Vec<String> = (0..y.len())
        .map(|i| {
            let features: Vec<f64> = cols.iter().map(|c| c[i]).collect();
            format_value(crate::logit::predict_proba(&fit.coefficients, &features).unwrap())
        })
        .collect();
    let shown: Vec<String> = probs.iter().take(12).cloned().collect();
    out.push_str(&format!("p = {}", shown.join(", ")));
    if probs.len() > 12 {
        out.push_str(", …");
    }
    Ok(Some(out))
}

fn do_qtile(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`qtile` needs: qtile <dist> <p> <params...>".into(),
        ));
    }
    let dist = tokens[0];
    let p: f64 = tokens[1].parse().map_err(|_| {
        crate::error::MathError::Eval("p must be a probability in (0, 1)".into())
    })?;
    let params = &tokens[2..];
    let need = |n: usize| -> Result<Vec<f64>> {
        if params.len() != n {
            return Err(crate::error::MathError::Eval(format!(
                "`{dist}` expects {n} parameter(s)"
            )));
        }
        params
            .iter()
            .map(|t| t.parse::<f64>().map_err(|_| crate::error::MathError::Eval("parameters must be numbers".into())))
            .collect()
    };
    let x = match dist {
        "normal" => {
            // optional mu sigma (default 0, 1)
            match params.len() {
                0 => crate::dists::normal_ppf(p),
                2 => {
                    let ps = need(2)?;
                    crate::dists::normal_ppf(p) * ps[1] + ps[0]
                }
                _ => return Err(crate::error::MathError::Eval("`normal` expects 0 or 2 parameters ([mu] [sigma])".into())),
            }
        }
        "uniform" => {
            let ps = need(2)?;
            crate::dists::uniform_ppf(p, ps[0], ps[1])
        }
        "t" | "student" => {
            let ps = need(1)?;
            crate::dists::student_t_ppf(p, ps[0])
        }
        "chi2" | "chisq" => {
            let ps = need(1)?;
            crate::dists::chi2_ppf(p, ps[0])
        }
        "f" => {
            let ps = need(2)?;
            crate::dists::f_ppf(p, ps[0], ps[1])
        }
        "binom" | "binomial" | "poisson" => {
            return Err(crate::error::MathError::Eval(
                "discrete distributions have no continuous quantile; use `dist` for pmf/cdf".into(),
            ));
        }
        _ => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown distribution '{}': use normal, uniform, t, chi2, or f",
                dist
            )));
        }
    };
    Ok(Some(format!("q = {}", format_value(x))))
}

fn do_ttest(rest: &str) -> Result<Option<String>> {
    let (paired, rest) = match rest.strip_prefix("paired ") {
        Some(r) => (true, r),
        None => (false, rest),
    };
    let groups = split_groups(rest)?;
    let r = match (paired, groups.len()) {
        (false, 1) => {
            let tokens: Vec<&str> = groups[0].split_whitespace().collect();
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval(
                    "`ttest` needs: ttest <mu> <values...>  or  ttest <values...> | <values...>".into(),
                ));
            }
            let mu: f64 = tokens[0].parse().map_err(|_| {
                crate::error::MathError::Eval(
                    "first token must be the tested mean mu (or use `|` for two-sample)".into(),
                )
            })?;
            let values = parse_floats(&tokens[1..])?;
            crate::dists::t_test_one(&values, mu)?
        }
        (false, 2) => {
            let (a, b) = (parse_group(groups[0])?, parse_group(groups[1])?);
            crate::dists::t_test_two(&a, &b)?
        }
        (true, 2) => {
            let (a, b) = (parse_group(groups[0])?, parse_group(groups[1])?);
            crate::dists::t_test_paired(&a, &b)?
        }
        _ => {
            return Err(crate::error::MathError::Eval(
                "`ttest` needs 1 group (with mu) or 2 groups separated by '|'".into(),
            ));
        }
    };
    Ok(Some(format_test_result(&r)))
}

fn do_chitest(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    let r = match groups.len() {
        1 => {
            let obs = parse_group(groups[0])?;
            crate::dists::chi_square_uniform(&obs)?
        }
        2 => {
            let obs = parse_group(groups[0])?;
            let exp = parse_group(groups[1])?;
            crate::dists::chi_square_gof(&obs, &exp)?
        }
        _ => {
            return Err(crate::error::MathError::Eval(
                "`chitest` needs: chitest <observed...> [| <expected...>]".into(),
            ));
        }
    };
    Ok(Some(format_test_result(&r)))
}

fn do_anova(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`anova` needs: anova <group1...> | <group2...> [| <group3...> ...]".into(),
        ));
    }
    let parsed: Vec<Vec<f64>> = groups.iter().map(|g| parse_group(g)).collect::<Result<_>>()?;
    let refs: Vec<&[f64]> = parsed.iter().map(|g| g.as_slice()).collect();
    let r = crate::dists::anova_oneway(&refs)?;
    Ok(Some(format_test_result(&r)))
}

/// Split a REPL argument string into groups on unquoted `|` separators.
fn split_groups(s: &str) -> Result<Vec<&str>> {
    let groups: Vec<&str> = s.split('|').collect();
    if groups.len() < 2 {
        return Ok(vec![s]);
    }
    Ok(groups)
}

fn parse_group(s: &str) -> Result<Vec<f64>> {
    let tokens: Vec<&str> = s.split_whitespace().collect();
    parse_floats(&tokens)
}

fn parse_floats(tokens: &[&str]) -> Result<Vec<f64>> {
    tokens
        .iter()
        .map(|t| {
            t.parse::<f64>()
                .map_err(|_| crate::error::MathError::Eval(format!("'{}' is not a number", t)))
        })
        .collect()
}

fn format_test_result(r: &crate::dists::TestResult) -> String {
    let df_str = match (r.df, r.df2) {
        (Some(d1), Some(d2)) => format!("df=({}, {})", format_value(d1), format_value(d2)),
        (Some(d1), None) => format!("df={}", format_value(d1)),
        (None, _) => String::new(),
    };
    let mut out = format!("{}: stat={}\n", r.name, format_value(r.statistic));
    if !df_str.is_empty() {
        out.push_str(&df_str);
        out.push('\n');
    }
    out.push_str(&format!("p={}", format_value(r.p_value)));
    out
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
                // Group repeated factors with exponents: 360 = 2^3 * 3^2 * 5
                let mut grouped: Vec<(u64, u32)> = Vec::new();
                for &p in &factors {
                    if let Some(last) = grouped.last_mut() {
                        if last.0 == p { last.1 += 1; continue; }
                    }
                    grouped.push((p, 1));
                }
                let strs: Vec<String> = grouped
                    .iter()
                    .map(|(p, e)| if *e == 1 { p.to_string() } else { format!("{}^{}", p, e) })
                    .collect();
                Ok(Some(strs.join(" · ")))
            }
        }
        "fib" => {
            let n = parse_u64_single(&tokens)?;
            // F(90) ≈ 2.88e18 fits in u64; F(91) overflows. Auto-upgrade to BigInt.
            if n <= 90 {
                Ok(Some(numtheory::fibonacci(n).to_string()))
            } else {
                Ok(Some(crate::bigint::fibonacci(n).to_string()))
            }
        }
        "binom" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("binom needs n k".into())); }
            let n: u64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
            let k: u64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad k".into()))?;
            // Try u64 first; fall back to BigInt on overflow (matching SymPy behavior).
            match numtheory::binomial(n, k) {
                Ok(r) => Ok(Some(r.to_string())),
                Err(_) => Ok(Some(crate::bigint::binomial(n, k).to_string())),
            }
        }
        "fact" => {
            let n = parse_u64_single(&tokens)?;
            // 20! fits in u64; 21! overflows. Auto-upgrade to BigInt.
            if n <= 20 {
                let r = numtheory::factorial(n)?;
                Ok(Some(r.to_string()))
            } else {
                Ok(Some(crate::bigint::factorial(n).to_string()))
            }
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

fn do_fast(rest: &str) -> Result<Option<String>> {
    // Format: "func x" — e.g. "fast sin 1.5"
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "fast needs: <func> <x>  (func: sin, cos, tan, exp, log, sqrt, pow)".into(),
        ));
    }
    let func = tokens[0];
    let x: f64 = tokens[1]
        .parse()
        .map_err(|_| crate::error::MathError::Eval("bad x".into()))?;
    let (fast_val, exact_val, name) = match func {
        "sin" => (crate::fastmath::fast_sin(x), x.sin(), "sin"),
        "cos" => (crate::fastmath::fast_cos(x), x.cos(), "cos"),
        "tan" => (crate::fastmath::fast_tan(x), x.tan(), "tan"),
        "exp" => (crate::fastmath::fast_exp(x), x.exp(), "exp"),
        "log" | "ln" => (crate::fastmath::fast_log(x), x.ln(), "ln"),
        "sqrt" => (crate::fastmath::fast_sqrt(x), x.sqrt(), "sqrt"),
        "pow" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("fast pow needs: <x> <y>".into()));
            }
            let y: f64 = tokens[2]
                .parse()
                .map_err(|_| crate::error::MathError::Eval("bad y".into()))?;
            let fv = crate::fastmath::fast_pow(x, y);
            let ev = x.powf(y);
            return Ok(Some(format!(
                "fast pow({}, {}) = {}  (exact: {}, err: {:e})",
                x, y, format_value(fv), format_value(ev), (fv - ev).abs()
            )));
        }
        _ => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown fast func '{}' (try: sin, cos, tan, exp, log, sqrt, pow)",
                func
            )));
        }
    };
    Ok(Some(format!(
        "fast {}({}) = {}  (exact: {}, err: {:e})",
        name,
        x,
        format_value(fast_val),
        format_value(exact_val),
        (fast_val - exact_val).abs()
    )))
}

fn do_big(rest: &str) -> Result<Option<String>> {
    use crate::bigint;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "big needs: <op> <args>  (op: prime, factor, gcd, lcm, modpow, totient)".into(),
        ));
    }
    let op = tokens[0];
    match op {
        "prime" | "is-prime" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big prime needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            let result = bigint::is_prime(&n, 20);
            Ok(Some(format!("{} is {}", n, if result { "prime" } else { "composite" })))
        }
        "factor" | "factorize" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big factor needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            let factors = bigint::factorize(&n);
            if factors.is_empty() {
                return Ok(Some(format!("{} = 1", n)));
            }
            let parts: Vec<String> = factors
                .iter()
                .map(|(p, e)| {
                    if *e == 1 {
                        p.to_string()
                    } else {
                        format!("{}^{}", p, e)
                    }
                })
                .collect();
            Ok(Some(format!("{} = {}", n, parts.join(" · "))))
        }
        "gcd" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("big gcd needs: <a> <b>".into()));
            }
            let a = bigint::parse(tokens[1])?;
            let b = bigint::parse(tokens[2])?;
            Ok(Some(format!("gcd({}, {}) = {}", a, b, bigint::gcd(&a, &b))))
        }
        "lcm" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("big lcm needs: <a> <b>".into()));
            }
            let a = bigint::parse(tokens[1])?;
            let b = bigint::parse(tokens[2])?;
            Ok(Some(format!("lcm({}, {}) = {}", a, b, bigint::lcm(&a, &b))))
        }
        "modpow" => {
            if tokens.len() < 4 {
                return Err(crate::error::MathError::Eval("big modpow needs: <base> <exp> <mod>".into()));
            }
            let base = bigint::parse(tokens[1])?;
            let exp = bigint::parse(tokens[2])?;
            let m = bigint::parse(tokens[3])?;
            let result = bigint::mod_pow(&base, &exp, &m)?;
            Ok(Some(format!("{}^{} ≡ {} (mod {})", base, exp, result, m)))
        }
        "totient" | "phi" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big totient needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            Ok(Some(format!("φ({}) = {}", n, bigint::totient(&n))))
        }
        _ => Err(crate::error::MathError::Eval(format!(
            "unknown big op '{}' (try: prime, factor, gcd, lcm, modpow, totient)",
            op
        ))),
    }
}

/// `dec <expr> [prec <n>] [with <var>=<val>,...]` — arbitrary-precision
/// decimal evaluation. `prec` is the number of significant digits (default
/// 30, max 1000); `pi`, `e`, and `tau` are computed at full precision.
fn do_dec(rest: &str) -> Result<Option<String>> {
    use crate::bigdec;
    const DEFAULT_PREC: usize = 30;
    let mut prec = DEFAULT_PREC;
    let mut work = rest.trim().to_string();
    if let Some(idx) = work.rfind(" prec ") {
        let tail = work[idx + 6..].trim();
        let n: usize = tail.parse().map_err(|_| {
            crate::error::MathError::Eval("dec: prec must be a positive integer".into())
        })?;
        if n == 0 || n > 1000 {
            return Err(crate::error::MathError::Eval(
                "dec: prec must be between 1 and 1000".into(),
            ));
        }
        prec = n;
        work.truncate(idx);
    }
    let mut assign_part: Option<String> = None;
    if let Some(idx) = work.find(" with ") {
        assign_part = Some(work[idx + 6..].to_string());
        work.truncate(idx);
    }
    let expr_src = work.trim();
    if expr_src.is_empty() {
        return Err(crate::error::MathError::Eval(
            "dec needs: <expr> [prec <n>] [with <var>=<val>,...]".into(),
        ));
    }
    let e = Parser::parse(expr_src)?;
    let mut vars: std::collections::HashMap<String, bigdecimal::BigDecimal> =
        std::collections::HashMap::new();
    if let Some(a) = assign_part {
        for assignment in split_assignments(&a) {
            let kv: Vec<&str> = assignment.splitn(2, '=').collect();
            if kv.len() != 2 {
                return Err(crate::error::MathError::Eval(format!(
                    "bad assignment `{}` (expected: var=value)",
                    assignment
                )));
            }
            let name = kv[0].trim().to_string();
            let val = bigdec::parse(kv[1].trim())?;
            vars.insert(name, val);
        }
    }
    let result = bigdec::eval_decimal_rounded(&e, &vars, prec)?;
    Ok(Some(result.to_string()))
}

/// `limit <expr> <var> <point>` or `limit <expr> <point>` — compute the
/// limit as the variable approaches the point (finite, `inf`, or `-inf`).
/// The variable may be omitted when the expression has exactly one.
fn parse_limit_target(rest: &str) -> Result<(crate::expr::Expr, String, f64)> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "limit needs: <expr> [<var>] <point>".into(),
        ));
    }
    let parse_point = |tok: &str| -> Result<f64> {
        match tok.to_ascii_lowercase().as_str() {
            "inf" | "+inf" | "infinity" | "+infinity" => Ok(f64::INFINITY),
            "-inf" | "-infinity" => Ok(f64::NEG_INFINITY),
            _ => tok.parse::<f64>().map_err(|_| {
                crate::error::MathError::Eval(format!("invalid limit point: {}", tok))
            }),
        }
    };
    if tokens.len() >= 3 {
        let point = parse_point(tokens[tokens.len() - 1])?;
        let var = tokens[tokens.len() - 2].to_string();
        let expr_src = tokens[..tokens.len() - 2].join(" ");
        let e = Parser::parse(&expr_src)?;
        Ok((e, var, point))
    } else if tokens.len() == 2 {
        let point = parse_point(tokens[1])?;
        let e = Parser::parse(tokens[0])?;
        let mut vars = e.variables();
        if vars.len() != 1 {
            return Err(crate::error::MathError::Eval(
                "limit needs: <expr> <var> <point> (expression has multiple or no variables)"
                    .into(),
            ));
        }
        Ok((e, vars.remove(0), point))
    } else {
        Err(crate::error::MathError::Eval(
            "limit needs: <expr> [<var>] <point>  e.g. `limit sin(x)/x x 0` or `limit 1/x inf`"
                .into(),
        ))
    }
}

fn do_limit(rest: &str) -> Result<Option<String>> {
    let (e, var, point) = parse_limit_target(rest)?;
    let v = crate::limit::limit(&e, &var, point)?;
    Ok(Some(format!(
        "lim {} as {} → {} = {}",
        crate::simplify::simplify(&e),
        var,
        crate::limit::fmt_point(point),
        v
    )))
}

/// Step-by-step variant used by `dispatch_steps` (notebook UI).
fn limit_steps_str(rest: &str) -> Result<Vec<String>> {
    let (e, var, point) = parse_limit_target(rest)?;
    crate::limit::limit_steps(&e, &var, point)
}

/// `apart <expr> [var]` — partial fraction decomposition of a rational
/// function. The variable may be omitted when the expression has exactly
/// one.
fn do_apart(rest: &str) -> Result<Option<String>> {
    let (e, var) = parse_single_var_target(rest)?;
    let out = crate::apart::apart(&e, &var)?;
    Ok(Some(out.to_string()))
}

/// Parse `"<expr> [var]"` — infer the variable from the expression when
/// omitted (requires exactly one free variable).
fn parse_single_var_target(rest: &str) -> Result<(crate::expr::Expr, String)> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "needs: <expr> [<var>]".into(),
        ));
    }
    if tokens.len() >= 2 {
        let var = tokens[tokens.len() - 1].to_string();
        let expr_src = tokens[..tokens.len() - 1].join(" ");
        let e = Parser::parse(&expr_src)?;
        Ok((e, var))
    } else {
        let e = Parser::parse(tokens[0])?;
        let mut vars = e.variables();
        if vars.len() != 1 {
            return Err(crate::error::MathError::Eval(
                "needs: <expr> <var> (expression has multiple or no variables)".into(),
            ));
        }
        Ok((e, vars.remove(0)))
    }
}

/// `mathml <expr>` → export to Presentation MathML
/// `mathml import <MathML>` → import from Presentation MathML
fn do_mathml(rest: &str) -> Result<Option<String>> {
    if let Some(ml) = rest.strip_prefix("import ") {
        let e = crate::mathml::from_mathml(ml.trim())?;
        return Ok(Some(e.to_string()));
    }
    let e = Parser::parse(rest)?;
    Ok(Some(crate::mathml::to_mathml_doc(&e)))
}

/// `serialize <fmt> <expr>` → export expression to `<fmt>` (sexpr | json | rpn)
/// `serialize <fmt> import <text>` → import from `<fmt>` and show the expression
fn do_serialize(rest: &str) -> Result<Option<String>> {
    let (fmt, tail) = rest
        .split_once(char::is_whitespace)
        .map(|(f, t)| (f, t.trim()))
        .unwrap_or((rest, ""));
    if tail.is_empty() {
        return Err(crate::error::MathError::InvalidArgument(format!(
            "serialize: expected `serialize <fmt> <expr>` or `serialize <fmt> import <text>` (fmt: sexpr, json, rpn); got `{}`",
            rest
        )));
    }
    match fmt {
        "sexpr" => {
            if let Some(s) = tail.strip_prefix("import ") {
                let e = crate::serialize::from_sexpr(s.trim())?;
                Ok(Some(e.to_string()))
            } else {
                let e = Parser::parse(tail)?;
                Ok(Some(crate::serialize::to_sexpr(&e)))
            }
        }
        "json" => {
            if let Some(s) = tail.strip_prefix("import ") {
                let e = crate::serialize::from_json(s.trim())?;
                Ok(Some(e.to_string()))
            } else {
                let e = Parser::parse(tail)?;
                Ok(Some(crate::serialize::to_json(&e)))
            }
        }
        "rpn" => {
            if let Some(s) = tail.strip_prefix("import ") {
                let e = crate::serialize::from_rpn(s.trim())?;
                Ok(Some(e.to_string()))
            } else {
                let e = Parser::parse(tail)?;
                Ok(Some(crate::serialize::to_rpn(&e)))
            }
        }
        other => Err(crate::error::MathError::InvalidArgument(format!(
            "serialize: unknown format `{}` (try: sexpr, json, rpn)",
            other
        ))),
    }
}

/// `interval <expr> with <var>=[lo,hi], <var>=[lo,hi], ...`
/// → rigorous bounds on the expression over the given variable intervals.
fn do_interval(rest: &str) -> Result<Option<String>> {
    let parts: Vec<&str> = rest.splitn(2, " with ").collect();
    if parts.len() != 2 {
        return Err(crate::error::MathError::Eval(
            "interval needs: <expr> with <var>=[lo,hi],...".into(),
        ));
    }
    let expr_src = parts[0].trim();
    let assignments = parts[1].trim();
    let e = Parser::parse(expr_src)?;
    let mut vars: std::collections::HashMap<String, crate::interval::Interval> =
        std::collections::HashMap::new();
    // Split assignments on commas, but respect `[lo,hi]` brackets so the comma
    // inside an interval bound is not treated as an assignment separator.
    for assignment in split_assignments(assignments) {
        let assignment = assignment.trim();
        if assignment.is_empty() {
            continue;
        }
        // Format: var=[lo, hi]
        let kv: Vec<&str> = assignment.splitn(2, '=').collect();
        if kv.len() != 2 {
            return Err(crate::error::MathError::Eval(format!(
                "bad interval assignment `{}` (expected: var=[lo,hi])",
                assignment
            )));
        }
        let var = kv[0].trim().to_string();
        let bracketed = kv[1].trim();
        let inner = bracketed
            .strip_prefix('[')
            .and_then(|s| s.strip_suffix(']'))
            .ok_or_else(|| {
                crate::error::MathError::Eval(format!(
                    "bad interval `{}` (expected: [lo,hi])",
                    bracketed
                ))
            })?;
        let bounds: Vec<&str> = inner.split(',').map(|s| s.trim()).collect();
        if bounds.len() != 2 {
            return Err(crate::error::MathError::Eval(format!(
                "bad interval `{}` (expected: [lo,hi])",
                bracketed
            )));
        }
        let lo: f64 = bounds[0]
            .parse()
            .map_err(|_| crate::error::MathError::Eval(format!("bad lo bound: {}", bounds[0])))?;
        let hi: f64 = bounds[1]
            .parse()
            .map_err(|_| crate::error::MathError::Eval(format!("bad hi bound: {}", bounds[1])))?;
        vars.insert(var, crate::interval::Interval::new(lo, hi));
    }
    let result = crate::interval::eval_interval(&e, &vars)?;
    Ok(Some(result.to_string()))
}

fn do_cval(rest: &str) -> Result<Option<String>> {
    let z = crate::ceval::eval_complex_str(rest)?;
    Ok(Some(crate::ceval::format_complex(z)))
}

/// Σ/∏ over integer bounds: `<expr> [<var>] <a> <b>` with bounds from the end.
fn do_sumprod(rest: &str, is_sum: bool) -> Result<Option<String>> {
    let (src, var, a, b) = crate::sumprod::parse_bounds(rest)?;
    let e = Parser::parse(src.trim())?;
    let var = if var.is_empty() {
        let mut vs = e.variables();
        match vs.len() {
            1 => vs.remove(0),
            0 => String::new(), // constant term — handled by the closed form
            _ => {
                return Err(crate::error::MathError::InvalidArgument(format!(
                    "multiple variables ({}) — specify one: {} <expr> <var> {} {}",
                    vs.join(", "),
                    if is_sum { "sum" } else { "prod" },
                    a,
                    b
                )))
            }
        }
    } else {
        var
    };
    let result = if is_sum {
        crate::sumprod::summation(&e, &var, a, b)?
    } else {
        crate::sumprod::product(&e, &var, a, b)?
    };
    let out = match &result {
        Expr::Num(v) => format_value(*v),
        other => other.to_string(),
    };
    Ok(Some(out))
}

/// Split `var=[lo,hi], var=[lo,hi], ...` into individual assignments, respecting
/// `[...]` brackets so the comma inside an interval is not treated as a separator.
fn split_assignments(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut buf = String::new();
    for c in s.chars() {
        match c {
            '[' => {
                depth += 1;
                buf.push(c);
            }
            ']' => {
                depth -= 1;
                buf.push(c);
            }
            ',' if depth == 0 => {
                let trimmed = buf.trim().to_string();
                if !trimmed.is_empty() {
                    out.push(trimmed);
                }
                buf.clear();
            }
            _ => buf.push(c),
        }
    }
    let trimmed = buf.trim().to_string();
    if !trimmed.is_empty() {
        out.push(trimmed);
    }
    out
}

fn do_ad(rest: &str, ctx: &Context) -> Result<Option<String>> {
    // Format: "ad <expr> at <var>=<val>" or "ad grad <expr> with <var>=<val>,..."
    // or "ad jacobian <f1>, <f2>, ... with <var>=<val>,..."
    let tokens: Vec<&str> = rest.splitn(2, " at ").collect();
    if tokens.len() == 2 {
        // Single-variable derivative: ad <expr> at <var>=<val>
        let expr_src = tokens[0].trim();
        let var_val = tokens[1].trim();
        let parts: Vec<&str> = var_val.splitn(2, '=').collect();
        if parts.len() != 2 {
            return Err(crate::error::MathError::Eval(
                "ad needs: <expr> at <var>=<val>".into(),
            ));
        }
        let var = parts[0].trim().to_string();
        let val: f64 = parts[1]
            .trim()
            .parse()
            .map_err(|_| crate::error::MathError::Eval("bad value".into()))?;
        let e = Parser::parse(expr_src)?;
        let d = crate::autodiff::derivative(&e, &var, val, ctx)?;
        return Ok(Some(format!(
            "f({}) = {},  f'({}) = {}",
            var, format_value(d.val), var, format_value(d.deriv)
        )));
    }

    // Gradient: ad grad <expr> with x=1,y=2,...
    if let Some(rest) = rest.strip_prefix("grad ") {
        let parts: Vec<&str> = rest.splitn(2, " with ").collect();
        if parts.len() != 2 {
            return Err(crate::error::MathError::Eval(
                "ad grad needs: <expr> with <var>=<val>,...".into(),
            ));
        }
        let expr_src = parts[0].trim();
        let assignments = parts[1].trim();
        let e = Parser::parse(expr_src)?;
        let mut point = ctx.clone();
        for assignment in assignments.split(',') {
            let kv: Vec<&str> = assignment.splitn(2, '=').collect();
            if kv.len() == 2 {
                let k = kv[0].trim();
                let v: f64 = kv[1]
                    .trim()
                    .parse()
                    .map_err(|_| crate::error::MathError::Eval("bad value".into()))?;
                point.set(k, v);
            }
        }
        let grad = crate::autodiff::gradient(&e, &point)?;
        let parts: Vec<String> = grad
            .iter()
            .map(|(name, val)| format!("∂f/∂{} = {}", name, format_value(*val)))
            .collect();
        return Ok(Some(format!("∇f = [{}]", parts.join(", "))));
    }

    // Jacobian: ad jacobian <f1>, <f2>, ... with x=1,y=2,...
    if let Some(rest) = rest.strip_prefix("jacobian ") {
        let parts: Vec<&str> = rest.splitn(2, " with ").collect();
        if parts.len() != 2 {
            return Err(crate::error::MathError::Eval(
                "ad jacobian needs: <f1>, <f2>, ... with <var>=<val>,...".into(),
            ));
        }
        let exprs_src = parts[0].trim();
        let assignments = parts[1].trim();
        let exprs: Result<Vec<Expr>> = exprs_src
            .split(',')
            .map(|s| Parser::parse(s.trim()))
            .collect();
        let exprs = exprs?;
        let mut point = ctx.clone();
        for assignment in assignments.split(',') {
            let kv: Vec<&str> = assignment.splitn(2, '=').collect();
            if kv.len() == 2 {
                let k = kv[0].trim();
                let v: f64 = kv[1]
                    .trim()
                    .parse()
                    .map_err(|_| crate::error::MathError::Eval("bad value".into()))?;
                point.set(k, v);
            }
        }
        let jac = crate::autodiff::jacobian(&exprs, &point)?;
        let rows: Vec<String> = jac
            .iter()
            .map(|row| {
                let vals: Vec<String> = row.iter().map(|v| format_value(*v)).collect();
                format!("[{}]", vals.join(", "))
            })
            .collect();
        return Ok(Some(format!("J = [{}]", rows.join(", "))));
    }

    Err(crate::error::MathError::Eval(
        "ad needs: <expr> at <var>=<val>  |  ad grad <expr> with <var>=<val>,...  |  ad jacobian <f1>,... with <var>=<val>,...".into(),
    ))
}

fn do_romberg(rest: &str) -> Result<Option<String>> {
    // Format: "<expr> a b" (expr may contain spaces; fixed 8 Romberg levels)
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval("romberg needs: <expr> a b".into()));
    }
    let b: f64 = tokens[tokens.len() - 1].parse()
        .map_err(|_| crate::error::MathError::Eval("bad b".into()))?;
    let a: f64 = tokens[tokens.len() - 2].parse()
        .map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
    let levels: usize = 8;
    let expr_src = tokens[..tokens.len() - 2].join(" ");
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

/// `qr <rows...>` — QR decomposition via Householder reflections.
/// Prints Q (orthogonal) and R (upper-triangular), plus reconstruction error.
fn do_qr(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let qr = m.qr()?;
    let q = qr.q();
    let r = qr.r();
    let recon = qr.reconstruct();
    let mut max_diff = 0.0_f64;
    for i in 0..m.rows {
        for j in 0..m.cols {
            let d = (m[(i, j)] - recon[(i, j)]).abs();
            if d > max_diff {
                max_diff = d;
            }
        }
    }
    Ok(Some(format!(
        "QR ok (max reconstruction error = {:.2e})\nQ =\n{}\nR =\n{}",
        max_diff, q, r
    )))
}

fn do_tikhonov(rest: &str) -> Result<Option<String>> {
    // Format: <matrix rows separated by |> | <b vector> <lambda>
    // The last token is lambda, the second-to-last starts the b vector.
    // We split on '|' — the last group is "b1 b2 ... lambda", the rest are matrix rows.
    let parts: Vec<&str> = rest.split('|').collect();
    if parts.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`tikhonov` needs: <rows...> | <b...> <lambda>".into(),
        ));
    }
    let b_lambda: Vec<&str> = parts[parts.len() - 1].split_whitespace().collect();
    if b_lambda.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "need at least b values and lambda after last |".into(),
        ));
    }
    let lambda: f64 = b_lambda[b_lambda.len() - 1]
        .parse()
        .map_err(|_| crate::error::MathError::Eval("lambda must be a number".into()))?;
    let b: Vec<f64> = b_lambda[..b_lambda.len() - 1]
        .iter()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| crate::error::MathError::Eval("b values must be numbers".into()))
        })
        .collect::<Result<Vec<_>>>()?;
    let rows: Vec<Vec<f64>> = parts[..parts.len() - 1]
        .iter()
        .map(|s| {
            s.split_whitespace()
                .map(|x| {
                    x.parse::<f64>().map_err(|_| {
                        crate::error::MathError::Eval("matrix entries must be numbers".into())
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let m = crate::matrix::Matrix::from_rows(&rows)?;
    let x = m.solve_tikhonov(&b, lambda)?;
    let x_strs: Vec<String> = x.iter().map(|v| format_value(*v)).collect();
    Ok(Some(format!("x = [{}]", x_strs.join(", "))))
}

/// `bspline x1 y1 x2 y2 ... x_at` — clamped B-spline interpolant (cubic for
/// 4+ points, reduced degree for 2-3) evaluated at x_at.
fn do_bspline(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 5 || tokens.len().is_multiple_of(2) {
        return Err(crate::error::MathError::Eval(
            "bspline expects: x1 y1 x2 y2 ... x_at (at least 3 points)".into(),
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
    let last = *tokens.last().unwrap();
    let x_at: f64 = last
        .parse()
        .map_err(|_| crate::error::MathError::Eval(format!("bad x_at: {last}")))?;
    let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let s = crate::bspline::cubic_bspline_interp(&xs, &ys)?;
    Ok(Some(format!(
        "bspline({}) = {} (degree {})",
        x_at,
        format_value(s.eval(x_at)),
        s.degree()
    )))
}

/// `hermite x1 y1 d1 x2 y2 d2 ... x_at` — piecewise cubic Hermite with
/// given slopes evaluated at x_at.
fn do_hermite(rest: &str) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 7 || !(tokens.len() - 1).is_multiple_of(3) {
        return Err(crate::error::MathError::Eval(
            "hermite expects: x1 y1 d1 x2 y2 d2 ... x_at (at least 2 points)"
                .into(),
        ));
    }
    let n = (tokens.len() - 1) / 3;
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    let mut ds = Vec::with_capacity(n);
    for i in 0..n {
        let x: f64 = tokens[3 * i]
            .parse()
            .map_err(|_| crate::error::MathError::Eval(format!("bad x at {}", tokens[3 * i])))?;
        let y: f64 = tokens[3 * i + 1].parse().map_err(|_| {
            crate::error::MathError::Eval(format!("bad y at {}", tokens[3 * i + 1]))
        })?;
        let d: f64 = tokens[3 * i + 2].parse().map_err(|_| {
            crate::error::MathError::Eval(format!("bad slope at {}", tokens[3 * i + 2]))
        })?;
        xs.push(x);
        ys.push(y);
        ds.push(d);
    }
    let last = *tokens.last().unwrap();
    let x_at: f64 = last
        .parse()
        .map_err(|_| crate::error::MathError::Eval(format!("bad x_at: {last}")))?;
    let h = crate::interpolate::CubicHermite::new(&xs, &ys, &ds)?;
    Ok(Some(format!(
        "hermite({}) = {}",
        x_at,
        format_value(h.eval(x_at))
    )))
}

/// `spcg [jacobi] <rows> | <b...>` — solve Ax = b via conjugate gradient on
/// the sparse (CSR) representation of A. A must be symmetric positive-definite.
/// The optional `jacobi` keyword enables diagonal preconditioning.
fn do_spcg(rest: &str) -> Result<Option<String>> {
    let (jacobi, rest) = match rest.strip_prefix("jacobi ") {
        Some(tail) => (true, tail),
        None => (false, rest),
    };
    let parts: Vec<&str> = rest.split('|').collect();
    if parts.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`spcg` needs: [jacobi] <rows...> | <b...>".into(),
        ));
    }
    let b: Vec<f64> = parts[parts.len() - 1]
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| crate::error::MathError::Eval("b values must be numbers".into()))
        })
        .collect::<Result<Vec<_>>>()?;
    let rows: Vec<Vec<f64>> = parts[..parts.len() - 1]
        .iter()
        .map(|s| {
            s.split_whitespace()
                .map(|x| {
                    x.parse::<f64>().map_err(|_| {
                        crate::error::MathError::Eval("matrix entries must be numbers".into())
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let m = crate::matrix::Matrix::from_rows(&rows)?;
    let a = crate::sparse::Csr::from_dense(&m);
    let res = if jacobi {
        crate::sparse::conjugate_gradient_jacobi(&a, &b, 1e-12, 1000)?
    } else {
        crate::sparse::conjugate_gradient(&a, &b, 1e-12, 1000)?
    };
    let x_strs: Vec<String> = res.x.iter().map(|v| format_value(*v)).collect();
    Ok(Some(format!(
        "x = [{}]  ({} iterations, residual {:.2e})",
        x_strs.join(", "),
        res.iterations,
        res.residual
    )))
}

/// `spbicg [ilu] <rows> | <b...>` — solve Ax = b via Jacobi-preconditioned
/// BiCGStab on the sparse (CSR) representation of A. No symmetry requirement.
/// The optional `ilu` keyword switches to ILU(0) preconditioning.
fn do_spbicg(rest: &str) -> Result<Option<String>> {
    let (ilu, rest) = match rest.strip_prefix("ilu ") {
        Some(tail) => (true, tail),
        None => (false, rest),
    };
    let parts: Vec<&str> = rest.split('|').collect();
    if parts.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "`spbicg` needs: <rows...> | <b...>".into(),
        ));
    }
    let b: Vec<f64> = parts[parts.len() - 1]
        .split_whitespace()
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| crate::error::MathError::Eval("b values must be numbers".into()))
        })
        .collect::<Result<Vec<_>>>()?;
    let rows: Vec<Vec<f64>> = parts[..parts.len() - 1]
        .iter()
        .map(|s| {
            s.split_whitespace()
                .map(|x| {
                    x.parse::<f64>().map_err(|_| {
                        crate::error::MathError::Eval("matrix entries must be numbers".into())
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let m = crate::matrix::Matrix::from_rows(&rows)?;
    let a = crate::sparse::Csr::from_dense(&m);
    let res = if ilu {
        let fact = crate::sparse::Ilu0::factorize(&a)?;
        crate::sparse::bicgstab_ilu(&a, &b, 1e-12, 1000, &fact)?
    } else {
        crate::sparse::bicgstab(&a, &b, 1e-12, 1000)?
    };
    let x_strs: Vec<String> = res.x.iter().map(|v| format_value(*v)).collect();
    Ok(Some(format!(
        "x = [{}]  ({} iterations, residual {:.2e})",
        x_strs.join(", "),
        res.iterations,
        res.residual
    )))
}

fn do_rank(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let r = m.rank(1e-10);
    Ok(Some(format!("rank = {}", r)))
}

fn do_cond(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let k = m.condition_number()?;
    if k.is_infinite() {
        return Ok(Some("cond = inf (singular)".into()));
    }
    Ok(Some(format!("cond = {}", format_value(k))))
}

fn do_null(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let basis = m.nullspace(0.0)?;
    if basis.is_empty() {
        return Ok(Some("nullspace: {} (full column rank)".into()));
    }
    let strs: Vec<String> = basis
        .iter()
        .map(|v| {
            let parts: Vec<String> = v.iter().map(|x| format_value(*x)).collect();
            format!("({})", parts.join(", "))
        })
        .collect();
    Ok(Some(format!(
        "nullspace basis (dim {}):\n  {}",
        basis.len(),
        strs.join("\n  ")
    )))
}

fn do_det(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let d = m.determinant()?;
    Ok(Some(format!("det = {}", format_value(d))))
}

fn do_spline(rest: &str) -> Result<Option<String>> {
    // Format: "x1 y1 x2 y2 ...  [eval x_value]"
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 || tokens.len().is_multiple_of(2) {
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

fn do_pchip(rest: &str) -> Result<Option<String>> {
    // Format: "x1 y1 x2 y2 ... x_at" (same shape as `spline`)
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 5 || tokens.len().is_multiple_of(2) {
        return Err(crate::error::MathError::Eval(
            "pchip expects: x1 y1 x2 y2 ... x_at (at least 3 points)".into(),
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
    let last = *tokens.last().unwrap();
    let x_at: f64 = last
        .parse()
        .map_err(|_| crate::error::MathError::Eval(format!("bad x_at: {last}")))?;
    let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
    let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
    let p = crate::pchip::Pchip::new(&xs, &ys)?;
    Ok(Some(format!(
        "pchip({}) = {}",
        x_at,
        format_value(p.eval(x_at))
    )))
}

fn do_minimize(rest: &str, ctx: &Context) -> Result<Option<String>> {
    // Format: "minimize <expr> [var] a b" — golden-section search on [a, b]
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`minimize` needs: minimize <expr> [var] a b".into(),
        ));
    }
    let b: f64 = tokens[tokens.len() - 1]
        .parse()
        .map_err(|_| crate::error::MathError::Eval("b must be a number".into()))?;
    let a: f64 = tokens[tokens.len() - 2]
        .parse()
        .map_err(|_| crate::error::MathError::Eval("a must be a number".into()))?;
    let mut var = "x".to_string();
    let mut expr_end = tokens.len() - 2;
    if expr_end > 1 {
        let candidate = tokens[expr_end - 1];
        if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
            var = candidate.to_string();
            expr_end -= 1;
        }
    }
    let expr_src = tokens[..expr_end].join(" ");
    let e = Parser::parse(&expr_src)?;
    // bind all other free variables from the REPL context
    let mctx = ctx.clone();
    for v in e.variables() {
        if v != var && !["pi", "e", "tau", "inf"].contains(&v.as_str()) && !mctx.vars.contains_key(&v)
        {
            return Err(crate::error::MathError::UnknownVariable(v));
        }
    }
    let var_c = var.clone();
    let mut call_ctx = mctx.clone();
    let obj = move |x: f64| -> f64 {
        call_ctx.set(var_c.clone(), x);
        eval(&e, &call_ctx).unwrap_or(f64::INFINITY)
    };
    let (x_star, fx) = crate::optim::golden_section(obj, a, b, 1e-10)?;
    Ok(Some(format!(
        "{}* = {}\nf({}) = {}",
        var,
        format_value(x_star),
        format_value(x_star),
        format_value(fx)
    )))
}

fn parse_f64_list(rest: &str) -> Result<Vec<f64>> {
    rest.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.trim().is_empty())
        .map(|t| t.trim().parse::<f64>().map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", t))))
        .collect()
}

pub fn guess_var(expr: &str) -> String {
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
  <expr>              evaluate (e.g. sin(pi/4), 5!, gamma(0.5), bessel_j0(5))
  let x = <expr>      bind a variable
  fn f(x) = <expr>    define a function
  diff <expr> [var]   symbolic derivative
  pdiff <expr> <var>  partial derivative
  gradient <expr>     gradient (all partials)
  idiff <expr> [= rhs] [<indep> <dep>]
                      implicit differentiation dy/dx = -F_x/F_y
  simplify <expr>     constant-fold & simplify
  int <expr> a b      numerical integral over [a, b]
  solve <expr> [var] [guess]
  plot <expr> a b [out.png]
  taylor <expr> [a] [order]
                      Taylor series around a (default 0, order 5)
  laurent <expr> [a] [pole_order] [n_positive]
                      Laurent series around a (default: a=0, k=1, N=5)
  rat <a> <op> <b>   exact rational arithmetic (a, b: int or n/d or decimal)
  fourier <expr> L N [x]
                      Fourier series on [-L, L] with N terms (eval at x)
  mc <expr> a b N [seed]
                      Monte Carlo integral over [a, b] with N samples
  sample <dist> <params...> N [seed]
                      random sampling (uniform/normal/exponential)
   dist <dist> <x> <params...>
                       PDF and CDF (normal/exponential/uniform/t/chi2/f/binom/poisson)
   qtile <dist> <p> <params...>
                       quantile / critical value (normal [mu sigma]/uniform/t/chi2/f)
   fit <expr> [in <var>] [with p=v,...] x1 y1 ...
                       Levenberg-Marquardt curve fit; params inferred or
                       initialized via `with` (e.g. fit a*exp(-b*x) with a=1, b=1)
   ttest <mu> <values...>
                       one-sample t-test (two-sided)
   ttest <values...> | <values...>
                       two-sample Welch t-test
   ttest paired <values...> | <values...>
                       paired t-test
   chitest <observed...> [| <expected...>]
                       chi-squared test (uniform expectation by default)
   anova <group1...> | <group2...> [| ...]
                       one-way ANOVA F-test
   mwu <group1...> | <group2...>
                       Mann-Whitney U test (nonparametric two-sample)
   wilcoxon <before...> | <after...>
                       Wilcoxon signed-rank test (paired)
   kw <group1...> | <group2...> [| ...]
                       Kruskal-Wallis H test (nonparametric ANOVA)
   boot <mean|median> <values...> [seed <n>]
                       bootstrap 95% confidence interval (10000 resamples)
   logit <y0 1...> with <x1...> [| <x2...>]
                       logistic regression (IRLS) with Wald standard errors
  fft <numbers...>    magnitude spectrum of a real signal
  conv <a...> x <b...>  convolution of two signals
  stats <numbers...>  descriptive statistics
  poly-roots <coeffs...>  polynomial roots (highest degree first)
  isolate-roots <ints...>  real root isolation (VAS, integer coefficients)
  lu <rows>           matrix LU decomposition (rows separated by '|')
                      computes determinant and inverse
  qr <rows>           matrix QR decomposition (Householder reflections)
                      prints Q (orthogonal) and R (upper-triangular)
   tikhonov <rows> | <b...> <lambda>
                       Tikhonov-regularised solve (Ax≈b with L2 penalty)
   spcg [jacobi] <rows> | <b...>
                       sparse conjugate-gradient solve of Ax = b (A must be
                       symmetric positive-definite); `jacobi` enables
                       diagonal preconditioning
   spbicg [ilu] <rows> | <b...>
                       sparse BiCGStab solve of Ax = b (no symmetry
                       requirement); `ilu` enables ILU(0) preconditioning
  rank <rows>         matrix rank
  cond <rows>         condition number (2-norm, from SVD)
  null <rows>         nullspace basis of a homogeneous system A·x = 0
  det <rows>          matrix determinant
   spline x1 y1 x2 y2 ... x_at
                       natural cubic spline interpolant at x_at
   pchip x1 y1 x2 y2 ... x_at
                        monotone piecewise cubic (no overshoot) at x_at
   bspline x1 y1 x2 y2 ... x_at
                        clamped B-spline interpolant at x_at (cubic for 4+
                        points, reduced degree for 2-3)
   hermite x1 y1 d1 x2 y2 d2 ... x_at
                        piecewise cubic Hermite with given slopes at x_at
   minimize <expr> [var] a b
                       golden-section minimum of expr on [a, b]
  cholesky <rows>     Cholesky decomposition of symmetric positive-definite
  eig <rows>          dominant eigenvalue + eigenvector (power iteration)
  svd <rows>          singular value decomposition
  chebyshev n [x]     Chebyshev T_n(x), or T_n nodes on [-1, 1]
  fast <func> <x> [y] Chebyshev fast approx (sin, cos, tan, exp, log, sqrt, pow)
  big <op> <args>        Big integer ops for inputs > u64::MAX (prime, factor, gcd, lcm, modpow, totient)
  dec <expr> [prec <n>] [with <var>=<val>,...]
                         arbitrary-precision evaluation (default prec 30)
  limit <expr> [<var>] <point>
                         limit as var approaches point (finite, inf, -inf)
  expand <expr>          distribute products/powers into a polynomial
  apart <expr> [<var>]   partial fraction decomposition
  ad <expr> at <var>=<val>   Automatic differentiation (dual numbers)
  ad grad <expr> with ...    Gradient of multivariate expression
  ad jacobian <f1>,... with ...  Jacobian matrix
  mathml <expr>      export expression to Presentation MathML
  mathml import <ml> import Presentation MathML to expression
  serialize <fmt> <expr>      export expression (fmt: sexpr | json | rpn)
  serialize <fmt> import <t>  import from format and show the expression
   interval <expr> with <var>=[lo,hi],...  rigorous bounds via interval arithmetic
   cval <expr> [with <var>=<val>,...]
                          complex evaluation (i = imaginary unit)
   qsolve <expr> [= rhs]  symbolic linear/quadratic solve (one variable)
   sum <expr> [<var>] <a> <b>
                          summation over integer bounds (exact when possible)
   prod <expr> [<var>] <a> <b>
                          product over integer bounds
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
  fact <n>            factorial (also: n! in expressions)
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
           gamma, erf, erfc, sinc, bessel_j0, bessel_j1, bessel_j,
           digamma, trigamma, polygamma(m, x), harmonic(n), zeta,
           hurwitz(s, a), elliptic_k, elliptic_e, elliptic_f(phi, k),
           elliptic_e_inc(phi, k)
";

/// Suppresses an unused warning for the `Cow` import in environments where
/// rustc may otherwise complain; the type is used implicitly through Helper.
#[allow(dead_code)]
fn _cow_silence<'a>(_: Cow<'a, str>) {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dispatch_steps_diff() {
        let steps = dispatch_steps("diff x^3", Context::standard()).unwrap();
        assert!(steps.len() >= 3);
        assert!(steps[0].contains("x^3"));
        assert!(steps[steps.len() - 1].contains("3"));
    }

    #[test]
    fn dispatch_steps_solve() {
        let steps = dispatch_steps("solve x^2 - 4", Context::standard()).unwrap();
        assert!(steps.len() >= 3);
        assert!(steps[0].contains("x^2 - 4"));
        assert!(steps[2].contains("x"));
    }

    #[test]
    fn dispatch_steps_taylor() {
        let steps = dispatch_steps("taylor exp(x) 0 3", Context::standard()).unwrap();
        assert!(steps.len() >= 3);
        assert!(steps[0].contains("exp(x)"));
        assert!(steps[1].contains("a = 0"));
        assert!(steps[2].contains("order = 3"));
    }

    #[test]
    fn dispatch_steps_simplify() {
        let steps = dispatch_steps("simplify x + x", Context::standard()).unwrap();
        assert_eq!(steps.len(), 2);
        assert!(steps[0].contains("simplify"));
    }

    #[test]
    fn dispatch_steps_rat() {
        let steps = dispatch_steps("rat 1/2 + 1/3", Context::standard()).unwrap();
        assert_eq!(steps.len(), 1);
        assert!(steps[0].contains("5/6"));
    }

    #[test]
    fn romberg_basic() {
        let out = dispatch_str("romberg x^2 0 1", Context::standard()).unwrap().unwrap();
        assert!(out.contains("0.3333"), "got: {out}");
    }

    #[test]
    fn romberg_multitoken_expr() {
        // regression: multi-token expressions were truncated at the old expr_end
        let out = dispatch_str("romberg x^2 + 1 0 1", Context::standard()).unwrap().unwrap();
        assert!(out.contains("1.3333"), "got: {out}");
    }

    #[test]
    fn romberg_trailing_tokens_bound_parse() {
        // regression: bounds ALWAYS come from the last two tokens (previously a
        // 4th token silently changed which numbers were treated as a and b);
        // extra tokens fold into the expr by implicit multiplication (x^2*0 = 0)
        let out = dispatch_str("romberg x^2 0 1 12", Context::standard()).unwrap().unwrap();
        assert!(out.contains("(1, 12;"), "got: {out}");
        assert!(out.contains("= 0"), "got: {out}");
    }

    #[test]
    fn dispatch_steps_laurent() {
        let steps = dispatch_steps("laurent 1/x 0 1 3", Context::standard()).unwrap();
        assert!(steps.len() >= 3);
        assert!(steps[0].contains("1/x"));
    }

    #[test]
    fn dispatch_steps_plain_eval() {
        let steps = dispatch_steps("sin(pi/4)", Context::standard()).unwrap();
        assert_eq!(steps.len(), 1);
        assert!(steps[0].contains("0.707"));
    }

    #[test]
    fn dispatch_dec_sqrt2() {
        let mut ctx = Context::standard();
        let result = dispatch("dec sqrt(2) prec 30", &mut ctx).unwrap().unwrap();
        assert!(result.starts_with("1.4142135623730950488016887242"));
    }

    #[test]
    fn dispatch_dec_pi_50() {
        let mut ctx = Context::standard();
        let result = dispatch("dec pi prec 50", &mut ctx).unwrap().unwrap();
        assert_eq!(
            result,
            "3.1415926535897932384626433832795028841971693993751"
        );
    }

    #[test]
    fn dispatch_dec_third() {
        let mut ctx = Context::standard();
        let result = dispatch("dec 1/3 prec 10", &mut ctx).unwrap().unwrap();
        assert_eq!(result, "0.3333333333");
    }

    #[test]
    fn dispatch_dec_with_vars_and_default_prec() {
        let mut ctx = Context::standard();
        let result = dispatch("dec x*2 + 1 with x=1.5", &mut ctx).unwrap().unwrap();
        assert_eq!(result, "4.00000000000000000000000000000");
    }

    #[test]
    fn dispatch_dec_prec_out_of_range() {
        let mut ctx = Context::standard();
        assert!(dispatch("dec 1 prec 0", &mut ctx).is_err());
        assert!(dispatch("dec 1 prec 2000", &mut ctx).is_err());
    }

    #[test]
    fn dispatch_limit_removable() {
        let mut ctx = Context::standard();
        let result = dispatch("limit (x^2 - 1)/(x - 1) x 1", &mut ctx)
            .unwrap()
            .unwrap();
        assert!(result.contains("2"), "got: {}", result);
    }

    #[test]
    fn dispatch_limit_var_inferred() {
        let mut ctx = Context::standard();
        let result = dispatch("limit sin(x)/x 0", &mut ctx).unwrap().unwrap();
        assert!(result.contains("= 1"), "got: {}", result);
    }

    #[test]
    fn dispatch_limit_infinity() {
        let mut ctx = Context::standard();
        let result = dispatch("limit 1/x inf", &mut ctx).unwrap().unwrap();
        assert!(result.contains("= 0"), "got: {}", result);
        let result = dispatch("limit 1/x^2 0", &mut ctx).unwrap().unwrap();
        assert!(result.contains("+∞"), "got: {}", result);
        let result = dispatch("limit 1/x 0", &mut ctx).unwrap().unwrap();
        assert!(result.contains("does not exist"), "got: {}", result);
    }

    #[test]
    fn dispatch_limit_steps() {
        let steps = dispatch_steps("limit sin(x)/x x 0", Context::standard()).unwrap();
        assert!(steps.len() >= 3);
        assert!(steps[0].contains("limit of sin(x)/x"));
        assert!(steps.iter().any(|s| s.contains("L'Hôpital")));
        assert!(steps.last().unwrap().contains("limit = 1"));
    }

    #[test]
    fn dispatch_expand_command() {
        let mut ctx = Context::standard();
        let result = dispatch("expand (x+1)^3", &mut ctx).unwrap().unwrap();
        assert!(result.contains("x^3"), "got: {}", result);
        assert!(result.contains("3*x^2"), "got: {}", result);
        assert!(result.ends_with("1"), "got: {}", result);
    }

    #[test]
    fn dispatch_expand_multivariate() {
        let mut ctx = Context::standard();
        let result = dispatch("expand (x+y)*(x-y)", &mut ctx).unwrap().unwrap();
        assert!(result.contains("x^2 - y^2"), "got: {}", result);
    }

    #[test]
    fn dispatch_expand_steps() {
        let steps = dispatch_steps("expand (x+2)*(x+3)", Context::standard()).unwrap();
        assert_eq!(steps.len(), 2);
        assert!(steps[0].contains("expand"));
        assert!(steps[1].contains("x^2 + 5*x + 6"));
    }

    #[test]
    fn dispatch_steps_integrate() {
        let steps = dispatch_steps("integrate x^2", Context::standard()).unwrap();
        assert!(steps.len() >= 2);
        assert!(steps[0].contains("x^2"));
    }

    #[test]
    fn dispatch_steps_empty() {
        let steps = dispatch_steps("", Context::standard()).unwrap();
        assert_eq!(steps.len(), 0);
    }

    #[test]
    fn dispatch_steps_error() {
        // Parse error should propagate
        assert!(dispatch_steps("diff @@@@", Context::standard()).is_err());
    }
}