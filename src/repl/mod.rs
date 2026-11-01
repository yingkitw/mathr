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

mod calc;
mod io;
mod numtheory;
mod solve;
mod stats;

use calc::{
    do_ad, do_apart, do_diff, do_gradient, do_idiff, do_integrate, do_integrate_sym, do_laurent,
    do_limit, do_pdiff, do_rat, do_romberg, do_sumprod, do_taylor, limit_steps_str,
};
use io::{do_cval, do_interval, do_latex, do_lean, do_mathml, do_serialize};
use numtheory::{do_base, do_big, do_dec, do_fast, do_numtheory};
use solve::{
    do_brent, do_bspline, do_chebyshev, do_cholesky, do_cond, do_conv, do_cross, do_det, do_dot,
    do_dsolve, do_eig, do_eigen, do_fft, do_fourier, do_hermite, do_hessenberg, do_inv,
    do_isolate_roots, do_legendre, do_linsolve, do_lu, do_mc, do_minimize, do_msolve, do_nsolve,
    do_null, do_pchip, do_plot, do_poly_roots, do_qr, do_rank, do_schur, do_spcg, do_spbicg,
    do_spline, do_solve, do_svd, do_symlig, do_tikhonov, do_trace, do_transpose_cmd,
};
use stats::{
    do_anova, do_boot, do_chitest, do_dist, do_fit, do_kw, do_logit, do_mwu, do_qtile, do_sample,
    do_stats, do_ttest, do_wilcoxon,
};

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

    // d/dx: unsolve-style alias (fixed variable x)
    if let Some(rest) = line.strip_prefix("d/dx ") {
        return diff_steps(rest.trim(), ctx, Some("x"));
    }

    // diff: show original, derivative (unsimplified), simplified
    if let Some(rest) = line.strip_prefix("diff ") {
        return diff_steps(rest.trim(), ctx, None);
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

    // expand-in: expand in one variable, others become exact coefficients
    if let Some(rest) = line.strip_prefix("expand-in ") {
        let ForVar { body, variable } = strip_for_var(rest.trim());
        let Some(var) = variable else {
            return Err(crate::error::MathError::Eval(
                "expand-in needs: <expr> for <var>".into(),
            ));
        };
        let e = Parser::parse(&body)?;
        let s = crate::poly::expand_in(&e, &var);
        return Ok(vec![
            format!("expand in {var}: {}", e),
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

    // nsolve: find all real roots by scanning for sign changes
    if let Some(rest) = line.strip_prefix("nsolve ") {
        return nsolve_steps(rest.trim(), ctx);
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
        "bspline ", "hermite ", "brent ",
        "sample ", "dist ", "qtile ", "fit ", "ttest ", "chitest ", "anova ", "mwu ", "wilcoxon ",
        "kw ", "boot ", "logit ", "pdiff ", "gradient ", "let ", "fn ",
        "gcd ", "lcm ", "is-prime ", "factor ", "fib ", "binom ",
        "fact ", "mr-prime ", "jacobi ", "cf ", "diophantine ", "dlog ", "base ",
        "det ",
        "inv ",
        "transpose ",
        "trace ",
        "eigen ",
        "linsolve ",
        "msolve ",
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
        "lean ",
        "sum ",
        "prod ",
        "spcg ",
        "spbicg ",
        "idiff ",
        "dot ",
        "cross ",
        "d/dx ",
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

struct DiffParse {
    expr_src: String,
    wrt: String,
    order: u32,
}

/// A command body split into its expression part and an optional
/// `for <var>` suffix.
struct ForVar {
    /// The body with the `for <var>` suffix removed.
    body: String,
    /// The variable named by the suffix, if present.
    variable: Option<String>,
}

/// Strip a trailing `for <name>` suffix from a command body.
///
/// The suffix names the free variable of the operation: the differentiation
/// variable for `diff`, the integration variable for `integrate`, the
/// expansion variable for `taylor`, the limit variable for `limit`, and the
/// unknown to isolate for `solve`.
///
/// Without this, `for` and the name are parsed as ordinary identifiers and
/// silently multiply into the expression (`diff y*x^2 for y` became
/// `x^2*for`). Cross-port pattern from unsolve `solve-problem.ts:119-124`.
///
/// The keyword is matched case-insensitively but the name keeps its case.
/// The suffix must be the final token group; `for` preceded by whitespace is
/// required so identifiers like `xfor` are untouched.
fn strip_for_var(body: &str) -> ForVar {
    let trimmed = body.trim();
    // Split off the trailing name.
    let mut tokens = trimmed.rsplitn(2, char::is_whitespace);
    let Some(last) = tokens.next() else {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    };
    let Some(head) = tokens.next() else {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    };
    if last.is_empty() || !last.chars().all(|c| c.is_ascii_alphabetic()) {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    }
    // The word before the name must be the keyword `for`.
    let head_trimmed = head.trim_end();
    let Some(keyword_end) = head_trimmed.rfind(char::is_whitespace) else {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    };
    let keyword = &head_trimmed[keyword_end + 1..];
    if !keyword.eq_ignore_ascii_case("for") {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    }
    let expr_src = head_trimmed[..keyword_end].trim();
    if expr_src.is_empty() {
        return ForVar {
            body: trimmed.to_string(),
            variable: None,
        };
    }
    ForVar {
        body: expr_src.to_string(),
        variable: Some(last.to_string()),
    }
}

/// A pair of integration/minimization bounds.
struct Bounds {
    /// The expression to operate on.
    expr_src: String,
    /// The lower bound, as source text (an expression, e.g. `pi/2`).
    a_src: String,
    /// The upper bound, as source text.
    b_src: String,
}

/// Split a command body into an expression and a pair of bounds.
///
/// Two forms are accepted:
///
/// 1. **Positional** — `<expr> <a> <b>`, the historical form. The bounds
///    are the last two tokens.
/// 2. **Suffix** — `<expr> from <a> to <b>`, where `<a>` and `<b>` are
///    full *expressions* and may contain spaces (`from pi/2 to pi`).
///
/// The suffix form is unambiguous, so it is preferred when present; it also
/// retires the trailing-token ambiguity class, where a stray extra token
/// silently changed which numbers were read as the bounds.
///
/// When no explicit `from … to …` is given, the positional form is used and
/// the *last two tokens* are the bounds — matching the long-standing
/// behaviour that `romberg x^2 0 1 12` integrates `x^2*0` over `[1, 12]`.
///
/// Cross-port pattern from unsolve `solve-problem.ts:115-133`.
fn split_bounds(body: &str) -> Result<Bounds> {
    // The `from … to …` suffix: find the *last* `from` that is followed by a
    // `to`, so an expression containing the word "from" does not confuse it.
    let lower = body.to_ascii_lowercase();
    if let Some(from_pos) = lower.rfind(" from ") {
        let after_from = &body[from_pos + " from ".len()..];
        let after_lower = &lower[from_pos + " from ".len()..];
        if let Some(to_pos) = after_lower.rfind(" to ") {
            let a_src = after_from[..to_pos].trim().to_string();
            let b_src = after_from[to_pos + " to ".len()..].trim().to_string();
            let expr_src = body[..from_pos].trim().to_string();
            if !a_src.is_empty() && !b_src.is_empty() && !expr_src.is_empty() {
                return Ok(Bounds {
                    expr_src,
                    a_src,
                    b_src,
                });
            }
        }
        // A `from` was present but the suffix is malformed — say so rather
        // than silently treating "from" as part of the expression.
        return Err(crate::error::MathError::Eval(
            "bounds must be given as `from <a> to <b>`".into(),
        ));
    }
    // Positional: the last two tokens are the bounds.
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "needs: <expr> <a> <b>  or  <expr> from <a> to <b>".into(),
        ));
    }
    Ok(Bounds {
        expr_src: tokens[..tokens.len() - 2].join(" "),
        a_src: tokens[tokens.len() - 2].to_string(),
        b_src: tokens[tokens.len() - 1].to_string(),
    })
}

/// Evaluate a bound given as source text.
///
/// The text is parsed as an expression, so `pi`, `pi/2`, `1+1` and `2*pi`
/// all work. Constant names (`pi`, `e`, `tau`) resolve from the standard
/// context; other free variables resolve from `ctx`.
fn eval_bound(src: &str, ctx: &Context) -> Result<f64> {
    let e = Parser::parse(src)?;
    crate::eval::eval(&e, ctx)
}

/// Parse `diff <expr> [var] [order]` or `d/dx <expr> [order]` when `force_wrt` is set.
///
/// A trailing `for <var>` suffix names the differentiation variable and wins
/// over the trailing single-letter heuristic.
fn parse_diff_tail(rest: &str, force_wrt: Option<&str>) -> Result<DiffParse> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`diff` needs an expression".into()));
    }
    let mut order = 1u32;
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        if let Ok(n) = tokens[expr_end - 1].parse::<u32>() {
            if n == 0 {
                return Err(crate::error::MathError::Eval(
                    "derivative order must be at least 1".into(),
                ));
            }
            if n > 64 {
                return Err(crate::error::MathError::Eval(
                    "derivative order capped at 64".into(),
                ));
            }
            order = n;
            expr_end -= 1;
        }
    }
    let mut wrt = force_wrt.map(String::from).unwrap_or_else(|| "x".to_string());
    if force_wrt.is_none() {
        if let Some(v) = variable {
            // An explicit `for <var>` names the differentiation variable.
            wrt = v;
        } else if expr_end > 1 {
            let candidate = tokens[expr_end - 1];
            if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
                wrt = candidate.to_string();
                expr_end -= 1;
            }
        }
    }
    if expr_end == 0 {
        return Err(crate::error::MathError::Eval(
            "`diff` needs an expression".into(),
        ));
    }
    Ok(DiffParse {
        expr_src: tokens[..expr_end].join(" "),
        wrt,
        order,
    })
}

fn diff_expr_prepared(parse: &DiffParse, ctx: &Context) -> Result<Expr> {
    let e = Parser::parse(&parse.expr_src)?;
    let bound: Vec<(String, Expr)> = ctx
        .vars
        .iter()
        .filter(|(k, _)| **k != parse.wrt && !["pi", "e", "tau", "inf"].contains(&k.as_str()))
        .map(|(k, v)| (k.clone(), Expr::num(*v)))
        .collect();
    let mut e = e;
    for (k, v) in &bound {
        e = e.substitute(k, v);
    }
    Ok(e)
}

fn diff_steps(rest: &str, ctx: Context, force_wrt: Option<&str>) -> Result<Vec<String>> {
    let parse = parse_diff_tail(rest, force_wrt)?;
    let e = diff_expr_prepared(&parse, &ctx)?;
    let d = if parse.order == 1 {
        differentiate(&e, &parse.wrt)?
    } else {
        crate::symbolic::differentiate_n(&e, &parse.wrt, parse.order)?
    };
    let s = simplify(&d);
    let mut steps = vec![
        format!("f({}) = {}", parse.wrt, e),
        if parse.order == 1 {
            format!("d/d{} f({})", parse.wrt, parse.wrt)
        } else {
            format!(
                "d^{}/d{}^{} f({})",
                parse.order, parse.wrt, parse.order, parse.wrt
            )
        },
        format!("= {}", d),
        format!("simplified = {}", s),
    ];
    Ok(steps)
}

fn solve_steps(rest: &str, ctx: Context) -> Result<Vec<String>> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`solve` needs an expression".into()));
    }
    // An explicit `for <var>` names the unknown to isolate.
    let explicit_var = variable.is_some();
    for consume in [2usize, 1, 0] {
        if tokens.len() <= consume {
            continue;
        }
        let mut guess = 1.0;
        let mut wrt = variable.clone().unwrap_or_else(|| "x".to_string());
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
            if explicit_var {
                // The variable is already fixed by `for <var>`; the trailing
                // token slot is consumed by it, so a single-letter candidate
                // here would be ambiguous — skip that slot.
                if expr_end > 0
                    && tokens[expr_end - 1].len() == 1
                    && tokens[expr_end - 1].chars().all(|c| c.is_ascii_alphabetic())
                {
                    expr_end -= 1;
                }
            } else {
                let candidate = tokens[expr_end - 1];
                if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
                    wrt = candidate.to_string();
                    expr_end -= 1;
                } else {
                    continue;
                }
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

/// `nsolve <expr> [var] [lo hi]` — step-by-step output for the sign-scan
/// root finder.
fn nsolve_steps(rest: &str, ctx: Context) -> Result<Vec<String>> {
    let mut ctx = ctx;
    let out = do_nsolve(rest, &mut ctx)?
        .ok_or_else(|| crate::error::MathError::Eval("`nsolve` produced no output".into()))?;
    let mut steps = vec![format!("nsolve: {}", out)];
    if !out.starts_with("no real roots") {
        // Show each root on its own line for a readable step board.
        let list = out.strip_prefix("roots: ").unwrap_or(&out);
        for r in list.split(", ") {
            steps.push(format!("root ≈ {}", r.trim()));
        }
    }
    Ok(steps)
}

fn taylor_steps(rest: &str) -> Result<Vec<String>> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
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
    // An explicit `for <var>` names the expansion variable.
    let wrt = variable.unwrap_or_else(|| "x".to_string());
    let series = crate::taylor::taylor_series_str(&expr_src, &wrt, around, order)?;
    Ok(vec![
        format!("f({}) = {}", wrt, expr_src),
        format!("Taylor expansion around a = {}", format_value(around)),
        format!("order = {}", order),
        format!("f({}) ≈ {}", wrt, series),
    ])
}

fn integrate_sym_steps(rest: &str) -> Result<Vec<String>> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("integrate needs: <expr> [var]".into()));
    }
    // An explicit `for <var>` names the integration variable and wins over
    // the trailing single-letter heuristic.
    let explicit_var = variable.is_some();
    let var = match &variable {
        Some(v) => v.clone(),
        None => {
            if tokens.len() >= 2
                && tokens[tokens.len() - 1].len() == 1
                && tokens[tokens.len() - 1].chars().all(|c| c.is_ascii_alphabetic())
            {
                tokens[tokens.len() - 1].to_string()
            } else {
                "x".to_string()
            }
        }
    };
    let expr_end = if !explicit_var
        && tokens.len() >= 2
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
    if line == "base" {
        return Err(crate::error::MathError::Eval(
            "base needs: <n> <target>   or   <n> from <base> to <base>".into(),
        ));
    }

    if let Some(rest) = line.strip_prefix("let ") {
        return bind_var(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("fn ") {
        return define_fn(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("d/dx ") {
        return do_diff(rest.trim(), Some("x"), ctx);
    }
    if let Some(rest) = line.strip_prefix("diff ") {
        return do_diff(rest.trim(), None, ctx);
    }
    if let Some(rest) = line.strip_prefix("dot ") {
        return do_dot(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("cross ") {
        return do_cross(rest.trim());
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
    if let Some(rest) = line.strip_prefix("nsolve ") {
        return do_nsolve(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("brent ") {
        return do_brent(rest.trim(), ctx);
    }
    if let Some(rest) = line.strip_prefix("dsolve ") {
        return do_dsolve(rest.trim());
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
    if let Some(rest) = line.strip_prefix("base ") {
        return do_base(rest.trim());
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
    if let Some(rest) = line.strip_prefix("inv ") {
        return do_inv(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("transpose ") {
        return do_transpose_cmd(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("trace ") {
        return do_trace(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("eigen ") {
        return do_eigen(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("linsolve ") {
        return do_linsolve(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("msolve ") {
        return do_msolve(rest.trim());
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
    if let Some(rest) = line.strip_prefix("expand-in ") {
        // `expand-in <expr> for <var>` — expand in one variable only,
        // treating everything else as an opaque exact coefficient.
        let ForVar { body, variable } = strip_for_var(rest.trim());
        let Some(var) = variable else {
            return Err(crate::error::MathError::Eval(
                "expand-in needs: <expr> for <var>".into(),
            ));
        };
        let e = Parser::parse(&body)?;
        return Ok(Some(crate::poly::expand_in(&e, &var).to_string()));
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
    if let Some(rest) = line.strip_prefix("latex ") {
        return do_latex(rest.trim());
    }
    if let Some(rest) = line.strip_prefix("lean ") {
        return do_lean(rest.trim());
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

    // A bare equation (`x^2 - 4 = 0`) routes to `solve` — the user clearly
    // wants the roots, not a parse error. The `=` must be at the top level
    // (not inside parentheses) so `f(x) = ...` style input is not hijacked.
    if let Some(eq_pos) = top_level_equals(line) {
        let lhs = line[..eq_pos].trim();
        let rhs = line[eq_pos + 1..].trim();
        if !lhs.is_empty() && !rhs.is_empty() {
            // Prefer the symbolic solver: it handles polynomials exactly
            // (complex roots, radicals). Fall back to the numeric solver
            // for equations that are not polynomials.
            let equation = format!("{} = {}", lhs, rhs);
            match crate::qsolve::qsolve_steps_str(&equation) {
                // The last step carries the roots; the earlier ones are the
                // working (equation, coefficients, discriminant).
                Ok(steps) => {
                    return Ok(Some(steps.last().cloned().unwrap_or_default()));
                }
                Err(_) => return do_solve(&equation, ctx),
            }
        }
    }

    // Default: evaluate the expression and print the value
    let e = Parser::parse(line)?;
    let v = eval(&e, ctx)?;
    Ok(Some(format_value(v)))
}

/// Find the position of a top-level `=` in `s`, if there is one.
///
/// "Top level" means not nested in parentheses or brackets — so
/// `f(x) = 1` is found, but `(a = b)` is not.
fn top_level_equals(s: &str) -> Option<usize> {
    let mut depth = 0i32;
    for (i, c) in s.char_indices() {
        match c {
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth -= 1,
            '=' if depth == 0 => {
                // `==` is not an equation.
                if s[i + 1..].starts_with('=') {
                    return None;
                }
                // A `!=` is a comparison, not an equation.
                if i > 0 && s[..i].ends_with('!') {
                    return None;
                }
                return Some(i);
            }
            _ => {}
        }
    }
    None
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

fn parse_f64_list(rest: &str) -> Result<Vec<f64>> {
    rest.split(|c: char| c == ',' || c.is_whitespace())
        .filter(|t| !t.trim().is_empty())
        .map(|t| t.trim().parse::<f64>().map_err(|_| crate::error::MathError::Eval(format!("not a number: {}", t))))
        .collect()
}

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

/// Format a floating-point value for display.
///
/// Integral values print without a fractional part; others use 10
/// significant decimals with trailing zeros trimmed. `NaN` and the
/// infinities get their conventional spellings.
///
/// Public so external surfaces (e.g. the MCP server) format numbers the
/// same way the REPL does.
pub fn format_value(v: f64) -> String {
    if v.is_nan() {
        "NaN".into()
    } else if v.is_infinite() {
        if v > 0.0 { "inf".into() } else { "-inf".into() }
    } else if v == v.trunc() && v.abs() < 1e15 {
        // `-0.0` would otherwise print as "-0" — a signed zero carries no
        // meaning in a result, so normalise it.
        let i = v as i64;
        if i == 0 { "0".into() } else { format!("{}", i) }
    } else if v != 0.0 && v.abs() < 1e-10 {
        // A value below the 10-decimal format's resolution would round away
        // to "0", so use scientific notation to keep every digit.
        format!("{:e}", v)
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
  diff <expr> [var] [order]  symbolic derivative (order ≤ 64)
  d/dx <expr> [order]        derivative with respect to x (unsolve alias)
  dot a1 a2 ... | b1 b2 ...  dot product
  cross a1 a2 a3 | b1 b2 b3  3-D cross product
  pdiff <expr> <var>  partial derivative
  gradient <expr>     gradient (all partials)
  idiff <expr> [= rhs] [<indep> <dep>]
                      implicit differentiation dy/dx = -F_x/F_y
  simplify <expr>     constant-fold & simplify
  int <expr> a b      numerical integral over [a, b]
  solve <expr> [var] [guess]
  brent <expr> [var] <lo> <hi>
                      Brent root finding in the bracket [lo, hi]
                      (keeps the sign change, converges faster than bisection)
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
  inv <rows>          matrix inverse
  transpose <rows>    matrix transpose
  trace <rows>        sum of diagonal entries
  eigen <rows>        symmetric eigendecomposition (λ ascending, vectors as columns)
  linsolve [A|b]      solve A·x = b from the augmented matrix (last column = b)
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
  lean <expr>        export expression to a Lean 4 file (needs Mathlib)
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

    // ---- `from a to b` bounds with expression values ----------------------

    #[test]
    fn split_bounds_positional_form() {
        let b = split_bounds("sin(x) 0 pi").unwrap();
        assert_eq!(b.expr_src, "sin(x)");
        assert_eq!(b.a_src, "0");
        assert_eq!(b.b_src, "pi");
    }

    #[test]
    fn split_bounds_suffix_form() {
        let b = split_bounds("sin(x) from 0 to pi").unwrap();
        assert_eq!(b.expr_src, "sin(x)");
        assert_eq!(b.a_src, "0");
        assert_eq!(b.b_src, "pi");
    }

    #[test]
    fn split_bounds_suffix_allows_spaces_in_bounds() {
        // The bounds are full expressions, so they may contain spaces.
        let b = split_bounds("sin(x) from pi / 2 to 2 * pi").unwrap();
        assert_eq!(b.expr_src, "sin(x)");
        assert_eq!(b.a_src, "pi / 2");
        assert_eq!(b.b_src, "2 * pi");
    }

    #[test]
    fn split_bounds_suffix_needs_both_endpoints() {
        // A malformed suffix is an error, not a silent fall-through.
        assert!(split_bounds("sin(x) from 0").is_err());
        assert!(split_bounds("sin(x) from to pi").is_err());
    }

    #[test]
    fn bounds_accept_expression_values() {
        let ctx = Context::standard();
        assert_eq!(eval_bound("pi", &ctx).unwrap(), std::f64::consts::PI);
        assert_eq!(eval_bound("pi/2", &ctx).unwrap(), std::f64::consts::PI / 2.0);
        assert_eq!(eval_bound("2*pi", &ctx).unwrap(), 2.0 * std::f64::consts::PI);
        assert_eq!(eval_bound("1+1", &ctx).unwrap(), 2.0);
    }

    #[test]
    fn int_accepts_expression_bounds() {
        // `int sin(x) 0 pi/2` = 1
        let out = dispatch_str("int sin(x) 0 pi/2", Context::standard()).unwrap().unwrap();
        assert!(out.contains("= 1"), "got: {out}");
    }

    #[test]
    fn romberg_accepts_expression_bounds() {
        // romberg used to reject `pi` outright ("bad b").
        let out = dispatch_str("romberg sin(x) 0 pi", Context::standard()).unwrap().unwrap();
        assert!(out.contains("= 2"), "got: {out}");
    }

    #[test]
    fn from_to_suffix_works_for_all_three() {
        let ctx = Context::standard();
        let out = dispatch_str("int sin(x) from 0 to pi", ctx.clone()).unwrap().unwrap();
        assert!(out.contains("= 2"), "got: {out}");
        let out = dispatch_str("romberg sin(x) from 0 to pi", ctx.clone()).unwrap().unwrap();
        assert!(out.contains("= 2"), "got: {out}");
        let out = dispatch_str("minimize (x-3)^2 from 0 to 10", ctx).unwrap().unwrap();
        assert!(out.contains("x* = 3"), "got: {out}");
    }

    #[test]
    fn from_to_suffix_is_equivalent_to_positional() {
        let ctx = Context::standard();
        let a = dispatch_str("int x^2 from 0 to 3", ctx.clone()).unwrap().unwrap();
        let b = dispatch_str("int x^2 0 3", ctx).unwrap().unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn from_to_suffix_handles_fraction_bounds() {
        let out = dispatch_str("int sin(x) from 0 to pi/2", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("= 1"), "got: {out}");
    }

    #[test]
    fn positional_bounds_still_take_the_last_two_tokens() {
        // The historical behaviour is unchanged: the last two tokens are the
        // bounds, and extra tokens fold into the expression.
        let out = dispatch_str("romberg x^2 0 1 12", Context::standard()).unwrap().unwrap();
        assert!(out.contains("(1, 12;"), "got: {out}");
    }

    // ---- BigInt auto-upgrade for is-prime / factor --------------------------

    #[test]
    fn is_prime_auto_upgrades_to_bigint() {
        // 10^20 is composite — it used to fail with "could not parse as integer".
        let out = dispatch_str("is-prime 100000000000000000000", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "false");
    }

    #[test]
    fn factor_auto_upgrades_to_bigint() {
        // 10^20 = 2^20 * 5^20
        let out = dispatch_str("factor 100000000000000000000", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "2^20 · 5^20");
    }

    #[test]
    fn is_prime_handles_a_big_mersenne_prime() {
        // 2^127 - 1 is prime.
        let out = dispatch_str(
            "is-prime 170141183460469231731687303715884105727",
            Context::standard(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(out, "true");
    }

    #[test]
    fn small_values_still_use_the_u64_path() {
        let ctx = Context::standard();
        let out = dispatch_str("is-prime 97", ctx.clone()).unwrap().unwrap();
        assert_eq!(out, "true");
        let out = dispatch_str("factor 360", ctx).unwrap().unwrap();
        assert_eq!(out, "2^3 · 3^2 · 5");
    }

    // ---- number-base conversion ---------------------------------------------

    #[test]
    fn base_converts_to_named_bases() {
        let ctx = Context::standard();
        let out = dispatch_str("base 255 hex", ctx.clone()).unwrap().unwrap();
        assert_eq!(out, "FF_16");
        let out = dispatch_str("base 255 bin", ctx.clone()).unwrap().unwrap();
        assert_eq!(out, "11111111_2");
        let out = dispatch_str("base 8 oct", ctx).unwrap().unwrap();
        assert_eq!(out, "10_8");
    }

    #[test]
    fn base_converts_to_numeric_base() {
        let ctx = Context::standard();
        let out = dispatch_str("base 255 16", ctx.clone()).unwrap().unwrap();
        assert_eq!(out, "FF_16");
        let out = dispatch_str("base 42 36", ctx).unwrap().unwrap();
        assert_eq!(out, "16_36");
    }

    #[test]
    fn base_converts_between_bases() {
        let ctx = Context::standard();
        let out = dispatch_str("base 255 from 10 to 2", ctx.clone()).unwrap().unwrap();
        assert_eq!(out, "11111111_2");
        let out = dispatch_str("base FF from 16 to 10", ctx).unwrap().unwrap();
        assert_eq!(out, "255_10");
    }

    #[test]
    fn base_handles_negative_values() {
        let out = dispatch_str("base -255 hex", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "-FF_16");
    }

    #[test]
    fn base_rejects_bad_input() {
        let ctx = Context::standard();
        assert!(dispatch_str("base 255 base99", ctx.clone()).is_err());
        assert!(dispatch_str("base 255 nonsense", ctx.clone()).is_err());
        assert!(dispatch_str("base", ctx.clone()).is_err());
        assert!(dispatch_str("base 255", ctx).is_err());
    }

    // ---- bare equations and LaTeX output ------------------------------------

    #[test]
    fn bare_equation_routes_to_solve() {
        // `x^2 - 4 = 0` used to be a parse error.
        let out = dispatch_str("x^2 - 4 = 0", Context::standard()).unwrap().unwrap();
        assert!(out.contains('2'), "got: {out}");
        assert!(out.contains("-2"), "got: {out}");
    }

    #[test]
    fn bare_equation_solves_exactly() {
        // A polynomial equation goes through the symbolic solver.
        let out = dispatch_str("x^2+1=0", Context::standard()).unwrap().unwrap();
        assert!(out.contains('i'), "got: {out}");
        // Radicals stay exact.
        let out = dispatch_str("x^2-2=0", Context::standard()).unwrap().unwrap();
        assert!(out.contains("sqrt(2)"), "got: {out}");
    }

    #[test]
    fn bare_equation_falls_back_to_numeric_solve() {
        // A non-polynomial equation uses the numeric solver.
        let out = dispatch_str("exp(x) = 5", Context::standard()).unwrap().unwrap();
        assert!(out.contains("1.609"), "got: {out}");
    }

    #[test]
    fn leftover_equals_suggests_solve() {
        // The error names the offending character and points at `solve`.
        let err = dispatch_str("diff x^2 = 4", Context::standard()).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains('='), "got: {msg}");
        assert!(msg.contains("solve"), "got: {msg}");
    }

    #[test]
    fn latex_output_round_trips() {
        // The point of the `latex` command: the TeX it emits parses back to
        // an equal expression.
        for src in ["x^2 + 1", "x^2/2", "sin(x)", "x^2 - 2*y", "1/2 + 3/4"] {
            let out = dispatch_str(&format!("latex {src}"), Context::standard())
                .unwrap()
                .unwrap();
            let back = Parser::parse(&out).unwrap();
            let orig = Parser::parse(src).unwrap();
            assert!(orig.equals(&back), "{src} -> {out} -> {}", back);
        }
    }

    #[test]
    fn top_level_equals_finds_the_equation_sign() {
        assert!(top_level_equals("x = 1").is_some());
        assert!(top_level_equals("x^2 - 4 = 0").is_some());
        // Not nested in parentheses.
        assert!(top_level_equals("(a = b)").is_none());
        // A comparison, not an equation.
        assert!(top_level_equals("x != 1").is_none());
        assert!(top_level_equals("x == 1").is_none());
        // None at all.
        assert!(top_level_equals("x + 1").is_none());
    }

    #[test]
    fn latex_command_renders_tex() {
        let out = dispatch_str("latex x^2/2", Context::standard()).unwrap().unwrap();
        assert_eq!(out, r"\frac{x^{2}}{2}");
        let out = dispatch_str("latex sin(x)", Context::standard()).unwrap().unwrap();
        assert!(out.contains(r"\sin"), "got: {out}");
    }

    // ---- the `at` keyword for limits ----------------------------------------

    #[test]
    fn limit_at_keyword_infers_the_variable() {
        // Regression: `at` was read as the *variable*, so `limit sin(x)/x at 0`
        // computed the limit of sin(x)/x with respect to a variable named
        // "at" and reported "does not exist".
        let out = dispatch_str("limit sin(x)/x at 0", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("= 1"), "got: {out}");
        assert!(out.contains("as x"), "got: {out}");
    }

    #[test]
    fn limit_at_keyword_matches_the_positional_form() {
        let ctx = Context::standard();
        let a = dispatch_str("limit sin(x)/x at 0", ctx.clone()).unwrap().unwrap();
        let b = dispatch_str("limit sin(x)/x x 0", ctx).unwrap().unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn limit_at_keyword_handles_a_pole() {
        // 1/x at 0: the two sides diverge oppositely.
        let out = dispatch_str("limit 1/x at 0", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("does not exist"), "got: {out}");
    }

    #[test]
    fn limit_at_keyword_with_explicit_variable() {
        // `for <var>` still wins over the inference.
        let out = dispatch_str("limit sin(t)/t at 0 for t", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("= 1"), "got: {out}");
        assert!(out.contains("as t"), "got: {out}");
    }

    // ---- numeric domain probes --------------------------------------------

    #[test]
    fn minimize_rejects_undefined_endpoint() {
        // sqrt(x) is not defined at x = -1 — the search must not silently
        // report a meaningless minimum there.
        let err = dispatch_str("minimize sqrt(x) -1 5", Context::standard()).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("not defined"), "got: {msg}");
        assert!(msg.contains("-1"), "got: {msg}");
    }

    #[test]
    fn minimize_rejects_undefined_log_endpoint() {
        let err = dispatch_str("minimize log(x) -5 5", Context::standard()).unwrap_err();
        assert!(err.to_string().contains("not defined"), "got: {err}");
    }

    #[test]
    fn minimize_accepts_a_defined_interval() {
        // 1/x is defined on [-5, -1] (the pole at 0 is outside).
        let out = dispatch_str("minimize 1/x -5 -1", Context::standard()).unwrap().unwrap();
        assert!(out.contains("x* = -1"), "got: {out}");
    }

    #[test]
    fn minimize_still_finds_the_minimum() {
        let out = dispatch_str("minimize (x-2)^2 0 10", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("x* = 2"), "got: {out}");
    }

    #[test]
    fn bounds_error_message_lists_both_forms() {
        let err = dispatch_str("int sin(x)", Context::standard()).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("from"), "got: {msg}");
        assert!(msg.contains("<a> <b>"), "got: {msg}");
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

    // ---- `for <var>` suffix ------------------------------------------------

    #[test]
    fn strip_for_var_extracts_name() {
        let fv = strip_for_var("y*x^2 for y");
        assert_eq!(fv.body, "y*x^2");
        assert_eq!(fv.variable.as_deref(), Some("y"));
    }

    #[test]
    fn strip_for_var_ignores_non_suffix() {
        // `for` must be the final token group.
        let fv = strip_for_var("x^2 for x extra");
        assert_eq!(fv.body, "x^2 for x extra");
        assert!(fv.variable.is_none());
    }

    #[test]
    fn strip_for_var_ignores_bare_for() {
        // `for` with no name is not a suffix.
        let fv = strip_for_var("x^2 for");
        assert_eq!(fv.body, "x^2 for");
        assert!(fv.variable.is_none());
    }

    #[test]
    fn strip_for_var_ignores_non_identifier_name() {
        // The name must be letters only.
        let fv = strip_for_var("x^2 for 3");
        assert_eq!(fv.body, "x^2 for 3");
        assert!(fv.variable.is_none());
    }

    #[test]
    fn strip_for_var_is_case_insensitive_keyword() {
        let fv = strip_for_var("y^2 FOR y");
        assert_eq!(fv.body, "y^2");
        assert_eq!(fv.variable.as_deref(), Some("y"));
    }

    #[test]
    fn strip_for_var_leaves_expressions_alone() {
        // A variable that merely contains "for" is not a suffix.
        let fv = strip_for_var("xfor*y");
        assert_eq!(fv.body, "xfor*y");
        assert!(fv.variable.is_none());
    }

    #[test]
    fn for_var_diff_names_variable() {
        // `for` must not be parsed as letters multiplying into the expression.
        let out = dispatch_str("diff y*x^2 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "x^2");
    }

    #[test]
    fn for_var_diff_does_not_leak_for_as_variable() {
        // Regression: `diff y*x^2 for y` used to yield "x^2*for".
        let out = dispatch_str("diff y*x^2 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert!(!out.contains("for"), "got {}", out);
    }

    #[test]
    fn for_var_integrate_names_variable() {
        let out = dispatch_str("integrate y^2 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("dy"), "got {}", out);
        assert!(out.contains("y^3"), "got {}", out);
    }

    #[test]
    fn for_var_taylor_names_variable() {
        let out = dispatch_str("taylor exp(y) for y", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains('y'), "got {}", out);
        assert!(!out.contains("unknown variable"), "got {}", out);
    }

    #[test]
    fn for_var_limit_names_variable() {
        let out = dispatch_str("limit sin(y)/y at 0 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("→ 0"), "got {}", out);
        assert!(out.ends_with("= 1"), "got {}", out);
    }

    #[test]
    fn for_var_solve_names_unknown() {
        // `for y` must target `y`, so the numeric solver differentiates
        // with respect to `y` (not `x`): y^2 - 4 = 0 has roots ±2.
        let out = dispatch_str("solve y^2 - 4 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("root ≈ 2"), "got {}", out);
        // The same equation solved for `x` would not converge (x is unbound).
        let steps = dispatch_steps("solve y^2 - 4 for y", Context::standard()).unwrap();
        assert!(steps[0].contains("y"), "got {:?}", steps);
    }

    #[test]
    fn for_var_diff_steps() {
        let steps = dispatch_steps("diff y*x^2 for y", Context::standard()).unwrap();
        assert!(steps[0].contains("y*x^2"), "got {:?}", steps);
        assert!(steps.join(" ").contains("x^2"), "got {:?}", steps);
        assert!(!steps.join(" ").contains("*for"), "got {:?}", steps);
    }

    #[test]
    fn for_var_integrate_steps() {
        let steps = dispatch_steps("integrate y^2 for y", Context::standard()).unwrap();
        assert!(steps[0].contains("y^2"), "got {:?}", steps);
        assert!(steps[1].contains("dy"), "got {:?}", steps);
    }

    #[test]
    fn for_var_absent_keeps_existing_behaviour() {
        // No `for` suffix — the trailing single-letter heuristic still applies.
        let out = dispatch_str("diff x^2", Context::standard()).unwrap().unwrap();
        assert_eq!(out, "2*x");
        let out = dispatch_str("integrate x^2", Context::standard())
            .unwrap()
            .unwrap();
        assert!(out.contains("dx"), "got {}", out);
    }

    #[test]
    fn for_var_diff_order_still_works() {
        // `for` and the order argument compose: d²/dy² (y*x²) = 0.
        let out = dispatch_str("diff y*x^2 for y 2", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "0");
        // The first derivative is still x².
        let out = dispatch_str("diff y*x^2 for y", Context::standard())
            .unwrap()
            .unwrap();
        assert_eq!(out, "x^2");
    }
}