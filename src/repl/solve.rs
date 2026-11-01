//! Solving, numerics, linear-algebra and signal command handlers.
//!
//! Handlers moved verbatim from the original monolithic `repl.rs`.

use super::*;
pub(super) fn do_dot(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() != 2 {
        return Err(crate::error::MathError::Eval(
            "`dot` needs two vectors: `dot a1 a2 ... | b1 b2 ...`".into(),
        ));
    }
    let a = parse_group(groups[0])?;
    let b = parse_group(groups[1])?;
    let v = crate::vector::dot(&a, &b)?;
    Ok(Some(format_value(v)))
}

pub(super) fn do_cross(rest: &str) -> Result<Option<String>> {
    let groups = split_groups(rest)?;
    if groups.len() != 2 {
        return Err(crate::error::MathError::Eval(
            "`cross` needs two 3-vectors: `cross a1 a2 a3 | b1 b2 b3`".into(),
        ));
    }
    let a = parse_group(groups[0])?;
    let b = parse_group(groups[1])?;
    let c = crate::vector::cross(&a, &b)?;
    Ok(Some(format!(
        "[{}, {}, {}]",
        format_value(c[0]),
        format_value(c[1]),
        format_value(c[2])
    )))
}

pub(super) fn do_solve(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let ForVar { body, variable } = strip_for_var(rest);
    // An equation `lhs = rhs` is solved as `lhs − rhs = 0`.
    let body = match top_level_equals(&body) {
        Some(eq) => {
            let lhs = body[..eq].trim();
            let rhs = body[eq + 1..].trim();
            if lhs.is_empty() || rhs.is_empty() {
                return Err(crate::error::MathError::Eval(
                    "an equation needs both sides, e.g. `x^2 - 4 = 0`".into(),
                ));
            }
            let l = Parser::parse(lhs)?;
            let r = Parser::parse(rhs)?;
            crate::simplify::simplify(&Expr::sub(l, r)).to_string()
        }
        None => body,
    };
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval("`solve` needs an expression".into()));
    }
    // An explicit `for <var>` names the unknown to isolate.
    let explicit_var = variable.is_some();
    // Try parsing from the right: optional guess (number), optional var (single alpha).
    // If the expression fails to parse with trailing tokens consumed, retry with fewer consumed.
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
                // The unknown is already fixed by `for <var>`; the single-letter
                // slot would be redundant, so consume it only if present.
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

/// `brent <expr> [var] <lo> <hi>` — find a root in a bracket with Brent's
/// method.
///
/// Brent's method keeps the sign-changing bracket while interpolating, so it
/// is as robust as bisection but usually converges much faster. Unlike
/// `solve` (which needs only a guess) and `nsolve` (which scans for every
/// root), `brent` needs a bracket.
pub(super) fn do_brent(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 3 {
        return Err(crate::error::MathError::Eval(
            "`brent` needs: <expr> [var] <lo> <hi>".into(),
        ));
    }
    // The bounds are the last two tokens, so they are parsed first; that
    // leaves the optional `for <var>` suffix at the end, where
    // `strip_for_var` can see it.
    let lo = tokens[tokens.len() - 2].parse::<f64>().map_err(|_| {
        crate::error::MathError::Eval(format!("bad lower bound: {}", tokens[tokens.len() - 2]))
    })?;
    let hi = tokens[tokens.len() - 1].parse::<f64>().map_err(|_| {
        crate::error::MathError::Eval(format!("bad upper bound: {}", tokens[tokens.len() - 1]))
    })?;
    let head = tokens[..tokens.len() - 2].join(" ");
    let ForVar { body, variable } = strip_for_var(&head);
    let tokens: Vec<&str> = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`brent` needs: <expr> [var] <lo> <hi>".into(),
        ));
    }
    // A trailing single letter is the variable (unless `for <var>` named it).
    let mut wrt = variable.clone();
    let mut expr_end = tokens.len();
    if wrt.is_none() && expr_end >= 2 {
        let candidate = tokens[expr_end - 1];
        if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
            wrt = Some(candidate.to_string());
            expr_end -= 1;
        }
    }
    let expr_src = tokens[..expr_end].join(" ");
    if expr_src.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`brent` needs: <expr> [var] <lo> <hi>".into(),
        ));
    }
    let e = Parser::parse(&expr_src)?;

    // Infer the variable when it was not named: exactly one free variable.
    let wrt = match wrt {
        Some(v) => v,
        None => {
            let mut vars = e.variables();
            vars.sort();
            vars.dedup();
            if vars.len() != 1 {
                return Err(crate::error::MathError::Eval(
                    "`brent` needs a variable: the expression has multiple or no free variables"
                        .into(),
                ));
            }
            vars.remove(0)
        }
    };

    let ctx2 = ctx.clone();
    let wrt_clone = wrt.clone();
    let f = move |x: f64| {
        let mut cx = ctx2.clone();
        cx.set(&wrt_clone, x);
        crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
    };
    let (root, fval) = crate::solver::brent(f, lo, hi, crate::solver::SolveOptions::default())?;
    Ok(Some(format!(
        "root ≈ {} (f = {})",
        format_value(root),
        format_value(fval)
    )))
}

/// `nsolve <expr> [var] [lo hi]` — find **all** real roots by scanning for
/// sign changes.
///
/// Unlike `solve`, which runs one Newton–Raphson iteration from a guess and
/// therefore returns a single root, `nsolve` samples an interval, bisects
/// every sign change, and returns them all. `nsolve sin(x)` finds
/// `0, ±π, ±2π, …`.
///
/// The optional `lo hi` bounds the scan (default `-30 30`).
pub(super) fn do_nsolve(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    // Parse the bounds *first*: they come last in the command, and
    // stripping them first leaves the `for <var>` suffix at the end,
    // where `strip_for_var` can see it.
    let mut tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`nsolve` needs: <expr> [var] [lo hi]".into(),
        ));
    }
    let mut lo: Option<f64> = None;
    let mut hi: Option<f64> = None;
    if tokens.len() >= 2 {
        if let (Ok(a), Ok(b)) = (
            tokens[tokens.len() - 2].parse::<f64>(),
            tokens[tokens.len() - 1].parse::<f64>(),
        ) {
            lo = Some(a);
            hi = Some(b);
            tokens.truncate(tokens.len() - 2);
        }
    }
    // Now the `for <var>` suffix (if any) is the final token group.
    let ForVar { body, variable } = strip_for_var(&tokens.join(" "));
    tokens = body.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`nsolve` needs: <expr> [var] [lo hi]".into(),
        ));
    }
    // A trailing single letter is the variable (unless `for <var>` named it).
    let mut wrt = variable.clone();
    let mut expr_end = tokens.len();
    if wrt.is_none() && expr_end >= 2 {
        let candidate = tokens[expr_end - 1];
        if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
            wrt = Some(candidate.to_string());
            expr_end -= 1;
        }
    }
    let expr_src = tokens[..expr_end].join(" ");
    if expr_src.is_empty() {
        return Err(crate::error::MathError::Eval(
            "`nsolve` needs: <expr> [var] [lo hi]".into(),
        ));
    }
    let e = Parser::parse(&expr_src)?;

    // Infer the variable when it was not named: the expression must have
    // exactly one free variable.
    let wrt = match wrt {
        Some(v) => v,
        None => {
            let vars = e.variables();
            let mut vars = vars;
            vars.sort();
            vars.dedup();
            if vars.len() != 1 {
                return Err(crate::error::MathError::Eval(
                    "`nsolve` needs a variable: the expression has multiple or no free variables"
                        .into(),
                ));
            }
            vars.remove(0)
        }
    };

    let mut opts = crate::solver::ScanOptions::default();
    if let Some(l) = lo {
        opts.lo = l;
    }
    if let Some(h) = hi {
        opts.hi = h;
    }
    if opts.lo >= opts.hi {
        return Err(crate::error::MathError::Eval(
            "`nsolve`: lo must be less than hi".into(),
        ));
    }

    let ctx2 = ctx.clone();
    let wrt_clone = wrt.clone();
    let f = move |x: f64| {
        let mut cx = ctx2.clone();
        cx.set(&wrt_clone, x);
        crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
    };
    let roots = crate::solver::find_all_roots(f, &opts);
    if roots.is_empty() {
        return Ok(Some(format!(
            "no real roots in [{}, {}]",
            format_value(opts.lo),
            format_value(opts.hi)
        )));
    }
    let list: Vec<String> = roots.iter().map(|r| format_value(*r)).collect();
    Ok(Some(format!("roots: {}", list.join(", "))))
}

pub(super) fn do_plot(rest: &str, _ctx: &mut Context) -> Result<Option<String>> {
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

pub(super) fn do_fft(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_fourier(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
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

pub(super) fn do_mc(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
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

pub(super) fn do_cholesky(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_eig(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let result = m.power_iteration(crate::matrix::PowerIterOptions::default())?;
    let v_strs: Vec<String> = result.vector.iter().map(|x| format!("{:.6}", x)).collect();
    Ok(Some(format!(
        "λ ≈ {} (dominant eigenvalue)\nv = [{}]",
        format_value(result.value),
        v_strs.join(", ")
    )))
}

pub(super) fn do_symlig(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_hessenberg(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_schur(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_chebyshev(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_svd(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let svd = m.svd()?;
    let s_strs: Vec<String> = svd.singular_values.iter().map(|s| format!("{:.6}", s)).collect();
    Ok(Some(format!(
        "σ = [{}]\nU = {} rows × {} cols\nV = {}",
        s_strs.join(", "),
        svd.u.rows, svd.u.cols, svd.v
    )))
}

pub(super) fn do_legendre(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_conv(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_poly_roots(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_isolate_roots(rest: &str) -> Result<Option<String>> {
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
pub(super) fn parse_matrix(rest: &str) -> Result<crate::matrix::Matrix> {
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

pub(super) fn do_lu(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_qr(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_tikhonov(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_bspline(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_hermite(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_spcg(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_spbicg(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_rank(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let r = m.rank(1e-10);
    Ok(Some(format!("rank = {}", r)))
}

pub(super) fn do_cond(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let k = m.condition_number()?;
    if k.is_infinite() {
        return Ok(Some("cond = inf (singular)".into()));
    }
    Ok(Some(format!("cond = {}", format_value(k))))
}

pub(super) fn do_null(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_det(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let d = m.determinant()?;
    Ok(Some(format!("det = {}", format_value(d))))
}

/// `inv <rows>` — matrix inverse (singular matrices error).
pub(super) fn do_inv(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    let inv = m.inverse()?;
    Ok(Some(format!("A⁻¹ =\n{inv}")))
}

/// `transpose <rows>` — matrix transpose.
pub(super) fn do_transpose_cmd(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    Ok(Some(format!("Aᵀ =\n{}", m.transpose())))
}

/// `trace <rows>` — sum of diagonal entries.
pub(super) fn do_trace(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    Ok(Some(format!("tr = {}", format_value(m.trace()?))))
}

/// `eigen <rows>` — symmetric eigenvalue decomposition (ascending λ,
/// eigenvectors as columns). Symmetry is checked here (the library
/// `symmetric_eig` documents its input as unchecked).
pub(super) fn do_eigen(rest: &str) -> Result<Option<String>> {
    let m = parse_matrix(rest)?;
    if !m.is_square() {
        return Err(crate::error::MathError::Eval(
            "eigen needs a square matrix".into(),
        ));
    }
    for i in 0..m.rows {
        for j in (i + 1)..m.cols {
            if (m[(i, j)] - m[(j, i)]).abs() > 1e-12 {
                return Err(crate::error::MathError::Eval(format!(
                    "eigen needs a symmetric matrix (a[{i}][{j}] = {} ≠ a[{j}][{i}] = {})",
                    m[(i, j)],
                    m[(j, i)]
                )));
            }
        }
    }
    let (values, vectors) = m.symmetric_eig()?;
    let vals: Vec<String> = values.iter().map(|v| format_value(*v)).collect();
    Ok(Some(format!(
        "eigenvalues (ascending) = [{}]\neigenvectors (columns) =\n{}",
        vals.join(", "),
        vectors
    )))
}

/// `linsolve <rows of [A|b]>` — solve A·x = b from the augmented matrix
/// (last column is b). Convention shared with the sparse `spcg` solver.
pub(super) fn do_linsolve(rest: &str) -> Result<Option<String>> {
    let aug = parse_matrix(rest)?;
    if aug.cols < 2 || aug.rows != aug.cols - 1 {
        return Err(crate::error::MathError::Eval(
            "linsolve needs the augmented matrix [A | b]: rows × (columns − 1) = A must be square"
                .into(),
        ));
    }
    let n = aug.rows;
    let mut a_rows: Vec<Vec<f64>> = Vec::with_capacity(n);
    let mut b: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let mut row = Vec::with_capacity(n);
        for j in 0..n {
            row.push(aug[(i, j)]);
        }
        b.push(aug[(i, n)]);
        a_rows.push(row);
    }
    let a = crate::matrix::Matrix::from_rows(&a_rows)?;
    let x = a.solve(&b)?;
    let strs: Vec<String> = x.iter().map(|v| format_value(*v)).collect();
    Ok(Some(format!("x = [{}]", strs.join(", "))))
}

/// `dsolve <eq>` — solve a first-order differential equation symbolically.
///
/// The equation is `y' = <rhs>`, `dy/dx = <rhs>` or `dy/dt = <rhs>`. The
/// dependent and independent variables are taken from the equation; `y` and
/// `x` are the defaults.
pub(super) fn do_dsolve(rest: &str) -> Result<Option<String>> {
    let text = rest.trim();
    // Split on the top-level `=`.
    let Some(eq_pos) = top_level_equals(text) else {
        return Err(crate::error::MathError::Eval(
            "dsolve needs an equation, e.g. dsolve y' = x*y".into(),
        ));
    };
    let lhs = text[..eq_pos].trim();
    let rhs = text[eq_pos + 1..].trim();
    if lhs.is_empty() || rhs.is_empty() {
        return Err(crate::error::MathError::Eval(
            "dsolve needs an equation with both sides, e.g. dsolve y' = x*y".into(),
        ));
    }
    // The left-hand side names the dependent and independent variables:
    // `y'`, `dy/dx`, `dy/dt`.
    let (dep, indep) = parse_ode_lhs(lhs)?;
    let e = Parser::parse(rhs)?;
    let solution = crate::dsolve::solve_ode(&e, &indep, &dep)?;
    let mut out = solution.steps.join("\n");
    if !out.ends_with('\n') {
        out.push('\n');
    }
    Ok(Some(out))
}

/// Parse the left-hand side of a differential equation.
///
/// Accepts `y'`, `dy/dx`, `dy/dt`, `y''` — returning the dependent and
/// independent variable names.
pub(super) fn parse_ode_lhs(lhs: &str) -> Result<(String, String)> {
    let t = lhs.trim();
    // `dy/dx` — a quotient form.
    if let Some((num, den)) = t.split_once('/') {
        let num = num.trim();
        let den = den.trim();
        // The numerator may carry a derivative marker: `dy'`, `d(y)`.
        let dep = num
            .strip_prefix('d')
            .unwrap_or(num)
            .trim()
            .trim_end_matches('\'')
            .to_string();
        // The denominator may carry a `d` prefix: `dx` → `x`.
        let indep = den
            .strip_prefix('d')
            .unwrap_or(den)
            .trim()
            .to_string();
        if dep.is_empty() || indep.is_empty() {
            return Err(crate::error::MathError::Eval(format!(
                "dsolve: cannot read '{}' as a differential equation",
                lhs
            )));
        }
        return Ok((dep, indep));
    }
    // `y'` — a prime form: the independent variable defaults to x.
    if let Some(base) = t.strip_suffix('\'') {
        let base = base.trim();
        if base.is_empty() {
            return Err(crate::error::MathError::Eval(format!(
                "dsolve: cannot read '{}' as a differential equation",
                lhs
            )));
        }
        return Ok((base.to_string(), "x".to_string()));
    }
    Err(crate::error::MathError::Eval(format!(
        "dsolve: cannot read '{}' as a differential equation (expected y' or dy/dx)",
        lhs
    )))
}

/// `msolve <eq>; <eq>; …` — solve a system of linear equations exactly.
///
/// The equations may be given with an explicit `=` (`2x + y = 5`) or with the
/// right-hand side implicit (`2x + y - 5`). The solution is exact rational
/// arithmetic, not a float approximation.
pub(super) fn do_msolve(rest: &str) -> Result<Option<String>> {
    let solution = crate::linsolve::solve_system(rest.trim())?;
    let parts: Vec<String> = solution
        .iter()
        .map(|(v, r)| format!("{} = {}", v, r.to_string()))
        .collect();
    Ok(Some(parts.join(", ")))
}

pub(super) fn do_spline(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_pchip(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_minimize(rest: &str, ctx: &Context) -> Result<Option<String>> {
    // Format: "minimize <expr> [var] a b" or "minimize <expr> from a to b".
    // The bounds may be expressions ("pi", "pi/2").
    let bounds = split_bounds(rest)?;
    let a = eval_bound(&bounds.a_src, ctx)?;
    let b = eval_bound(&bounds.b_src, ctx)?;
    let mut var = "x".to_string();
    let mut tokens: Vec<&str> = bounds.expr_src.split_whitespace().collect();
    let mut expr_end = tokens.len();
    if expr_end > 1 {
        let candidate = tokens[expr_end - 1];
        if candidate.len() == 1 && candidate.chars().all(|c| c.is_ascii_alphabetic()) {
            var = candidate.to_string();
            expr_end -= 1;
        }
    }
    tokens.truncate(expr_end);
    let expr_src = tokens.join(" ");
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
    let mut obj = move |x: f64| -> f64 {
        call_ctx.set(var_c.clone(), x);
        eval(&e, &call_ctx).unwrap_or(f64::INFINITY)
    };
    // Domain probe: the function must be defined at both ends of the bracket.
    // Without this, an undefined endpoint is silently mapped to +∞ and the
    // search reports a meaningless "minimum" there.
    for x in [a, b] {
        let v = obj(x);
        if !v.is_finite() {
            return Err(crate::error::MathError::Domain(format!(
                "minimize: the function is not defined at {} = {} \
                 (it must be finite on the interval)",
                var,
                format_value(x)
            )));
        }
    }
    let (x_star, fx) = crate::optim::golden_section(obj, a, b, 1e-10)?;
    Ok(Some(format!(
        "{}* = {}\nf({}) = {}",
        var,
        format_value(x_star),
        format_value(x_star),
        format_value(fx)
    )))
}

