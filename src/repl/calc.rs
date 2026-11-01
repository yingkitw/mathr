//! Calculus command handlers.
//!
//! Handlers moved verbatim from the original monolithic `repl.rs`.

use super::*;
pub(super) fn do_diff(rest: &str, force_wrt: Option<&str>, ctx: &mut Context) -> Result<Option<String>> {
    let parse = parse_diff_tail(rest, force_wrt)?;
    let e = diff_expr_prepared(&parse, ctx)?;
    let s = if parse.order == 1 {
        differentiate(&e, &parse.wrt)?
    } else {
        crate::symbolic::differentiate_n(&e, &parse.wrt, parse.order)?
    };
    Ok(Some(s.to_string()))
}

pub(super) fn do_pdiff(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
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

pub(super) fn do_gradient(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
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

pub(super) fn do_idiff(rest: &str) -> Result<Option<String>> {
    let steps = crate::symbolic::idiff_steps_str(rest)?;
    // Final step: `dy/dx = -(∂F/∂x)/(∂F/∂y) = <expr>` — return the result expr.
    let last = steps.last().cloned().unwrap_or_default();
    Ok(Some(
        last.rsplit_once(" = ")
            .map(|(_, r)| r.trim().to_string())
            .unwrap_or(last),
    ))
}

pub(super) fn do_integrate(rest: &str, ctx: &mut Context) -> Result<Option<String>> {
    // Format: "<expr> <a> <b>" or "<expr> from <a> to <b>". The bounds may be
    // expressions ("pi", "pi/2").
    let bounds = split_bounds(rest)?;
    let a = eval_bound(&bounds.a_src, ctx)?;
    let b = eval_bound(&bounds.b_src, ctx)?;
    let e = Parser::parse(&bounds.expr_src)?;

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

pub(super) fn do_taylor(rest: &str) -> Result<Option<String>> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
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
    // An explicit `for <var>` names the expansion variable.
    let wrt = variable.unwrap_or_else(|| "x".to_string());
    let series = crate::taylor::taylor_series_str(&expr_src, &wrt, around, order)?;
    Ok(Some(series.to_string()))
}

pub(super) fn do_laurent(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_rat(rest: &str) -> Result<Option<String>> {
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

pub(super) fn parse_limit_target(rest: &str) -> Result<(crate::expr::Expr, String, f64)> {
    let ForVar { body, variable } = strip_for_var(rest);
    let tokens: Vec<&str> = body.split_whitespace().collect();
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
        // The `at` keyword separates the expression from the point:
        // `limit sin(x)/x at 0`. When present, the variable is inferred
        // from the expression instead of being read positionally.
        let has_at = tokens[tokens.len() - 2].eq_ignore_ascii_case("at");
        let var = match &variable {
            Some(v) => v.clone(),
            None if has_at => {
                // Infer the variable from the expression: exactly one free
                // variable is required.
                let e = Parser::parse(&tokens[..tokens.len() - 2].join(" "))?;
                let mut vars = e.variables();
                vars.sort();
                vars.dedup();
                if vars.len() != 1 {
                    return Err(crate::error::MathError::Eval(
                        "limit needs exactly one variable (found: none); use `limit <expr> <var> <point>`"
                            .into(),
                    ));
                }
                vars.remove(0)
            }
            None => tokens[tokens.len() - 2].to_string(),
        };
        // The expression ends before the point, and also before the
        // variable slot (which may be the literal `at` or a name).
        let expr_end = if has_at {
            tokens.len() - 2
        } else if variable.is_some() {
            // An explicit `for <var>`: the variable was not a positional
            // token, so only the point is excluded.
            tokens.len() - 1
        } else {
            tokens.len() - 2
        };
        let expr_src = tokens[..expr_end].join(" ");
        let e = Parser::parse(&expr_src)?;
        Ok((e, var, point))
    } else if tokens.len() == 2 {
        let point = parse_point(tokens[1])?;
        let e = Parser::parse(tokens[0])?;
        let var = match variable {
            Some(v) => v,
            None => {
                let mut vars = e.variables();
                if vars.len() != 1 {
                    return Err(crate::error::MathError::Eval(
                        "limit needs: <expr> <var> <point> (expression has multiple or no variables)"
                            .into(),
                    ));
                }
                vars.remove(0)
            }
        };
        Ok((e, var, point))
    } else {
        Err(crate::error::MathError::Eval(
            "limit needs: <expr> [<var>] <point>  e.g. `limit sin(x)/x x 0` or `limit 1/x inf`"
                .into(),
        ))
    }
}

pub(super) fn do_limit(rest: &str) -> Result<Option<String>> {
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
pub(super) fn limit_steps_str(rest: &str) -> Result<Vec<String>> {
    let (e, var, point) = parse_limit_target(rest)?;
    crate::limit::limit_steps(&e, &var, point)
}

/// `apart <expr> [var]` — partial fraction decomposition of a rational
/// function. The variable may be omitted when the expression has exactly
/// one.
pub(super) fn do_apart(rest: &str) -> Result<Option<String>> {
    let (e, var) = parse_single_var_target(rest)?;
    let out = crate::apart::apart(&e, &var)?;
    Ok(Some(out.to_string()))
}

/// Parse `"<expr> [var]"` — infer the variable from the expression when
/// omitted (requires exactly one free variable).
pub(super) fn parse_single_var_target(rest: &str) -> Result<(crate::expr::Expr, String)> {
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
/// Σ/∏ over integer bounds: `<expr> [<var>] <a> <b>` with bounds from the end.
pub(super) fn do_sumprod(rest: &str, is_sum: bool) -> Result<Option<String>> {
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
pub(super) fn do_ad(rest: &str, ctx: &Context) -> Result<Option<String>> {
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

pub(super) fn do_romberg(rest: &str) -> Result<Option<String>> {
    // Format: "<expr> a b" or "<expr> from a to b" (fixed 8 Romberg levels).
    // The bounds may be expressions ("pi", "pi/2").
    let bounds = split_bounds(rest)?;
    let ctx = Context::standard();
    let a = eval_bound(&bounds.a_src, &ctx)?;
    let b = eval_bound(&bounds.b_src, &ctx)?;
    let levels: usize = 8;
    let e = Parser::parse(&bounds.expr_src)?;
    let ctx2 = Context::standard();
    let f = move |x: f64| {
        let mut cx = ctx2.clone();
        cx.set("x", x);
        crate::eval::eval(&e, &cx).unwrap_or(f64::NAN)
    };
    let v = crate::calculus::integrate_romberg(f, a, b, levels)?;
    Ok(Some(format!(
        "∫ ({}, {}; levels = {}) = {}",
        format_value(a),
        format_value(b),
        levels,
        format_value(v)
    )))
}

pub(super) fn do_integrate_sym(rest: &str, _ctx: &mut Context) -> Result<Option<String>> {
    // Format: "<expr> [var]" — symbolic integration.
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
    Ok(Some(format!("∫ {} d{} = {}", expr_src, var, result)))
}

