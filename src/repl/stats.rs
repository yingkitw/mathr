//! Statistics and hypothesis-testing command handlers.
//!
//! Handlers moved verbatim from the original monolithic `repl.rs`.

use super::*;
pub(super) fn do_sample(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_dist(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_fit(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_mwu(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_wilcoxon(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_kw(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_boot(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_logit(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_qtile(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_ttest(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_chitest(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_anova(rest: &str) -> Result<Option<String>> {
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
pub(super) fn format_test_result(r: &crate::dists::TestResult) -> String {
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

pub(super) fn do_stats(rest: &str) -> Result<Option<String>> {
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

