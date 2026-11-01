//! Export/import and alternate-representation command handlers.
//!
//! Handlers moved verbatim from the original monolithic `repl.rs`.

use super::*;
pub(super) fn do_mathml(rest: &str) -> Result<Option<String>> {
    if let Some(ml) = rest.strip_prefix("import ") {
        let e = crate::mathml::from_mathml(ml.trim())?;
        return Ok(Some(e.to_string()));
    }
    let e = Parser::parse(rest)?;
    Ok(Some(crate::mathml::to_mathml_doc(&e)))
}

/// `latex <expr>` → render the expression as LaTeX.
///
/// The counterpart of mathr's LaTeX *input* support: an expression written
/// in plain syntax comes back as TeX that can be pasted into a document.
pub(super) fn do_latex(rest: &str) -> Result<Option<String>> {
    let e = Parser::parse(rest)?;
    Ok(Some(e.to_tex()))
}

/// `lean <expr>` → export to a self-contained Lean 4 file (needs Mathlib)
pub(super) fn do_lean(rest: &str) -> Result<Option<String>> {    let e = Parser::parse(rest)?;
    Ok(Some(crate::lean::to_lean_doc(&e)?))
}

/// `serialize <fmt> <expr>` → export expression to `<fmt>` (sexpr | json | rpn)
/// `serialize <fmt> import <text>` → import from `<fmt>` and show the expression
pub(super) fn do_serialize(rest: &str) -> Result<Option<String>> {
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
pub(super) fn do_interval(rest: &str) -> Result<Option<String>> {
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

pub(super) fn do_cval(rest: &str) -> Result<Option<String>> {
    let z = crate::ceval::eval_complex_str(rest)?;
    Ok(Some(crate::ceval::format_complex(z)))
}

