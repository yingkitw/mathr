//! Symbolic differentiation.
//!
//! Given an [`Expr`] and the name of the variable to differentiate with
//! respect to, returns a derivative expression. The result is run through
//! [`simplify`] so identities like `x*0 = 0`, `-1*x = -x`, and
//! `x*1 = x` are applied.

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::simplify::simplify;

/// Differentiate `expr` with respect to `var`. Returns the simplified form
/// so the output is friendlier to read in a REPL.
pub fn differentiate(expr: &Expr, var: &str) -> Result<Expr> {
    Ok(simplify(&diff(expr, var)?))
}

fn diff(expr: &Expr, var: &str) -> Result<Expr> {
    match expr {
        Expr::Num(_) => Ok(Expr::num(0.0)),
        Expr::Var(v) => {
            if v == var {
                Ok(Expr::num(1.0))
            } else {
                Ok(Expr::num(0.0))
            }
        }
        Expr::Neg(e) => Ok(Expr::neg(diff(e, var)?)),
        Expr::Add(a, b) => Ok(Expr::add(diff(a, var)?, diff(b, var)?)),
        Expr::Sub(a, b) => Ok(Expr::sub(diff(a, var)?, diff(b, var)?)),
        Expr::Mul(a, b) => {
            // (fg)' = f'g + fg'
            let da = diff(a, var)?;
            let db = diff(b, var)?;
            Ok(Expr::add(
                Expr::mul(da, (**b).clone()),
                Expr::mul((**a).clone(), db),
            ))
        }
        Expr::Div(a, b) => {
            // (f/g)' = (f'g - fg') / g^2
            let da = diff(a, var)?;
            let db = diff(b, var)?;
            Ok(Expr::div(
                Expr::sub(
                    Expr::mul(da, (**b).clone()),
                    Expr::mul((**a).clone(), db),
                ),
                Expr::pow((**b).clone(), Expr::num(2.0)),
            ))
        }
        Expr::Pow(base, exp) => match (base.as_ref(), exp.as_ref()) {
            // f(x)^c
            (_, Expr::Num(c)) => {
                let c_val = *c;
                let inner = Expr::pow((**base).clone(), Expr::num(c_val - 1.0));
                Ok(Expr::mul(
                    Expr::mul(Expr::num(c_val), inner),
                    diff(base, var)?,
                ))
            }
            // c^g(x)
            (Expr::Num(_), _) => Ok(Expr::mul(
                expr.clone(),
                Expr::mul(
                    Expr::func("ln", vec![(**base).clone()]),
                    diff(exp, var)?,
                ),
            )),
            // f^g in general
            _ => Ok(Expr::mul(
                expr.clone(),
                Expr::add(
                    Expr::mul(diff(exp, var)?, Expr::func("ln", vec![(**base).clone()])),
                    Expr::mul(
                        Expr::div((**exp).clone(), (**base).clone()),
                        diff(base, var)?,
                    ),
                ),
            )),
        },
        Expr::Func(name, args) => {
            // Chain rule: d/dx f(g) = f'(g) * g'
            if args.len() != 1 {
                return Err(MathError::Eval(format!(
                    "cannot differentiate multi-arg function {}",
                    name
                )));
            }
            let arg = &args[0];
            let arg_diff = diff(arg, var)?;
            let deriv = derivative_of_builtin(name, arg)?;
            Ok(Expr::mul(deriv, arg_diff))
        }
    }
}

/// Return the derivative of the named elementary function applied to `arg`.
fn derivative_of_builtin(name: &str, arg: &Expr) -> Result<Expr> {
    let x = arg.clone();
    Ok(match name {
        "sin" => Expr::func("cos", vec![x]),
        "cos" => Expr::neg(Expr::func("sin", vec![x])),
        "tan" => Expr::div(
            Expr::num(1.0),
            Expr::pow(Expr::func("cos", vec![x]), Expr::num(2.0)),
        ),
        "asin" => Expr::div(
            Expr::num(1.0),
            Expr::func(
                "sqrt",
                vec![Expr::sub(Expr::num(1.0), Expr::pow(x, Expr::num(2.0)))],
            ),
        ),
        "acos" => Expr::neg(Expr::div(
            Expr::num(1.0),
            Expr::func(
                "sqrt",
                vec![Expr::sub(Expr::num(1.0), Expr::pow(x, Expr::num(2.0)))],
            ),
        )),
        "atan" => Expr::div(
            Expr::num(1.0),
            Expr::add(Expr::num(1.0), Expr::pow(x, Expr::num(2.0))),
        ),
        "sinh" => Expr::func("cosh", vec![x]),
        "cosh" => Expr::func("sinh", vec![x]),
        "tanh" => Expr::div(
            Expr::num(1.0),
            Expr::pow(Expr::func("cosh", vec![x]), Expr::num(2.0)),
        ),
        "exp" => Expr::func("exp", vec![x]),
        "ln" | "log" => Expr::div(Expr::num(1.0), x),
        "log2" => Expr::div(
            Expr::num(1.0),
            Expr::mul(x, Expr::func("ln", vec![Expr::num(2.0)])),
        ),
        "log10" => Expr::div(
            Expr::num(1.0),
            Expr::mul(x, Expr::func("ln", vec![Expr::num(10.0)])),
        ),
        "sqrt" => Expr::div(
            Expr::num(1.0),
            Expr::mul(Expr::num(2.0), Expr::func("sqrt", vec![x.clone()])),
        ),
        "abs" => Expr::div(x.clone(), Expr::func("abs", vec![x])),
        "floor" | "ceil" | "round" | "sign" | "fract" => Expr::num(0.0),
        _ => {
            return Err(MathError::Eval(format!(
                "no symbolic derivative for function '{}'",
                name
            )))
        }
    })
}

/// Compute the gradient of a multivariate expression.
///
/// Returns a vector of `(variable_name, partial_derivative)` pairs, one for
/// each free variable in `expr` (sorted alphabetically). Each partial
/// derivative is simplified.
pub fn gradient(expr: &Expr) -> Result<Vec<(String, Expr)>> {
    let vars = expr.variables();
    let mut result = Vec::with_capacity(vars.len());
    for v in &vars {
        result.push((v.clone(), differentiate(expr, v)?));
    }
    Ok(result)
}

/// Implicit differentiation for `F(indep, dep) = 0`.
///
/// Returns `d(dep)/d(indep) = −(∂F/∂indep) / (∂F/∂dep)`, simplified.
/// Errors when the two variable names coincide or when `∂F/∂dep` is
/// identically zero (vertical tangent / not a function of `indep`).
pub fn implicit_diff(expr: &Expr, indep: &str, dep: &str) -> Result<Expr> {
    if indep == dep {
        return Err(MathError::InvalidArgument(
            "implicit differentiation needs distinct independent and dependent variables".into(),
        ));
    }
    let f_indep = differentiate(expr, indep)?;
    let f_dep = differentiate(expr, dep)?;
    if matches!(&f_dep, Expr::Num(n) if *n == 0.0) {
        return Err(MathError::Eval(format!(
            "∂F/∂{dep} is zero — cannot solve for d{dep}/d{indep}"
        )));
    }
    // Cancel shared numeric coefficients so REPL output round-trips
    // (e.g. −(2x)/(2y) → −x/y rather than ambiguous 2*x/2*y).
    Ok(simplify(&cancel_numeric_ratio(Expr::neg(f_indep), f_dep)))
}

/// Build `num/den` after cancelling shared numeric coefficients peeled from
/// multiply/negate chains (local to idiff — does not change `simplify`).
fn cancel_numeric_ratio(num: Expr, den: Expr) -> Expr {
    let (cn, rn) = peel_numeric_coeff(&num);
    let (cd, rd) = peel_numeric_coeff(&den);
    if !cn.is_finite() || !cd.is_finite() || cd == 0.0 {
        return Expr::div(num, den);
    }
    let ratio = cn / cd;
    // Denominator is 1 after peeling → scaled numerator only.
    if matches!(&rd, Expr::Num(d) if (*d - 1.0).abs() < 1e-15) {
        return scale_by(ratio, rn);
    }
    // Fold the coefficient into the numerator: (ratio·rn)/rd.
    let numer = if matches!(&rn, Expr::Num(n) if (*n - 1.0).abs() < 1e-15) {
        Expr::num(ratio)
    } else {
        scale_by(ratio, rn)
    };
    Expr::div(numer, rd)
}

fn scale_by(ratio: f64, e: Expr) -> Expr {
    if (ratio - 1.0).abs() < 1e-12 {
        e
    } else if (ratio + 1.0).abs() < 1e-12 {
        Expr::neg(e)
    } else {
        Expr::mul(Expr::num(ratio), e)
    }
}

fn peel_numeric_coeff(e: &Expr) -> (f64, Expr) {
    match e {
        Expr::Num(n) => (*n, Expr::num(1.0)),
        Expr::Neg(inner) => {
            let (c, r) = peel_numeric_coeff(inner);
            (-c, r)
        }
        Expr::Mul(a, b) => match (a.as_ref(), b.as_ref()) {
            (Expr::Num(n), rest) => (*n, rest.clone()),
            (rest, Expr::Num(n)) => (*n, rest.clone()),
            _ => {
                let (ca, ra) = peel_numeric_coeff(a);
                let (cb, rb) = peel_numeric_coeff(b);
                (ca * cb, simplify(&Expr::mul(ra, rb)))
            }
        },
        other => (1.0, other.clone()),
    }
}

fn is_bound_name(v: &str) -> bool {
    matches!(v, "pi" | "e" | "tau" | "inf" | "i")
}

/// Free variables of `expr`, excluding built-in constants / the `i` literal.
fn free_vars(expr: &Expr) -> Vec<String> {
    expr.variables()
        .into_iter()
        .filter(|v| !is_bound_name(v))
        .collect()
}

/// Infer `(indep, dep)` for implicit differentiation.
///
/// Prefer `(x, other)` when `x` is present among exactly two free variables;
/// otherwise use alphabetical order. Callers with other conventions should
/// pass the pair explicitly.
pub fn infer_implicit_vars(expr: &Expr) -> Result<(String, String)> {
    let mut vars = free_vars(expr);
    vars.sort();
    vars.dedup();
    match vars.as_slice() {
        [a, b] if a == "x" => Ok((a.clone(), b.clone())),
        [a, b] if b == "x" => Ok((b.clone(), a.clone())),
        [a, b] => Ok((a.clone(), b.clone())),
        _ => Err(MathError::InvalidArgument(format!(
            "implicit differentiation needs exactly two free variables (found {:?}); pass <indep> <dep>",
            vars
        ))),
    }
}

/// REPL-facing steps for `idiff <expr> [= rhs] [<indep> <dep>]`.
pub fn idiff_steps_str(input: &str) -> Result<Vec<String>> {
    let (body, indep, dep) = parse_idiff_input(input)?;
    let f_indep = differentiate(&body, &indep)?;
    let f_dep = differentiate(&body, &dep)?;
    let result = implicit_diff(&body, &indep, &dep)?;
    Ok(vec![
        format!("F({}, {}) = {}", indep, dep, body),
        format!("∂F/∂{} = {}", indep, f_indep),
        format!("∂F/∂{} = {}", dep, f_dep),
        format!(
            "d{}/d{} = -(∂F/∂{})/(∂F/∂{}) = {}",
            dep, indep, indep, dep, result
        ),
    ])
}

/// Parse `idiff` input into `(F, indep, dep)`. Supports optional `= rhs`
/// (moved to the left) and optional trailing `<indep> <dep>`.
fn parse_idiff_input(input: &str) -> Result<(Expr, String, String)> {
    let input = input.trim();
    if input.is_empty() {
        return Err(MathError::Eval(
            "`idiff` needs: <expr> [= rhs] [<indep> <dep>]".into(),
        ));
    }

    let tokens: Vec<&str> = input.split_whitespace().collect();
    // Try trailing `<indep> <dep>` when the head still forms a valid relation.
    let (body, indep, dep) = if tokens.len() >= 3
        && is_ident(tokens[tokens.len() - 2])
        && is_ident(tokens[tokens.len() - 1])
    {
        let indep = tokens[tokens.len() - 2];
        let dep = tokens[tokens.len() - 1];
        let head = tokens[..tokens.len() - 2].join(" ");
        match parse_relation(&head) {
            Ok(body) => (body, indep.to_string(), dep.to_string()),
            Err(_) => {
                let body = parse_relation(input)?;
                let (a, b) = infer_implicit_vars(&body)?;
                (body, a, b)
            }
        }
    } else {
        let body = parse_relation(input)?;
        let (a, b) = infer_implicit_vars(&body)?;
        (body, a, b)
    };

    if indep == dep {
        return Err(MathError::InvalidArgument(
            "implicit differentiation needs distinct independent and dependent variables".into(),
        ));
    }
    Ok((body, indep, dep))
}

fn is_ident(s: &str) -> bool {
    let mut chars = s.chars();
    match chars.next() {
        Some(c) if c.is_ascii_alphabetic() || c == '_' => {
            chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
        }
        _ => false,
    }
}

fn parse_relation(src: &str) -> Result<Expr> {
    let src = src.trim();
    if src.is_empty() {
        return Err(MathError::Eval("empty expression".into()));
    }
    match src.split_once('=') {
        Some((lhs, rhs)) => {
            let l = crate::parser::Parser::parse(lhs.trim())?;
            let r = crate::parser::Parser::parse(rhs.trim())?;
            Ok(simplify(&Expr::sub(l, r)))
        }
        None => crate::parser::Parser::parse(src),
    }
}

/// Symbolic indefinite integration for the common elementary rules.
///
/// Handles:
/// - constants: ∫ c dx = c·x
/// - variable: ∫ x dx = x²/2
/// - powers: ∫ x^n dx = x^(n+1)/(n+1) for `n ≠ -1`, otherwise `ln(x)`
/// - sums / differences: linearity
/// - constant multiples: ∫ c·f dx = c·∫f dx
/// - exponentials: ∫ e^x dx = e^x,  ∫ a^x dx = a^x / ln(a)
/// - trig: ∫ sin(x) dx = −cos(x),  ∫ cos(x) dx = sin(x),  ∫ sec²(x) dx = tan(x)
/// - inverses: ∫ 1/x dx = ln(x),  ∫ 1/(1+x²) dx = atan(x),  ∫ 1/√(1−x²) dx = asin(x)
///
/// Returns `Err` for unsupported integrands (e.g., products of non-constants).
pub fn integrate(expr: &Expr, var: &str) -> Result<Expr> {
    Ok(simplify(&int_step(expr, var)?))
}

fn int_step(expr: &Expr, var: &str) -> Result<Expr> {
    match expr {
        Expr::Num(n) => Ok(Expr::mul(Expr::num(*n), Expr::var(var))),
        Expr::Var(v) if v == var => Ok(Expr::div(Expr::pow(Expr::var(var), Expr::num(2.0)), Expr::num(2.0))),
        Expr::Var(_) => Ok(Expr::mul(expr.clone(), Expr::var(var))),
        Expr::Neg(e) => Ok(Expr::neg(int_step(e, var)?)),
        Expr::Add(a, b) => Ok(Expr::add(int_step(a, var)?, int_step(b, var)?)),
        Expr::Sub(a, b) => Ok(Expr::sub(int_step(a, var)?, int_step(b, var)?)),
        Expr::Mul(a, b) => {
            if a.is_constant() {
                Ok(Expr::mul((**a).clone(), int_step(b, var)?))
            } else if b.is_constant() {
                Ok(Expr::mul(int_step(a, var)?, (**b).clone()))
            } else {
                Err(MathError::Eval(format!(
                    "integrate: cannot integrate non-linear product: {}",
                    expr
                )))
            }
        }
        Expr::Div(a, b) => {
            if b.is_constant() {
                Ok(Expr::mul(Expr::div(Expr::num(1.0), (**b).clone()), int_step(a, var)?))
            } else if a.is_constant() {
                integrate_constant_over(a, b, var)
            } else {
                Err(MathError::Eval(format!(
                    "integrate: cannot integrate non-constant numerator over non-constant denominator: {}",
                    expr
                )))
            }
        }
        Expr::Pow(base, exp) => match (base.as_ref(), exp.as_ref()) {
            (Expr::Var(v), Expr::Num(n)) if v == var && *n == -1.0 => {
                // ∫ x⁻¹ dx = ln|x| — abs makes it correct on both branches.
                Ok(Expr::func("ln", vec![Expr::func("abs", vec![Expr::var(var)])]))
            }
            (Expr::Var(v), Expr::Num(n)) if v == var => {
                let np1 = *n + 1.0;
                Ok(Expr::div(
                    Expr::pow(Expr::var(var), Expr::num(np1)),
                    Expr::num(np1),
                ))
            }
            // (a·var + b)ⁿ for numeric n: linear-substitution power rule.
            // Cross-port pattern from unsolve `symbolic.ts:331-342`.
            (base_e, Expr::Num(n)) if linear_in(base_e, var).is_some() => {
                if let Some((a, _)) = linear_in(base_e, var) {
                    if *n == -1.0 {
                        // ∫ (a·var+b)⁻¹ = ln|a·var+b| / a
                        Ok(Expr::div(
                            Expr::func("ln", vec![Expr::func("abs", vec![base_e.clone()])]),
                            Expr::num(a),
                        ))
                    } else {
                        let np1 = *n + 1.0;
                        Ok(Expr::div(
                            Expr::pow(base_e.clone(), Expr::num(np1)),
                            Expr::mul(Expr::num(np1), Expr::num(a)),
                        ))
                    }
                } else {
                    unreachable!("guard above ensures success")
                }
            }
            (Expr::Num(_), _) | (_, _) if (**base).is_constant() => {
                Ok(Expr::div(
                    expr.clone(),
                    Expr::func("ln", vec![(**base).clone()]),
                ))
            }
            _ => Err(MathError::Eval(format!(
                "integrate: cannot integrate power: {}",
                expr
            ))),
        },
        Expr::Func(name, args) if args.len() == 1 => {
            let inner = &args[0];
            // Linear-substitution fast path.  For inner = a·var + b (with a, b
            // constant in `var` and a != 0), any table integrand f(inner)
            // integrates as F(inner)/a.  This subsumes the bare-Var cases
            // below and additionally covers `sin(3*x+1)`, `exp(2*x)`,
            // `(x+1)^n`, `1/(2*x+1) → ln|2*x+1|/2`, `sqrt(2*x+1)`, etc.
            // (Cross-port pattern from unsolve `symbolic.ts:228-391`.)
            if let Some((a, b)) = linear_in(inner, var) {
                let shifted = shift_var(var, a, b);
                match name.as_str() {
                    "exp" => return Ok(Expr::div(
                        Expr::func("exp", vec![shifted]),
                        Expr::num(a),
                    )),
                    "sin" => return Ok(Expr::div(Expr::neg(shifted.cos()), Expr::num(a))),
                    "cos" => return Ok(Expr::div(shifted.sin(), Expr::num(a))),
                    "sinh" => return Ok(Expr::div(shifted.cosh(), Expr::num(a))),
                    "cosh" => return Ok(Expr::div(shifted.sinh(), Expr::num(a))),
                    "tan" => return Ok(Expr::neg(Expr::div(
                        Expr::func("ln", vec![Expr::func("abs", vec![shifted.cos()])]),
                        Expr::num(a),
                    ))),
                    // ∫ sqrt(a·x + b) dx = (2/3a) · (a·x + b)^(3/2)
                    "sqrt" => return Ok(Expr::div(
                        Expr::mul(
                            Expr::num(2.0),
                            Expr::pow(shifted, Expr::num(1.5)),
                        ),
                        Expr::mul(Expr::num(3.0), Expr::num(a)),
                    )),
                    // ∫ cbrt(a·x + b) dx = (3/4a) · (a·x + b)^(4/3)
                    "cbrt" => return Ok(Expr::div(
                        Expr::mul(
                            Expr::num(3.0),
                            Expr::pow(shifted, Expr::num(4.0 / 3.0)),
                        ),
                        Expr::mul(Expr::num(4.0), Expr::num(a)),
                    )),
                    _ => {}
                }
            }
            let inner_is_var = matches!(inner, Expr::Var(v) if v == var);
            match (name.as_str(), inner_is_var) {
                ("exp", true) => Ok(Expr::func("exp", vec![Expr::var(var)])),
                ("sin", true) => Ok(Expr::neg(Expr::func("cos", vec![Expr::var(var)]))),
                ("cos", true) => Ok(Expr::func("sin", vec![Expr::var(var)])),
                ("sec", true) => Ok(Expr::func("ln", vec![Expr::add(
                    Expr::func("sec", vec![Expr::var(var)]),
                    Expr::func("tan", vec![Expr::var(var)]),
                )])),
                ("tan", true) => Ok(Expr::neg(Expr::func("ln", vec![Expr::func("cos", vec![Expr::var(var)])]))),
                _ => Err(MathError::Eval(format!(
                    "integrate: unsupported integrand `{}{}`", name, inner
                ))),
            }
        }
        Expr::Func(name, args) => Err(MathError::Eval(format!(
            "integrate: cannot integrate multi-argument function {} with {} args",
            name,
            args.len()
        ))),
    }
}

/// Decompose an expression as `a*var + b` with `a, b` constant in `var`.
///
/// Returns `Some((a, b))` only when the expression is *affine-linear* in
/// `var` with a non-zero slope.  `Some((1.0, 0.0))` for plain `Var(v)`,
/// `Some((2.0, 1.0))` for `2*x + 1`, `Some((-1.0, 3.0))` for `3 - x`, etc.
/// Returns `None` for constants, non-linear polynomials, function calls,
/// or anything else that isn't a clean linear-in-`var` form.
fn linear_in(e: &Expr, var: &str) -> Option<(f64, f64)> {
    match e {
        Expr::Var(v) if v == var => Some((1.0, 0.0)),
        Expr::Neg(inner) => {
            let (a, b) = linear_in(inner, var)?;
            Some((-a, -b))
        }
        // Pure constants are a=0, b=c (degenerate linear).
        Expr::Num(c) => Some((0.0, *c)),
        Expr::Mul(a, b) => {
            // c·var or var·c (with c constant)
            let (lhs, rhs) = (a.as_ref(), b.as_ref());
            if let (Expr::Var(v), Expr::Num(c)) = (lhs, rhs) {
                if v == var {
                    return Some((*c, 0.0));
                }
            }
            if let (Expr::Num(c), Expr::Var(v)) = (lhs, rhs) {
                if v == var {
                    return Some((*c, 0.0));
                }
            }
            None
        }
        Expr::Add(a, b) => {
            let (a1, b1) = linear_in(a, var)?;
            let (a2, b2) = linear_in(b, var)?;
            // Two coincident var terms = quadratic → not linear.
            if a1 != 0.0 && a2 != 0.0 {
                return None;
            }
            Some((a1 + a2, b1 + b2))
        }
        Expr::Sub(a, b) => {
            let (a1, b1) = linear_in(a, var)?;
            let (a2, b2) = linear_in(b, var)?;
            if a1 != 0.0 && a2 != 0.0 {
                return None;
            }
            Some((a1 - a2, b1 - b2))
        }
        _ => None,
    }
}

/// Reconstruct the shifted argument `a*var + b` from linear decomposition.
/// Caller must have already verified `linear_in` succeeded.
fn shift_var(var: &str, a: f64, b: f64) -> Expr {
    if (a - 1.0).abs() < 1e-15 && b == 0.0 {
        Expr::var(var)
    } else if (a - 1.0).abs() < 1e-15 {
        Expr::add(Expr::var(var), Expr::num(b))
    } else if b == 0.0 {
        Expr::mul(Expr::num(a), Expr::var(var))
    } else {
        Expr::add(Expr::mul(Expr::num(a), Expr::var(var)), Expr::num(b))
    }
}

/// Trait extension for building the function-call Exprs in the linear-substitution
/// fast path. Keeps the call sites compact.
trait FuncExpr {
    fn cos(self) -> Expr;
    fn sin(self) -> Expr;
    fn cosh(self) -> Expr;
    fn sinh(self) -> Expr;
}
impl FuncExpr for Expr {
    fn cos(self) -> Expr { Expr::func("cos", vec![self]) }
    fn sin(self) -> Expr { Expr::func("sin", vec![self]) }
    fn cosh(self) -> Expr { Expr::func("cosh", vec![self]) }
    fn sinh(self) -> Expr { Expr::func("sinh", vec![self]) }
}

/// Handle integrands of the form `c / g(x)` where `c` is constant.
/// Recognises `1/x` and `1/(1+x²)` and `1/√(1−x²)`.
fn integrate_constant_over(num: &Expr, den: &Expr, var: &str) -> Result<Expr> {
    let num_val = if let Expr::Num(n) = num { *n } else { 1.0 };
    let _ = num_val;
    // c / (a·var + b) → (c/a) · ln|a·var+b|.  Linear-substitution pattern.
    if let Some((a, _b)) = linear_in(den, var) {
        return Ok(Expr::mul(
            Expr::num(num_val / a),
            Expr::func("ln", vec![Expr::func("abs", vec![den.clone()])]),
        ));
    }
    match den {
        Expr::Var(v) if v == var => Ok(Expr::mul(Expr::num(num_val), Expr::func("ln", vec![Expr::var(var)]))),
        Expr::Add(a, b) | Expr::Sub(a, b) => {
            // 1 / (1 ± x²) → atan or -atan
            let is_one = matches!(a.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
            let is_xsq = match b.as_ref() {
                Expr::Pow(p, e) => matches!(p.as_ref(), Expr::Var(v) if v == var)
                    && matches!(e.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12),
                _ => false,
            };
            if is_one && is_xsq {
                let sign = if matches!(den, Expr::Sub(..)) { -1.0 } else { 1.0 };
                Ok(Expr::mul(
                    Expr::num(sign * num_val),
                    Expr::func("atan", vec![Expr::var(var)]),
                ))
            } else {
                Err(MathError::Eval(format!(
                    "integrate: unsupported integrand 1/{}", den
                )))
            }
        }
        Expr::Func(name, args) if name == "sqrt" && args.len() == 1 => {
            // 1 / sqrt(1 - x²) → asin
            if let Expr::Sub(a_inner, b_inner) = &args[0] {
                let is_one = matches!(a_inner.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
                let is_xsq = matches!(
                    b_inner.as_ref(),
                    Expr::Pow(pp, ee)
                        if matches!(pp.as_ref(), Expr::Var(v) if v == var)
                            && matches!(ee.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12)
                );
                if is_one && is_xsq {
                    return Ok(Expr::mul(
                        Expr::num(num_val),
                        Expr::func("asin", vec![Expr::var(var)]),
                    ));
                }
            }
            Err(MathError::Eval(format!(
                "integrate: unsupported integrand 1/{}",
                den
            )))
        }
        Expr::Pow(p, e) => {
            // 1 / (1 - x²)^{1/2} (as Pow rather than sqrt) → asin
            if let Expr::Sub(a_inner, b_inner) = p.as_ref() {
                let is_one = matches!(a_inner.as_ref(), Expr::Num(n) if (*n - 1.0).abs() < 1e-12);
                let is_xsq = matches!(
                    b_inner.as_ref(),
                    Expr::Pow(pp, ee)
                        if matches!(pp.as_ref(), Expr::Var(v) if v == var)
                            && matches!(ee.as_ref(), Expr::Num(n) if (*n - 2.0).abs() < 1e-12)
                );
                let is_half = matches!(e.as_ref(), Expr::Num(n) if (*n - 0.5).abs() < 1e-12);
                if is_one && is_xsq && is_half {
                    return Ok(Expr::mul(
                        Expr::num(num_val),
                        Expr::func("asin", vec![Expr::var(var)]),
                    ));
                }
            }
            Err(MathError::Eval(format!(
                "integrate: unsupported integrand 1/{}",
                den
            )))
        }
        _ => Err(MathError::Eval(format!(
            "integrate: unsupported integrand {}/{}", num, den
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eval::{eval, Context};
    use crate::parser::Parser;

    /// Compare two expressions by evaluating them at the supplied `x` values.
    /// Simpler and more robust than AST equality once you go past the simplest
    /// of cases — the symbolic output is rarely byte-identical to the "nice"
    /// textbook form.
    fn agrees(got: &Expr, want: &Expr, xs: &[f64]) {
        let mut ctx = Context::standard();
        for &x in xs {
            ctx.set("x", x);
            let g = eval(got, &ctx).unwrap();
            let w = eval(want, &ctx).unwrap();
            assert!(
                (g - w).abs() < 1e-9,
                "disagree at x={}: got {} want {}",
                x,
                g,
                w
            );
        }
    }

    fn d(s: &str) -> Expr {
        let e = Parser::parse(s).unwrap();
        differentiate(&e, "x").unwrap()
    }

    #[test]
    fn polynomial() {
        // d/dx (x^3 + 2x^2 + x + 5) = 3x^2 + 4x + 1
        let got = d("x^3 + 2*x^2 + x + 5");
        let want = Parser::parse("3*x^2 + 4*x + 1").unwrap();
        agrees(&got, &want, &[0.0, 1.0, -2.0, 3.5, 0.7]);
    }

    #[test]
    fn product_rule() {
        // d/dx (x * sin(x)) = sin(x) + x*cos(x)
        let got = d("x*sin(x)");
        let want = Parser::parse("sin(x) + x*cos(x)").unwrap();
        agrees(&got, &want, &[0.1, 0.5, 1.0, 2.0, -0.7]);
    }

    #[test]
    fn quotient_rule() {
        // d/dx (x / (x+1)) = 1/(x+1)^2
        let got = d("x/(x+1)");
        let want = Parser::parse("1/(x+1)^2").unwrap();
        agrees(&got, &want, &[0.5, 1.5, 2.0, -0.5]);
    }

    #[test]
    fn chain_rule() {
        // d/dx sin(x^2) = 2*x*cos(x^2)
        let got = d("sin(x^2)");
        let want = Parser::parse("2*x*cos(x^2)").unwrap();
        agrees(&got, &want, &[0.1, 0.5, 1.0, -0.7]);
    }

    #[test]
    fn exp_ln() {
        // d/dx exp(x) = exp(x)
        let got = d("exp(x)");
        let want = Parser::parse("exp(x)").unwrap();
        agrees(&got, &want, &[0.0, 1.0, -2.0]);
    }

    fn integrate_agrees(src: &str, want: &str, xs: &[f64]) {
        // Differentiate the symbolic integral and verify it matches the
        // original integrand at sample points.
        let e = Parser::parse(src).unwrap();
        let antideriv = integrate(&e, "x").unwrap();
        let derived = differentiate(&antideriv, "x").unwrap();
        let want_e = Parser::parse(want).unwrap();
        agrees(&derived, &want_e, xs);
    }

    #[test]
    fn integrate_constant() {
        // ∫ 3 dx = 3x
        let e = Parser::parse("3").unwrap();
        let result = integrate(&e, "x").unwrap();
        let want = Parser::parse("3*x").unwrap();
        agrees(&result, &want, &[1.0, 2.0, -5.0]);
    }

    #[test]
    fn integrate_polynomial() {
        integrate_agrees("x", "x", &[0.5, 1.0, 2.0]);
        integrate_agrees("x^2", "x^2", &[0.5, 1.0, 2.0]);
        integrate_agrees("x^3 - 2*x + 1", "x^3 - 2*x + 1", &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn integrate_reciprocal() {
        // ∫ 1/x dx = ln(x)
        integrate_agrees("1/x", "1/x", &[0.5, 1.5, 3.0]);
    }

#[test]
    fn integrate_exp() {
        integrate_agrees("2^x", "2^x", &[0.0, 1.0, 2.0]);
    }

    // ----- Linear-substitution integration table (cross-port from
    // unsolve `symbolic.ts:228-391`).  For inner = a·x + b, ∫f(inner)dx = F(inner)/a.
    // Each test verifies via d/dx(antideriv) = integrand at sample points.

    #[test]
    fn integrate_linear_subst_sin() {
        // ∫ sin(3x + 1) dx = −cos(3x + 1)/3
        integrate_agrees("sin(3*x + 1)", "sin(3*x + 1)", &[-0.3, 0.1, 0.5, 0.9]);
    }

    #[test]
    fn integrate_linear_subst_cos() {
        integrate_agrees("cos(2*x - 1)", "cos(2*x - 1)", &[-0.5, 0.0, 0.5, 1.2]);
    }

    #[test]
    fn integrate_linear_subst_exp() {
        integrate_agrees("exp(2*x)", "exp(2*x)", &[-0.5, 0.0, 0.5, 1.0]);
    }

    #[test]
    fn integrate_linear_subst_sinh() {
        integrate_agrees("sinh(2*x)", "sinh(2*x)", &[-0.5, 0.0, 0.5, 1.0]);
    }

    #[test]
    fn integrate_linear_subst_cosh() {
        integrate_agrees("cosh(x + 1)", "cosh(x + 1)", &[-0.7, 0.0, 0.5, 1.5]);
    }

    #[test]
    fn integrate_linear_subst_power() {
        // ∫ (x + 1)^3 dx = (x + 1)^4 / 4
        integrate_agrees("(x + 1)^3", "(x + 1)^3", &[-0.5, 0.0, 0.5, 1.5]);
    }

    #[test]
    fn integrate_linear_subst_power_negative() {
        // ∫ (x + 1)^(-1) dx = ln|x + 1|
        integrate_agrees("1/(x + 1)", "1/(x + 1)", &[0.0, 0.5, 1.5, 2.5]);
    }

    #[test]
    fn integrate_linear_subst_reciprocal_with_coefficient() {
        // ∫ 5/(2*x + 3) dx = (5/2)·ln|2*x + 3|
        integrate_agrees("5/(2*x + 3)", "5/(2*x + 3)", &[-1.0, 0.0, 0.5, 1.5]);
    }

    #[test]
    fn integrate_linear_subst_reciprocal_negative_var() {
        // ∫ 1/(1 - x) dx = −ln|1 - x|
        integrate_agrees("1/(1 - x)", "1/(1 - x)", &[-1.0, -0.5, 0.5, 0.9]);
    }

    #[test]
    fn integrate_linear_subst_ln_x_now_uses_abs() {
        // ∫ 1/x dx = ln|x| — verify by evaluating at a NEGATIVE x that
        // the new abs-wrapped form matches the integrand.  Plain ln(x)
        // would silently disagree there (returns NaN for x<0).
        let e = Parser::parse("1/x").unwrap();
        let antideriv = integrate(&e, "x").unwrap();
        let derived = differentiate(&antideriv, "x").unwrap();
        let want_e = Parser::parse("1/x").unwrap();
        // Use a non-zero mix including negative x.
        agrees(&derived, &want_e, &[-2.0, -0.5, 0.5, 1.0, 3.0]);
    }

    #[test]
    fn integrate_linear_subst_quarter() {
        // ∫ (2x)^4 dx = (2x)^5 / 10.  D/dx should give (2x)^4.
        integrate_agrees("(2*x)^4", "(2*x)^4", &[-0.5, 0.0, 0.5, 1.0]);
    }

    #[test]
    fn integrate_linear_subst_sqrt() {
        // sqrt(2*x + 1) parses as Func("sqrt", ...) — verify it's recognised
        // via the linear-substitution path in the Func arm.
        // ∫ sqrt(2x+1) dx = (2x+1)^(3/2) / 3  (since chain rule gives ·2/2=1).
        let e = Parser::parse("sqrt(2*x + 1)").unwrap();
        let antideriv = integrate(&e, "x").unwrap();
        let derived = differentiate(&antideriv, "x").unwrap();
        // Verify at x values where 2x+1 > 0
        let mut ctx = Context::standard();
        for &x in &[0.0, 0.5, 1.0, 2.5] {
            ctx.set("x", x);
            let d = eval(&derived, &ctx).unwrap();
            let expected = (2.0 * x + 1.0).sqrt();
            assert!((d - expected).abs() < 1e-9, "at x={x}: got {d} want {expected}");
        }
    }

    #[test]
    fn integrate_trig() {
        integrate_agrees("sin(x)", "sin(x)", &[0.5, 1.0, 2.0]);
        integrate_agrees("cos(x)", "cos(x)", &[0.5, 1.0, 2.0]);
    }

    #[test]
    fn integrate_atan_arcsin() {
        integrate_agrees("1/(1+x^2)", "1/(1+x^2)", &[0.5, 1.0, 2.0]);
        integrate_agrees("1/sqrt(1-x^2)", "1/sqrt(1-x^2)", &[0.0, 0.3, 0.5]);
    }

    #[test]
    fn partial_derivative() {
        // ∂/∂x (x^2 * y + y^3) = 2*x*y
        let e = Parser::parse("x^2 * y + y^3").unwrap();
        let got = differentiate(&e, "x").unwrap();
        let want = Parser::parse("2*x*y").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at x={},y={}: got {} want {}", x, y, g, w);
        }
    }

    #[test]
    fn partial_derivative_other_var() {
        // ∂/∂y (x^2 * y + y^3) = x^2 + 3*y^2
        let e = Parser::parse("x^2 * y + y^3").unwrap();
        let got = differentiate(&e, "y").unwrap();
        let want = Parser::parse("x^2 + 3*y^2").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at x={},y={}: got {} want {}", x, y, g, w);
        }
    }

    #[test]
    fn gradient_multivariate() {
        // ∇(x^2 + x*y + y^2) = [∂/∂x = 2*x + y, ∂/∂y = x + 2*y]
        let e = Parser::parse("x^2 + x*y + y^2").unwrap();
        let grad = gradient(&e).unwrap();
        assert_eq!(grad.len(), 2);
        assert_eq!(grad[0].0, "x");
        assert_eq!(grad[1].0, "y");

        let want_dx = Parser::parse("2*x + y").unwrap();
        let want_dy = Parser::parse("x + 2*y").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(1.0, 2.0), (0.5, -1.0), (3.0, 0.7)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let gdx = eval(&grad[0].1, &ctx).unwrap();
            let wdx = eval(&want_dx, &ctx).unwrap();
            assert!((gdx - wdx).abs() < 1e-9, "dx at x={},y={}: got {} want {}", x, y, gdx, wdx);
            let gdy = eval(&grad[1].1, &ctx).unwrap();
            let wdy = eval(&want_dy, &ctx).unwrap();
            assert!((gdy - wdy).abs() < 1e-9, "dy at x={},y={}: got {} want {}", x, y, gdy, wdy);
        }
    }

    #[test]
    fn gradient_three_vars() {
        // ∇(x*y*z) = [∂/∂x = y*z, ∂/∂y = x*z, ∂/∂z = x*y]
        let e = Parser::parse("x*y*z").unwrap();
        let grad = gradient(&e).unwrap();
        assert_eq!(grad.len(), 3);

        let wants = [
            Parser::parse("y*z").unwrap(),
            Parser::parse("x*z").unwrap(),
            Parser::parse("x*y").unwrap(),
        ];
        let mut ctx = Context::standard();
        for &(x, y, z) in &[(1.0, 2.0, 3.0), (0.5, -1.0, 2.0)] {
            ctx.set("x", x);
            ctx.set("y", y);
            ctx.set("z", z);
            for (i, want) in wants.iter().enumerate() {
                let g = eval(&grad[i].1, &ctx).unwrap();
                let w = eval(want, &ctx).unwrap();
                assert!((g - w).abs() < 1e-9, "var {} at ({},{},{}): got {} want {}", grad[i].0, x, y, z, g, w);
            }
        }
    }

    #[test]
    fn gradient_constant() {
        // ∇(42) = [] (no variables)
        let e = Parser::parse("42").unwrap();
        let grad = gradient(&e).unwrap();
        assert!(grad.is_empty());
    }

    #[test]
    fn implicit_circle() {
        // x² + y² − 1 = 0 → dy/dx = −x/y
        let e = Parser::parse("x^2 + y^2 - 1").unwrap();
        let got = implicit_diff(&e, "x", "y").unwrap();
        let want = Parser::parse("-x/y").unwrap();
        let mut ctx = Context::standard();
        for &(x, y) in &[(0.6, 0.8), (-0.3, 0.5), (0.1, -0.9)] {
            ctx.set("x", x);
            ctx.set("y", y);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at ({},{}): got {} want {}", x, y, g, w);
        }
    }

    #[test]
    fn implicit_infer_xy() {
        let e = Parser::parse("x^2 + y^2 - 1").unwrap();
        let (indep, dep) = infer_implicit_vars(&e).unwrap();
        assert_eq!((indep.as_str(), dep.as_str()), ("x", "y"));
    }

    #[test]
    fn implicit_vertical_rejected() {
        // F = x − 1 → ∂F/∂y = 0
        let e = Parser::parse("x - 1").unwrap();
        assert!(implicit_diff(&e, "x", "y").is_err());
    }

    #[test]
    fn implicit_same_var_rejected() {
        let e = Parser::parse("x^2 + y").unwrap();
        assert!(implicit_diff(&e, "x", "x").is_err());
    }

    #[test]
    fn idiff_steps_equation_form() {
        let steps = idiff_steps_str("x^2 + y^2 = 1").unwrap();
        assert_eq!(steps.len(), 4);
        assert!(steps[0].contains("F(x, y)"), "{:?}", steps[0]);
        assert!(steps[3].starts_with("dy/dx ="), "last step: {}", steps[3]);
        assert!(steps[3].contains("/"), "expected a quotient in {}", steps[3]);
    }

    #[test]
    fn idiff_steps_explicit_vars() {
        // Treat t as independent, s as dependent: s^2 − t = 0 → ds/dt = 1/(2s)
        let steps = idiff_steps_str("s^2 - t t s").unwrap();
        assert!(steps[0].contains("F(t, s)"), "{:?}", steps[0]);
        assert!(steps[3].starts_with("ds/dt ="), "last step: {}", steps[3]);
        let got = implicit_diff(&Parser::parse("s^2 - t").unwrap(), "t", "s").unwrap();
        let want = Parser::parse("1/(2*s)").unwrap();
        let mut ctx = Context::standard();
        for &s in &[0.5, 1.0, 2.0] {
            ctx.set("s", s);
            ctx.set("t", s * s);
            let g = eval(&got, &ctx).unwrap();
            let w = eval(&want, &ctx).unwrap();
            assert!((g - w).abs() < 1e-9, "at s={}: got {} want {}", s, g, w);
        }
    }
}