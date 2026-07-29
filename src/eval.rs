use std::collections::HashMap;

use crate::error::{MathError, Result};
use crate::expr::Expr;

/// A binding of variables and functions used during evaluation.
#[derive(Debug, Clone, Default)]
pub struct Context {
    pub vars: HashMap<String, f64>,
    pub funcs: HashMap<String, Func>,
}

/// A built-in or user-defined single- or multi-argument math function.
#[derive(Debug, Clone)]
pub enum Func {
    Builtin(fn(&[f64]) -> Result<f64>),
    /// A user-supplied closure (used when an expression is assigned to a name).
    User(Expr, Vec<String>),
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }

    /// Standard math context: constants `pi`, `e`, `tau`, `inf`, and the full
    /// set of elementary functions (`sin`, `cos`, `tan`, `exp`, `log`, ...).
    pub fn standard() -> Self {
        let mut ctx = Self::new();
        ctx.vars.insert("pi".into(), std::f64::consts::PI);
        ctx.vars.insert("e".into(), std::f64::consts::E);
        ctx.vars.insert("tau".into(), std::f64::consts::TAU);
        for (name, f) in builtins() {
            ctx.funcs.insert(name.into(), Func::Builtin(f));
        }
        ctx
    }

    pub fn set<S: Into<String>>(&mut self, name: S, value: f64) {
        self.vars.insert(name.into(), value);
    }

    pub fn with<S: Into<String>>(mut self, name: S, value: f64) -> Self {
        self.set(name, value);
        self
    }

    pub fn define<S: Into<String>>(&mut self, name: S, expr: Expr, params: Vec<String>) {
        self.funcs.insert(name.into(), Func::User(expr, params));
    }
}

/// Evaluate `expr` in the supplied `ctx`, returning `f64`.
pub fn eval(expr: &Expr, ctx: &Context) -> Result<f64> {
    match expr {
        Expr::Num(n) => Ok(*n),
        Expr::Var(name) => ctx
            .vars
            .get(name)
            .copied()
            .ok_or_else(|| MathError::UnknownVariable(name.clone())),
        Expr::Neg(e) => Ok(-eval(e, ctx)?),
        Expr::Add(a, b) => Ok(eval(a, ctx)? + eval(b, ctx)?),
        Expr::Sub(a, b) => Ok(eval(a, ctx)? - eval(b, ctx)?),
        Expr::Mul(a, b) => Ok(eval(a, ctx)? * eval(b, ctx)?),
        Expr::Div(a, b) => {
            let av = eval(a, ctx)?;
            let bv = eval(b, ctx)?;
            if bv == 0.0 {
                if av == 0.0 {
                    Ok(f64::NAN)
                } else {
                    Ok(if av.signum() > 0.0 { f64::INFINITY } else { f64::NEG_INFINITY })
                }
            } else {
                Ok(av / bv)
            }
        }
        Expr::Pow(a, b) => {
            let av = eval(a, ctx)?;
            let bv = eval(b, ctx)?;
            Ok(av.powf(bv))
        }
        Expr::Func(name, args) => {
            let values: Result<Vec<f64>> = args.iter().map(|a| eval(a, ctx)).collect();
            let values = values?;
            let f = ctx
                .funcs
                .get(name)
                .ok_or_else(|| MathError::UnknownFunction(name.clone()))?;
            match f {
                Func::Builtin(g) => g(&values),
                Func::User(body, params) => {
                    if params.len() != values.len() {
                        return Err(MathError::Eval(format!(
                            "function {} expects {} args, got {}",
                            name,
                            params.len(),
                            values.len()
                        )));
                    }
                    // Build a child context that shares funcs but has its own vars
                    let mut inner_vars = ctx.vars.clone();
                    for (p, v) in params.iter().zip(values.iter()) {
                        inner_vars.insert(p.clone(), *v);
                    }
                    let inner = Context { vars: inner_vars, funcs: ctx.funcs.clone() };
                    eval(body, &inner)
                }
            }
        }
    }
}

/// Evaluate an expression string in the standard context, optionally
/// overriding variables through `vars` (e.g., `&[("x", "2.0")]`).
pub fn eval_str(src: &str, vars: &[(&str, f64)]) -> Result<f64> {
    let e = crate::parser::Parser::parse(src)?;
    let mut ctx = Context::standard();
    for (k, v) in vars {
        ctx.set(*k, *v);
    }
    eval(&e, &ctx)
}

/// The standard library of mathematical functions exposed by the parser.
fn builtins() -> Vec<(&'static str, fn(&[f64]) -> Result<f64>)> {
    vec![
        ("sin", |a| unary(a, |x| Ok(x.sin()))),
        ("cos", |a| unary(a, |x| Ok(x.cos()))),
        ("tan", |a| unary(a, |x| Ok(x.tan()))),
        ("asin", |a| unary(a, |x| Ok(x.asin()))),
        ("acos", |a| unary(a, |x| Ok(x.acos()))),
        ("atan", |a| unary(a, |x| Ok(x.atan()))),
        ("sinh", |a| unary(a, |x| Ok(x.sinh()))),
        ("cosh", |a| unary(a, |x| Ok(x.cosh()))),
        ("tanh", |a| unary(a, |x| Ok(x.tanh()))),
        ("exp", |a| unary(a, |x| Ok(x.exp()))),
        ("ln", |a| unary(a, |x| domain(x > 0.0, "ln", x).map(|x| x.ln()))),
        ("log", |a| match a {
            [x, b] if *x > 0.0 && *b > 0.0 && *b != 1.0 => Ok(x.log(*b)),
            _ => Err(MathError::Domain(format!("log({}, {})", a.get(0).copied().unwrap_or(0.0), a.get(1).copied().unwrap_or(0.0)))),
        }),
        ("log2", |a| unary(a, |x| domain(x > 0.0, "log2", x).map(|x| x.log2()))),
        ("log10", |a| unary(a, |x| domain(x > 0.0, "log10", x).map(|x| x.log10()))),
        ("sqrt", |a| unary(a, |x| domain(x >= 0.0, "sqrt", x).map(|x| x.sqrt()))),
        ("cbrt", |a| unary(a, |x| Ok(x.cbrt()))),
        ("abs", |a| unary(a, |x| Ok(x.abs()))),
        ("floor", |a| unary(a, |x| Ok(x.floor()))),
        ("ceil", |a| unary(a, |x| Ok(x.ceil()))),
        ("round", |a| unary(a, |x| Ok(x.round()))),
        ("sign", |a| unary(a, |x| Ok(x.signum()))),
        ("min", |a| {
            if a.is_empty() {
                return Err(MathError::Eval("min needs at least one arg".into()));
            }
            Ok(a.iter().cloned().fold(f64::INFINITY, f64::min))
        }),
        ("max", |a| {
            if a.is_empty() {
                return Err(MathError::Eval("max needs at least one arg".into()));
            }
            Ok(a.iter().cloned().fold(f64::NEG_INFINITY, f64::max))
        }),
        ("pow", |a| match a {
            [x, y] => Ok(x.powf(*y)),
            _ => Err(MathError::Eval("pow(x, y) takes two args".into())),
        }),
        ("mod", |a| match a {
            [x, y] if *y != 0.0 => Ok(x.rem_euclid(*y)),
            _ => Err(MathError::Eval("mod requires non-zero divisor".into())),
        }),
        ("fract", |a| unary(a, |x| Ok(x.fract()))),
        ("gamma", |a| unary(a, |x| Ok(crate::special::gamma(x)))),
        ("erf", |a| unary(a, |x| Ok(crate::special::erf(x)))),
        ("erfc", |a| unary(a, |x| Ok(crate::special::erfc(x)))),
        ("sinc", |a| unary(a, |x| Ok(crate::special::sinc(x)))),
    ]
}

fn unary<F: FnOnce(f64) -> Result<f64>>(a: &[f64], f: F) -> Result<f64> {
    match a {
        [x] => f(*x),
        _ => Err(MathError::Eval("function expects exactly one argument".into())),
    }
}

fn domain(ok: bool, name: &str, x: f64) -> Result<f64> {
    if ok {
        Ok(x)
    } else {
        Err(MathError::Domain(format!("{} domain error at {}", name, x)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn ev(s: &str) -> f64 {
        eval_str(s, &[]).unwrap()
    }

    #[test]
    fn basic_eval() {
        assert_eq!(ev("1 + 2*3"), 7.0);
        assert_eq!(ev("2^10"), 1024.0);
        assert_eq!(ev("pi"), std::f64::consts::PI);
        assert_eq!(ev("sin(0)"), 0.0);
        assert_eq!(ev("cos(0)"), 1.0);
        assert_eq!(ev("log(8, 2)"), 3.0);
        assert_eq!(ev("sqrt(16)"), 4.0);
    }

    #[test]
    fn variable_eval() {
        let e = Parser::parse("x^2 + y").unwrap();
        let mut ctx = Context::standard();
        ctx.set("x", 3.0);
        ctx.set("y", 4.0);
        assert_eq!(eval(&e, &ctx).unwrap(), 13.0);
    }

    #[test]
    fn user_function() {
        let body = Parser::parse("a^2 + b^2").unwrap();
        let mut ctx = Context::standard();
        ctx.define("hypot", body, vec!["a".into(), "b".into()]);
        assert_eq!(eval(&Parser::parse("hypot(3, 4)").unwrap(), &ctx).unwrap(), 25.0);
    }

    #[test]
    fn domain_errors() {
        assert!(eval_str("sqrt(-1)", &[]).is_err());
        assert!(eval_str("ln(0)", &[]).is_err());
    }
}