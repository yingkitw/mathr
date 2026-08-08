//! Automatic differentiation using dual numbers.
//!
//! A dual number `a + b·ε` where `ε² = 0` carries both the value and
//! the derivative of a function at a point. Evaluating `f(Dual(x, 1))`
//! yields `Dual(f(x), f'(x))` in a single forward pass — no symbolic
//! differentiation needed, and it works for arbitrary compositions.
//!
//! This module provides:
//! - [`Dual`] — a 1st-order dual number type with full arithmetic
//! - [`eval`] — evaluate an [`Expr`] AST with dual numbers
//! - [`derivative`] — compute `f'(x)` at a point for any parseable expression
//! - [`gradient`] — compute the gradient of a multivariate expression
//! - [`jacobian`] — compute the Jacobian of a system of expressions

use crate::error::{MathError, Result};
use crate::expr::Expr;
use crate::eval::Context;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A first-order dual number: `value + deriv * ε` where `ε² = 0`.
///
/// Arithmetic on dual numbers automatically computes derivatives via
/// the chain rule. For `f: R → R`, evaluating `f(Dual(x, 1.0))` gives
/// `Dual(f(x), f'(x))`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// The primal value (the function value at the point).
    pub val: f64,
    /// The tangent value (the derivative at the point).
    pub deriv: f64,
}

impl Dual {
    /// Create a dual number with the given value and zero derivative (a constant).
    pub fn constant(val: f64) -> Self {
        Dual { val, deriv: 0.0 }
    }

    /// Create a dual number representing a variable at `val` with derivative 1.
    pub fn var(val: f64) -> Self {
        Dual { val, deriv: 1.0 }
    }

    /// Create a dual number with explicit value and derivative.
    pub fn new(val: f64, deriv: f64) -> Self {
        Dual { val, deriv }
    }

    // --- Elementary functions ---

    /// sin(a + bε) = sin(a) + b·cos(a)·ε
    pub fn sin(self) -> Self {
        Dual {
            val: self.val.sin(),
            deriv: self.deriv * self.val.cos(),
        }
    }

    /// cos(a + bε) = cos(a) - b·sin(a)·ε
    pub fn cos(self) -> Self {
        Dual {
            val: self.val.cos(),
            deriv: -self.deriv * self.val.sin(),
        }
    }

    /// tan(a + bε) = tan(a) + b·sec²(a)·ε
    pub fn tan(self) -> Self {
        let t = self.val.tan();
        Dual {
            val: t,
            deriv: self.deriv * (1.0 + t * t),
        }
    }

    /// exp(a + bε) = exp(a) + b·exp(a)·ε
    pub fn exp(self) -> Self {
        let e = self.val.exp();
        Dual {
            val: e,
            deriv: self.deriv * e,
        }
    }

    /// ln(a + bε) = ln(a) + b/a·ε
    pub fn ln(self) -> Self {
        Dual {
            val: self.val.ln(),
            deriv: self.deriv / self.val,
        }
    }

    /// log(a + bε) = ln(a)/ln(base) + b/(a·ln(base))·ε
    pub fn log(self, base: f64) -> Self {
        let ln_base = base.ln();
        Dual {
            val: self.val.ln() / ln_base,
            deriv: self.deriv / (self.val * ln_base),
        }
    }

    /// sqrt(a + bε) = sqrt(a) + b/(2·sqrt(a))·ε
    pub fn sqrt(self) -> Self {
        let s = self.val.sqrt();
        Dual {
            val: s,
            deriv: self.deriv / (2.0 * s),
        }
    }

    /// (a + bε)^c = a^c + c·a^(c-1)·b·ε  (constant exponent c)
    pub fn powf(self, c: f64) -> Self {
        let v = self.val.powf(c);
        let d = if self.val != 0.0 {
            c * self.val.powf(c - 1.0) * self.deriv
        } else if c > 1.0 {
            0.0 // derivative of 0^c is 0 for c > 1
        } else {
            f64::INFINITY
        };
        Dual { val: v, deriv: d }
    }

    /// (a + bε)^(c + dε) — general dual exponent.
    /// Uses: a^c = exp(c·ln(a)), so d/da(a^c) = c·a^(c-1) and d/dc(a^c) = a^c·ln(a).
    pub fn pow_dual(self, other: Dual) -> Self {
        if other.deriv == 0.0 {
            // Constant exponent — use the simpler formula
            self.powf(other.val)
        } else {
            // a^c = exp(c * ln(a))
            // d(a^c) = a^c * (c' * ln(a) + c * a'/a)
            let ln_a = self.val.ln();
            let val = self.val.powf(other.val);
            let deriv = val * (other.deriv * ln_a + other.val * self.deriv / self.val);
            Dual { val, deriv }
        }
    }

    /// asin(a + bε) = asin(a) + b/sqrt(1-a²)·ε
    pub fn asin(self) -> Self {
        Dual {
            val: self.val.asin(),
            deriv: self.deriv / (1.0 - self.val * self.val).sqrt(),
        }
    }

    /// acos(a + bε) = acos(a) - b/sqrt(1-a²)·ε
    pub fn acos(self) -> Self {
        Dual {
            val: self.val.acos(),
            deriv: -self.deriv / (1.0 - self.val * self.val).sqrt(),
        }
    }

    /// atan(a + bε) = atan(a) + b/(1+a²)·ε
    pub fn atan(self) -> Self {
        Dual {
            val: self.val.atan(),
            deriv: self.deriv / (1.0 + self.val * self.val),
        }
    }

    /// sinh(a + bε) = sinh(a) + b·cosh(a)·ε
    pub fn sinh(self) -> Self {
        Dual {
            val: self.val.sinh(),
            deriv: self.deriv * self.val.cosh(),
        }
    }

    /// cosh(a + bε) = cosh(a) + b·sinh(a)·ε
    pub fn cosh(self) -> Self {
        Dual {
            val: self.val.cosh(),
            deriv: self.deriv * self.val.sinh(),
        }
    }

    /// tanh(a + bε) = tanh(a) + b·(1-tanh²(a))·ε
    pub fn tanh(self) -> Self {
        let t = self.val.tanh();
        Dual {
            val: t,
            deriv: self.deriv * (1.0 - t * t),
        }
    }

    /// abs(a + bε) = |a| + b·sign(a)·ε
    pub fn abs(self) -> Self {
        Dual {
            val: self.val.abs(),
            deriv: self.deriv * self.val.signum(),
        }
    }
}

// --- Arithmetic operator overloads ---

impl Add for Dual {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        Dual {
            val: self.val + other.val,
            deriv: self.deriv + other.deriv,
        }
    }
}

impl Sub for Dual {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        Dual {
            val: self.val - other.val,
            deriv: self.deriv - other.deriv,
        }
    }
}

impl Mul for Dual {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        // (a + bε)(c + dε) = ac + (ad + bc)ε
        Dual {
            val: self.val * other.val,
            deriv: self.val * other.deriv + self.deriv * other.val,
        }
    }
}

impl Div for Dual {
    type Output = Dual;
    fn div(self, other: Dual) -> Dual {
        // (a + bε)/(c + dε) = a/c + (bc - ad)/c² · ε
        let val = self.val / other.val;
        let deriv = (self.deriv * other.val - self.val * other.deriv) / (other.val * other.val);
        Dual { val, deriv }
    }
}

impl Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Dual {
        Dual {
            val: -self.val,
            deriv: -self.deriv,
        }
    }
}

// --- Expr evaluation with dual numbers ---

/// Evaluate an [`Expr`] AST using dual numbers, returning the value and derivative.
///
/// The `var` argument is the variable to differentiate with respect to.
/// All other variables are looked up in `ctx` as constants.
pub fn eval(expr: &Expr, var: &str, x: f64, ctx: &Context) -> Result<Dual> {
    match expr {
        Expr::Num(n) => Ok(Dual::constant(*n)),
        Expr::Var(name) => {
            if name == var {
                Ok(Dual::var(x))
            } else {
                let v = ctx
                    .vars
                    .get(name)
                    .ok_or_else(|| MathError::UnknownVariable(name.clone()))?;
                Ok(Dual::constant(*v))
            }
        }
        Expr::Neg(e) => Ok(-eval(e, var, x, ctx)?),
        Expr::Add(a, b) => Ok(eval(a, var, x, ctx)? + eval(b, var, x, ctx)?),
        Expr::Sub(a, b) => Ok(eval(a, var, x, ctx)? - eval(b, var, x, ctx)?),
        Expr::Mul(a, b) => Ok(eval(a, var, x, ctx)? * eval(b, var, x, ctx)?),
        Expr::Div(a, b) => Ok(eval(a, var, x, ctx)? / eval(b, var, x, ctx)?),
        Expr::Pow(a, b) => {
            let da = eval(a, var, x, ctx)?;
            let db = eval(b, var, x, ctx)?;
            Ok(da.pow_dual(db))
        }
        Expr::Func(name, args) => {
            if args.len() != 1 {
                // Multi-arg function: finite-difference derivative
                let h = 1e-8;
                let mut cx0 = ctx.clone();
                cx0.set(var, x);
                let v0 = crate::eval::eval(expr, &cx0)?;
                let mut cx1 = ctx.clone();
                cx1.set(var, x + h);
                let v1 = crate::eval::eval(expr, &cx1)?;
                return Ok(Dual::new(v0, (v1 - v0) / h));
            }
            let d = eval(&args[0], var, x, ctx)?;
            match name.as_str() {
                "sin" => Ok(d.sin()),
                "cos" => Ok(d.cos()),
                "tan" => Ok(d.tan()),
                "exp" => Ok(d.exp()),
                "ln" | "log" if name == "ln" => Ok(d.ln()),
                "log" => Ok(d.log(10.0)),
                "log2" => Ok(d.log(2.0)),
                "sqrt" => Ok(d.sqrt()),
                "asin" => Ok(d.asin()),
                "acos" => Ok(d.acos()),
                "atan" => Ok(d.atan()),
                "sinh" => Ok(d.sinh()),
                "cosh" => Ok(d.cosh()),
                "tanh" => Ok(d.tanh()),
                "abs" => Ok(d.abs()),
                // Unknown function: evaluate numerically, derivative via finite difference
                _ => {
                    let h = 1e-8;
                    let v0 = d.val;
                    let v1 = {
                        let mut cx = ctx.clone();
                        cx.set(var, x + h);
                        crate::eval::eval(expr, &cx)?
                    };
                    let deriv = (v1 - v0) / h;
                    Ok(Dual::new(v0, deriv))
                }
            }
        }
    }
}

/// Compute the derivative of `f(var)` at `x` using automatic differentiation.
///
/// # Example
/// ```ignore
/// use mathr::autodiff::derivative;
/// use mathr::parser::Parser;
/// use mathr::eval::Context;
///
/// let expr = Parser::parse("sin(x) * x^2").unwrap();
/// let ctx = Context::standard();
/// let d = derivative(&expr, "x", 1.0, &ctx).unwrap();
/// // f(1) = sin(1) * 1 = 0.8415...
/// // f'(1) = cos(1)*1 + sin(1)*2 = 2.4434...
/// ```
pub fn derivative(expr: &Expr, var: &str, x: f64, ctx: &Context) -> Result<Dual> {
    eval(expr, var, x, ctx)
}

/// Compute the gradient of a multivariate expression at a point.
///
/// Returns a vector of `(variable_name, derivative)` pairs, one for each
/// variable found in the expression.
pub fn gradient(expr: &Expr, point: &Context) -> Result<Vec<(String, f64)>> {
    // Collect all variables in the expression
    let vars = collect_vars(expr);
    let mut result = Vec::with_capacity(vars.len());
    for var in &vars {
        let x = point
            .vars
            .get(var)
            .ok_or_else(|| MathError::UnknownVariable(var.clone()))?;
        let d = eval(expr, var, *x, point)?;
        result.push((var.clone(), d.deriv));
    }
    Ok(result)
}

/// Compute the Jacobian matrix of a system of expressions.
///
/// Returns a matrix where `result[i][j]` is `∂f_i/∂x_j` evaluated at the point.
/// Variables are ordered alphabetically.
pub fn jacobian(exprs: &[Expr], point: &Context) -> Result<Vec<Vec<f64>>> {
    let mut all_vars: std::collections::HashSet<String> = std::collections::HashSet::new();
    for expr in exprs {
        for v in collect_vars(expr) {
            all_vars.insert(v);
        }
    }
    let mut all_vars: Vec<String> = all_vars.into_iter().collect();
    all_vars.sort();

    let mut jacobian = Vec::with_capacity(exprs.len());
    for expr in exprs {
        let mut row = Vec::with_capacity(all_vars.len());
        for var in &all_vars {
            let x = point
                .vars
                .get(var)
                .ok_or_else(|| MathError::UnknownVariable(var.clone()))?;
            let d = eval(expr, var, *x, point)?;
            row.push(d.deriv);
        }
        jacobian.push(row);
    }
    Ok(jacobian)
}

/// Collect all variable names appearing in an expression, sorted alphabetically.
fn collect_vars(expr: &Expr) -> Vec<String> {
    let mut vars = std::collections::HashSet::new();
    collect_vars_inner(expr, &mut vars);
    let mut v: Vec<_> = vars.into_iter().collect();
    v.sort();
    v
}

fn collect_vars_inner(expr: &Expr, vars: &mut std::collections::HashSet<String>) {
    match expr {
        Expr::Var(name) => {
            vars.insert(name.clone());
        }
        Expr::Num(_) => {}
        Expr::Neg(e) => collect_vars_inner(e, vars),
        Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) | Expr::Pow(a, b) => {
            collect_vars_inner(a, vars);
            collect_vars_inner(b, vars);
        }
        Expr::Func(_, args) => {
            for a in args {
                collect_vars_inner(a, vars);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn dual_arithmetic() {
        let a = Dual::new(3.0, 1.0);
        let b = Dual::new(2.0, 0.0);
        let c = a + b;
        assert!(close(c.val, 5.0, 1e-10) && close(c.deriv, 1.0, 1e-10));
        let d = a * b;
        assert!(close(d.val, 6.0, 1e-10) && close(d.deriv, 2.0, 1e-10));
        let e = a / b;
        assert!(close(e.val, 1.5, 1e-10) && close(e.deriv, 0.5, 1e-10));
    }

    #[test]
    fn dual_sin() {
        let d = Dual::var(0.0).sin();
        assert!(close(d.val, 0.0, 1e-10));
        assert!(close(d.deriv, 1.0, 1e-10)); // cos(0) = 1
    }

    #[test]
    fn dual_exp() {
        let d = Dual::var(0.0).exp();
        assert!(close(d.val, 1.0, 1e-10));
        assert!(close(d.deriv, 1.0, 1e-10)); // exp(0) = 1
    }

    #[test]
    fn dual_powf() {
        let d = Dual::var(3.0).powf(2.0);
        assert!(close(d.val, 9.0, 1e-10));
        assert!(close(d.deriv, 6.0, 1e-10)); // 2*3 = 6
    }

    #[test]
    fn dual_log() {
        let d = Dual::var(1.0).ln();
        assert!(close(d.val, 0.0, 1e-10));
        assert!(close(d.deriv, 1.0, 1e-10)); // 1/1 = 1
    }

    #[test]
    fn dual_sqrt() {
        let d = Dual::var(4.0).sqrt();
        assert!(close(d.val, 2.0, 1e-10));
        assert!(close(d.deriv, 0.25, 1e-10)); // 1/(2*2) = 0.25
    }

    #[test]
    fn dual_tan() {
        let d = Dual::var(0.0).tan();
        assert!(close(d.val, 0.0, 1e-10));
        assert!(close(d.deriv, 1.0, 1e-10)); // sec²(0) = 1
    }

    #[test]
    fn dual_atan() {
        let d = Dual::var(0.0).atan();
        assert!(close(d.val, 0.0, 1e-10));
        assert!(close(d.deriv, 1.0, 1e-10)); // 1/(1+0) = 1
    }

    #[test]
    fn dual_sinh_cosh_tanh() {
        let s = Dual::var(0.0).sinh();
        assert!(close(s.val, 0.0, 1e-10) && close(s.deriv, 1.0, 1e-10));
        let c = Dual::var(0.0).cosh();
        assert!(close(c.val, 1.0, 1e-10) && close(c.deriv, 0.0, 1e-10));
        let t = Dual::var(0.0).tanh();
        assert!(close(t.val, 0.0, 1e-10) && close(t.deriv, 1.0, 1e-10));
    }

    #[test]
    fn autodiff_polynomial() {
        // f(x) = x^3 + 2x^2 - x + 5, f'(x) = 3x^2 + 4x - 1
        let expr = Parser::parse("x^3 + 2*x^2 - x + 5").unwrap();
        let ctx = Context::standard();
        for x in [-2.0, -1.0, 0.0, 1.0, 3.5] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = x * x * x + 2.0 * x * x - x + 5.0;
            let expected_deriv = 3.0 * x * x + 4.0 * x - 1.0;
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_trig_composition() {
        // f(x) = sin(x^2), f'(x) = 2x*cos(x^2)
        let expr = Parser::parse("sin(x^2)").unwrap();
        let ctx = Context::standard();
        for x in [0.0, 0.5, 1.0, 2.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = (x * x).sin();
            let expected_deriv = 2.0 * x * (x * x).cos();
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_exp_composition() {
        // f(x) = exp(-(x^2)), f'(x) = -2x*exp(-x^2)
        // Note: parser treats -x^2 as (-x)^2, so we use explicit parens.
        let expr = Parser::parse("exp(-(x^2))").unwrap();
        let ctx = Context::standard();
        for x in [0.0, 0.5, 1.0, 2.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = (-x * x).exp();
            let expected_deriv = -2.0 * x * (-x * x).exp();
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_log_sqrt() {
        // f(x) = ln(sqrt(x)), f'(x) = 1/(2x)
        let expr = Parser::parse("ln(sqrt(x))").unwrap();
        let ctx = Context::standard();
        for x in [1.0, 4.0, 9.0, 100.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = x.sqrt().ln();
            let expected_deriv = 1.0 / (2.0 * x);
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_product_rule() {
        // f(x) = x * sin(x), f'(x) = sin(x) + x*cos(x)
        let expr = Parser::parse("x * sin(x)").unwrap();
        let ctx = Context::standard();
        for x in [0.0, 1.0, 2.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = x * x.sin();
            let expected_deriv = x.sin() + x * x.cos();
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_quotient_rule() {
        // f(x) = sin(x) / x, f'(x) = (x*cos(x) - sin(x)) / x^2
        let expr = Parser::parse("sin(x) / x").unwrap();
        let ctx = Context::standard();
        for x in [1.0, 2.0, 5.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = x.sin() / x;
            let expected_deriv = (x * x.cos() - x.sin()) / (x * x);
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_chain_rule_deep() {
        // f(x) = cos(sin(x^2 + 1))
        // f'(x) = -sin(sin(x^2+1)) * cos(x^2+1) * 2x
        let expr = Parser::parse("cos(sin(x^2 + 1))").unwrap();
        let ctx = Context::standard();
        for x in [0.5, 1.0, 1.5] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let inner = x * x + 1.0;
            let expected_val = inner.sin().cos();
            let expected_deriv = -inner.sin().sin() * inner.cos() * 2.0 * x;
            assert!(close(d.val, expected_val, 1e-10), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-10), "deriv at x={}", x);
        }
    }

    #[test]
    fn autodiff_gradient() {
        // f(x, y) = x^2 + y^3, ∇f = (2x, 3y^2)
        let expr = Parser::parse("x^2 + y^3").unwrap();
        let mut ctx = Context::standard();
        ctx.set("x", 2.0);
        ctx.set("y", 3.0);
        let grad = gradient(&expr, &ctx).unwrap();
        // Should have entries for x and y
        let x_grad = grad.iter().find(|(name, _)| name == "x").unwrap().1;
        let y_grad = grad.iter().find(|(name, _)| name == "y").unwrap().1;
        assert!(close(x_grad, 4.0, 1e-10)); // 2*2
        assert!(close(y_grad, 27.0, 1e-10)); // 3*9
    }

    #[test]
    fn autodiff_jacobian() {
        // f1 = x^2 + y, f2 = x * y^2
        // J = [[2x, 1], [y^2, 2xy]]
        let f1 = Parser::parse("x^2 + y").unwrap();
        let f2 = Parser::parse("x * y^2").unwrap();
        let mut ctx = Context::standard();
        ctx.set("x", 2.0);
        ctx.set("y", 3.0);
        let jac = jacobian(&[f1, f2], &ctx).unwrap();
        // J[0] = [2*2, 1] = [4, 1]
        assert!(close(jac[0][0], 4.0, 1e-10));
        assert!(close(jac[0][1], 1.0, 1e-10));
        // J[1] = [9, 2*2*3] = [9, 12]
        assert!(close(jac[1][0], 9.0, 1e-10));
        assert!(close(jac[1][1], 12.0, 1e-10));
    }

    #[test]
    fn autodiff_constant() {
        let expr = Parser::parse("42").unwrap();
        let ctx = Context::standard();
        let d = derivative(&expr, "x", 1.0, &ctx).unwrap();
        assert!(close(d.val, 42.0, 1e-10));
        assert!(close(d.deriv, 0.0, 1e-10));
    }

    #[test]
    fn autodiff_var_not_in_expr() {
        let expr = Parser::parse("y^2").unwrap();
        let mut ctx = Context::standard();
        ctx.set("y", 5.0);
        let d = derivative(&expr, "x", 1.0, &ctx).unwrap();
        // y is treated as a constant, so derivative w.r.t. x is 0
        assert!(close(d.val, 25.0, 1e-10));
        assert!(close(d.deriv, 0.0, 1e-10));
    }

    #[test]
    fn autodiff_powf_general() {
        // f(x) = x^x, f'(x) = x^x * (ln(x) + 1)
        let expr = Parser::parse("x^x").unwrap();
        let ctx = Context::standard();
        for x in [1.0, 2.0, 3.0] {
            let d = derivative(&expr, "x", x, &ctx).unwrap();
            let expected_val = x.powf(x);
            let expected_deriv = x.powf(x) * (x.ln() + 1.0);
            assert!(close(d.val, expected_val, 1e-8), "val at x={}", x);
            assert!(close(d.deriv, expected_deriv, 1e-8), "deriv at x={}", x);
        }
    }
}
