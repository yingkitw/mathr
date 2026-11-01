//! Polynomial expansion (distribution) over the `Expr` AST.
//!
//! [`expand`] rewrites products and non-negative integer powers of
//! polynomial sub-expressions into a collected sum of monomials, e.g.
//! `(x+1)^3 → x^3 + 3·x^2 + 3·x + 1` and `(x+y)·(x−y) → x^2 − y^2`
//! (multivariate, coefficients are `f64`).
//!
//! Sub-expressions that are not polynomial (functions, symbolic powers,
//! variables in denominators) are left intact but their children are still
//! expanded: `sin(x+1)·(x+2) → (x+2)·sin(x+1)` keeps the call opaque while
//! `(x+2)` factors distribute wherever the other side is polynomial.
//!
//! # Example
//!
//! ```
//! use mathr::parser::Parser;
//! use mathr::poly::expand;
//!
//! let e = Parser::parse("(x+1)^3").unwrap();
//! assert!(expand(&e).to_string().contains("x^3"));
//! ```

use crate::expr::Expr::{self, *};
use std::collections::BTreeMap;

/// One monomial: `coeff · ∏ var^exp`, with variables in sorted order.
#[derive(Debug, Clone, PartialEq)]
pub struct Term {
    pub coeff: f64,
    pub monomial: BTreeMap<String, u32>,
}

/// A polynomial: a sum of terms. Like terms are combined and zero
/// coefficients dropped by the constructors.
pub type Poly = Vec<Term>;

/// One term of a univariate expansion in `var`: `coeff · var^exp`.
///
/// The coefficient is an arbitrary [`Expr`] — it may contain other
/// variables, function calls, or anything else. Only the exponent of the
/// *target* variable is tracked.
#[derive(Debug, Clone, PartialEq)]
pub struct UniTerm {
    /// The opaque coefficient (may mention other variables).
    pub coeff: Expr,
    /// The exponent of the target variable.
    pub exp: u32,
}

/// A univariate expansion: a sum of `coeff · var^exp` terms.
pub type UniPoly = Vec<UniTerm>;

/// Safety caps so hostile inputs cannot explode memory.
const MAX_TERMS: usize = 20_000;
const MAX_EXPONENT: u32 = 64;
const MAX_MONOMIAL_EXP: u32 = 10_000;

/// Expand products and non-negative integer powers into a collected sum of
/// monomials wherever the sub-expression is polynomial. Non-polynomial
/// parts (functions, symbolic powers, denominators with variables) are
/// kept, but their children are still expanded.
pub fn expand(expr: &Expr) -> Expr {
    let expanded: Expr = match expr {
        Num(_) | Var(_) => expr.clone(),
        Neg(a) => Neg(Box::new(expand(a))),
        Add(a, b) => Add(Box::new(expand(a)), Box::new(expand(b))),
        Sub(a, b) => Sub(Box::new(expand(a)), Box::new(expand(b))),
        Mul(a, b) => Mul(Box::new(expand(a)), Box::new(expand(b))),
        Div(a, b) => Div(Box::new(expand(a)), Box::new(expand(b))),
        Pow(b, e) => Pow(Box::new(expand(b)), Box::new(expand(e))),
        Func(name, args) => Func(name.clone(), args.iter().map(expand).collect()),
    };
    match to_poly(&expanded) {
        Some(p) => poly_to_expr(&p),
        None => expanded,
    }
}

/// Expand an expression **in one variable**, treating every other
/// sub-expression as an opaque exact coefficient.
///
/// The multivariate [`expand`] can only distribute when both sides are
/// polynomials, so `(x + sin(y))^2` stays unexpanded. This mode expands in
/// `var` only, and anything that is not a power of `var` becomes part of
/// the coefficient — exactly, as an expression, not as a float.
///
/// ```text
/// (x + sin(y))^2  →  x^2 + 2·x·sin(y) + sin(y)^2
/// (x + 2·y)^3    →  x^3 + 6·x^2·y + 12·x·y^2 + 8·y^3
/// ```
///
/// This is the "termsOf" style the CASs use: the target variable is
/// expanded, the rest is carried along.
///
/// Cross-port pattern from unsolve `poly.ts:132-210`.
pub fn expand_in(expr: &Expr, var: &str) -> Expr {
    match to_uni_poly(expr, var) {
        Some(p) => uni_poly_to_expr(&p, var),
        None => expand(expr),
    }
}

/// Interpret an expression as a univariate polynomial in `var`.
///
/// Returns `None` when the expression cannot be expanded in `var` — for
/// example when `var` appears in a denominator, or under a symbolic power.
pub fn to_uni_poly(expr: &Expr, var: &str) -> Option<UniPoly> {
    match expr {
        Num(x) => Some(vec![UniTerm {
            coeff: Expr::num(*x),
            exp: 0,
        }]),
        Var(name) => {
            if name == var {
                Some(vec![UniTerm {
                    coeff: Expr::num(1.0),
                    exp: 1,
                }])
            } else {
                // Another variable: an opaque coefficient.
                Some(vec![UniTerm {
                    coeff: expr.clone(),
                    exp: 0,
                }])
            }
        }
        Neg(a) => Some(
            to_uni_poly(a, var)?
                .into_iter()
                .map(|mut t| {
                    t.coeff = Expr::neg(t.coeff);
                    t
                })
                .collect(),
        ),        Add(a, b) => Some(uni_add(&to_uni_poly(a, var)?, &to_uni_poly(b, var)?)),
        Sub(a, b) => Some(uni_add(
            &to_uni_poly(a, var)?,
            &uni_neg(&to_uni_poly(b, var)?),
        )),
        Mul(a, b) => Some(uni_mul(&to_uni_poly(a, var)?, &to_uni_poly(b, var)?)?),
        Div(a, b) => {
            // Only a constant denominator keeps the expression polynomial.
            let pb = to_uni_poly(b, var)?;
            if pb.len() != 1 || pb[0].exp != 0 {
                return None;
            }
            let pa = to_uni_poly(a, var)?;
            Some(
                pa.into_iter()
                    .map(|mut t| {
                        t.coeff = Expr::div(t.coeff, pb[0].coeff.clone());
                        t
                    })
                    .collect(),
            )
        }
        Pow(base, e) => {
            // Only a non-negative integer exponent keeps it polynomial.
            let (Expr::Num(n),) = (e.as_ref(),) else {
                return None;
            };
            if n.fract() != 0.0 || *n < 0.0 || *n > MAX_EXPONENT as f64 {
                return None;
            }
            let k = *n as u32;
            let pb = to_uni_poly(base, var)?;
            Some(uni_powi(&pb, k)?)
        }
        Func(_, _) => {
            // A function call is an opaque coefficient — but only if it does
            // not contain the target variable in a way that would need
            // expanding. Treat the whole call as opaque.
            Some(vec![UniTerm {
                coeff: expr.clone(),
                exp: 0,
            }])
        }
    }
}

/// Negate every coefficient of a univariate polynomial.
fn uni_neg(p: &UniPoly) -> UniPoly {
    p.iter()
        .map(|t| UniTerm {
            coeff: Expr::neg(t.coeff.clone()),
            exp: t.exp,
        })
        .collect()
}

/// Add two univariate polynomials, combining like exponents.
fn uni_add(a: &UniPoly, b: &UniPoly) -> UniPoly {
    let mut out: UniPoly = a.clone();
    for t in b {
        match out.iter_mut().find(|u| u.exp == t.exp) {
            Some(u) => {
                u.coeff = Expr::add(u.coeff.clone(), t.coeff.clone());
            }
            None => out.push(t.clone()),
        }
    }
    // Drop zero coefficients.
    out.retain(|t| !is_zero_expr(&t.coeff));
    out
}

/// Multiply two univariate polynomials.
fn uni_mul(a: &UniPoly, b: &UniPoly) -> Option<UniPoly> {
    let mut out: UniPoly = Vec::new();
    for ta in a {
        for tb in b {
            let exp = ta.exp + tb.exp;
            if exp > MAX_EXPONENT {
                return None;
            }
            let coeff = Expr::mul(ta.coeff.clone(), tb.coeff.clone());
            match out.iter_mut().find(|u| u.exp == exp) {
                Some(u) => u.coeff = Expr::add(u.coeff.clone(), coeff),
                None => out.push(UniTerm { coeff, exp }),
            }
            if out.len() > MAX_TERMS {
                return None;
            }
        }
    }
    out.retain(|t| !is_zero_expr(&t.coeff));
    Some(out)
}

/// Raise a univariate polynomial to a non-negative integer power.
fn uni_powi(p: &UniPoly, k: u32) -> Option<UniPoly> {
    if k == 0 {
        return Some(vec![UniTerm {
            coeff: Expr::num(1.0),
            exp: 0,
        }]);
    }
    let mut result: UniPoly = vec![UniTerm {
        coeff: Expr::num(1.0),
        exp: 0,
    }];
    for _ in 0..k {
        result = uni_mul(&result, p)?;
        if result.len() > MAX_TERMS {
            return None;
        }
    }
    Some(result)
}

/// Is this expression the number zero?
fn is_zero_expr(e: &Expr) -> bool {
    matches!(e, Expr::Num(n) if *n == 0.0)
}

/// Render a univariate polynomial back to an expression.
///
/// The coefficients are simplified, so the result reads `x^2 + 2·x·sin(y) +
/// sin(y)^2` rather than the raw product form.
fn uni_poly_to_expr(p: &UniPoly, var: &str) -> Expr {
    if p.is_empty() {
        return Expr::num(0.0);
    }
    // Sort by descending exponent for a readable result.
    let mut sorted = p.clone();
    sorted.sort_by_key(|t| std::cmp::Reverse(t.exp));
    let mut out: Option<Expr> = None;
    for t in sorted {
        // The coefficient may be a product of factors (`2*y*y`); expanding it
        // folds the repeated factors into powers (`2*y^2`).
        let coeff = crate::simplify::simplify(&expand(&t.coeff));
        let term = if t.exp == 0 {
            coeff
        } else {
            let x = Expr::var(var.to_string());
            // An exponent of 1 is redundant: `2*x`, not `2*x^1`.
            let pow = if t.exp == 1 {
                x
            } else {
                Expr::pow(x, Expr::num(t.exp as f64))
            };
            if matches!(coeff, Expr::Num(n) if n == 1.0) {
                pow
            } else {
                Expr::mul(coeff, pow)
            }
        };
        out = Some(match out {
            None => term,
            Some(acc) => Expr::add(acc, term),
        });
    }
    out.unwrap_or_else(|| Expr::num(0.0))
}

pub fn to_poly(expr: &Expr) -> Option<Poly> {
    match expr {
        Num(x) => Some(vec![Term {
            coeff: *x,
            monomial: BTreeMap::new(),
        }]),
        Var(name) => Some(vec![Term {
            coeff: 1.0,
            monomial: [(name.clone(), 1u32)].into_iter().collect(),
        }]),
        Neg(a) => to_poly(a).map(|p| p.into_iter().map(|mut t| { t.coeff = -t.coeff; t }).collect()),
        Add(a, b) => Some(poly_add(&to_poly(a)?, &to_poly(b)?)),
        Sub(a, b) => Some(poly_add(&to_poly(a)?, &poly_neg(&to_poly(b)?))),
        Mul(a, b) => Some(poly_mul(&to_poly(a)?, &to_poly(b)?)?),
        Div(a, b) => {
            let pa = to_poly(a)?;
            let pb = to_poly(b)?;
            // Only constant denominators stay polynomial.
            if pb.len() != 1 || !pb[0].monomial.is_empty() || pb[0].coeff == 0.0 {
                return None;
            }
            let inv = 1.0 / pb[0].coeff;
            Some(pa.into_iter().map(|mut t| { t.coeff *= inv; t }).collect())
        }
        Pow(base, exp) => {
            let Expr::Num(k) = &**exp else { return None };
            if *k < 0.0 || k.fract() != 0.0 || *k > MAX_EXPONENT as f64 {
                return None;
            }
            let pb = to_poly(base)?;
            poly_powi(&pb, *k as u32)
        }
        Func(_, _) => None,
    }
}

/// Convert a polynomial back into an `Expr` sum, ordered by descending
/// total degree (then lexicographic exponent order, so `x^2` precedes
/// `x*y`). Negative terms are joined with `Sub` for readable output.
pub fn poly_to_expr(poly: &Poly) -> Expr {
    if poly.is_empty() {
        return Num(0.0);
    }
    let mut terms: Vec<&Term> = poly.iter().collect();
    terms.sort_by(|a, b| {
        let da: u32 = a.monomial.values().sum();
        let db: u32 = b.monomial.values().sum();
        db.cmp(&da).then_with(|| lex_cmp(&a.monomial, &b.monomial))
    });

    let mut out: Option<Expr> = None;
    for t in terms {
        let mag = Term {
            coeff: t.coeff.abs(),
            monomial: t.monomial.clone(),
        };
        out = Some(match out {
            None => term_expr(t),
            Some(acc) => {
                if t.coeff < 0.0 {
                    Sub(Box::new(acc), Box::new(term_expr(&mag)))
                } else {
                    Add(Box::new(acc), Box::new(term_expr(&mag)))
                }
            }
        });
    }
    out.unwrap()
}

/// Build the expression for one term. A negative coefficient renders
/// through `Num(-c)` (`-2*x`) or `Neg` (`-x`).
fn term_expr(t: &Term) -> Expr {
    if t.monomial.is_empty() {
        return Num(t.coeff);
    }
    let mono = monomial_to_expr(&t.monomial);
    if t.coeff == 1.0 {
        return mono;
    }
    if t.coeff == -1.0 {
        return Neg(Box::new(mono));
    }
    Mul(Box::new(Num(t.coeff)), Box::new(mono))
}

/// Lexicographic comparison of monomials over the union of their variables
/// (exponents compared in descending order, so `x^2` sorts before `x*y`).
fn lex_cmp(a: &BTreeMap<String, u32>, b: &BTreeMap<String, u32>) -> std::cmp::Ordering {
    use std::cmp::Ordering;
    let mut names: Vec<&String> = a.keys().chain(b.keys()).collect();
    names.sort();
    names.dedup();
    for n in names {
        let ea = a.get(n).copied().unwrap_or(0);
        let eb = b.get(n).copied().unwrap_or(0);
        if ea != eb {
            return eb.cmp(&ea);
        }
    }
    Ordering::Equal
}

fn monomial_to_expr(m: &BTreeMap<String, u32>) -> Expr {
    let factors: Vec<Expr> = m
        .iter()
        .map(|(name, exp)| {
            if *exp == 1 {
                Var(name.clone())
            } else {
                Pow(Box::new(Var(name.clone())), Box::new(Num(*exp as f64)))
            }
        })
        .collect();
    match factors.len() {
        0 => Num(1.0),
        1 => factors.into_iter().next().unwrap(),
        _ => factors
            .into_iter()
            .reduce(|acc, f| Mul(Box::new(acc), Box::new(f)))
            .unwrap(),
    }
}

/// Combine like terms and drop exact zeros.
fn normalize(poly: Poly) -> Poly {
    let mut combined: BTreeMap<BTreeMap<String, u32>, f64> = BTreeMap::new();
    for t in poly {
        *combined.entry(t.monomial).or_insert(0.0) += t.coeff;
    }
    combined
        .into_iter()
        .filter(|(_, c)| *c != 0.0)
        .map(|(monomial, coeff)| Term { coeff, monomial })
        .collect()
}

fn poly_neg(p: &Poly) -> Poly {
    p.iter()
        .map(|t| Term {
            coeff: -t.coeff,
            monomial: t.monomial.clone(),
        })
        .collect()
}

fn poly_add(a: &Poly, b: &Poly) -> Poly {
    let mut out = a.clone();
    out.extend(b.clone());
    normalize(out)
}

fn poly_mul(a: &Poly, b: &Poly) -> Option<Poly> {
    if a.len().saturating_mul(b.len()) > MAX_TERMS {
        return None;
    }
    let mut out = Vec::with_capacity(a.len() * b.len());
    for ta in a {
        for tb in b {
            let mut monomial = ta.monomial.clone();
            for (name, e) in &tb.monomial {
                let sum = monomial.entry(name.clone()).or_insert(0);
                *sum = sum.saturating_add(*e);
                if *sum > MAX_MONOMIAL_EXP {
                    return None;
                }
            }
            out.push(Term {
                coeff: ta.coeff * tb.coeff,
                monomial,
            });
        }
    }
    Some(normalize(out))
}

fn poly_powi(p: &Poly, k: u32) -> Option<Poly> {
    let mut result = vec![Term {
        coeff: 1.0,
        monomial: BTreeMap::new(),
    }];
    let mut base = p.clone();
    let mut kk = k;
    while kk > 0 {
        if kk & 1 == 1 {
            result = poly_mul(&result, &base)?;
        }
        kk >>= 1;
        if kk > 0 {
            base = poly_mul(&base, &base)?;
        }
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn expand_str(src: &str) -> String {
        expand(&Parser::parse(src).unwrap()).to_string()
    }

    fn expands_to(src: &str, expected: &str) {
        let e = expand(&Parser::parse(src).unwrap());
        let target = Parser::parse(expected).unwrap();
        assert!(
            e.equals(&target),
            "expand({}) = {} but expected {}",
            src,
            e,
            target
        );
    }

    // ---- univariate expansion with opaque coefficients ----------------------

    fn expand_in_str(src: &str, var: &str) -> String {
        expand_in(&Parser::parse(src).unwrap(), var).to_string()
    }

    /// The expansion must agree with the original at sample points.
    fn expand_in_agrees(src: &str, var: &str, points: &[(f64, f64)]) {
        let e = Parser::parse(src).unwrap();
        let out = expand_in(&e, var);
        for &(x, y) in points {
            let mut c1 = crate::eval::Context::standard();
            c1.set(var, x);
            c1.set("y", y);
            let c2 = c1.clone();
            let a = crate::eval::eval(&e, &c1).unwrap();
            let b = crate::eval::eval(&out, &c2).unwrap();
            assert!(
                (a - b).abs() < 1e-9,
                "{} in {}: at x={} y={} got {} want {}",
                src,
                var,
                x,
                y,
                b,
                a
            );
        }
    }

    #[test]
    fn expand_in_keeps_function_calls_opaque() {
        // The point of the mode: a function call is an exact coefficient.
        let out = expand_in_str("(x + sin(y))^2", "x");
        assert!(out.contains("sin(y)"), "got {}", out);
        assert!(out.contains("x^2"), "got {}", out);
    }

    #[test]
    fn expand_in_produces_correct_exponents() {
        // (x + 2y)³ = x³ + 6yx² + 12y²x + 8y³
        let out = expand_in_str("(x + 2*y)^3", "x");
        assert!(out.contains("x^3"), "got {}", out);
        assert!(out.contains("8*y^3"), "got {}", out);
        assert!(out.contains("12*y^2*x"), "got {}", out);
    }

    #[test]
    fn expand_in_agrees_with_the_original() {
        let pts = [(0.7, 0.5), (1.3, 1.7), (-2.1, 0.3)];
        expand_in_agrees("(x + sin(y))^2", "x", &pts);
        expand_in_agrees("(x + 2*y)^3", "x", &pts);
        expand_in_agrees("(x+1)^2", "x", &pts);
        expand_in_agrees("(x-y)^3", "x", &pts);
    }

    #[test]
    fn expand_in_is_univariate() {
        // Expanding in y leaves x inside the coefficient.
        let out = expand_in_str("(x + sin(y))^2", "y");
        assert!(out.contains("sin(y)"), "got {}", out);
        // The target variable is y, so the result is a power of y.
        let e = Parser::parse(&out).unwrap();
        // Every term's exponent of y is tracked; x only appears in a coefficient.
        let uni = to_uni_poly(&e, "y").unwrap();
        assert!(uni.iter().all(|t| t.exp <= 2), "got {:?}", uni);
    }

    #[test]
    fn expand_in_handles_a_plain_polynomial() {
        // The ordinary case still works.
        let out = expand_in_str("(x+1)^2", "x");
        assert_eq!(out, "x^2 + 2*x + 1");
    }

    #[test]
    fn expand_in_rejects_a_symbolic_power() {
        // x^y is not a polynomial in x.
        let e = Parser::parse("x^y").unwrap();
        assert!(to_uni_poly(&e, "x").is_none());
    }

    #[test]
    fn expand_in_rejects_a_variable_in_a_denominator() {
        let e = Parser::parse("1/x").unwrap();
        assert!(to_uni_poly(&e, "x").is_none());
    }

    #[test]
    fn expand_cube() {
        expands_to("(x+1)^3", "x^3 + 3*x^2 + 3*x + 1");
    }

    #[test]
    fn expand_square_of_sum() {
        expands_to("(x + y)^2", "x^2 + 2*x*y + y^2");
    }

    #[test]
    fn expand_difference_of_squares() {
        expands_to("(x+y)*(x-y)", "x^2 - y^2");
    }

    #[test]
    fn expand_collect_like_terms() {
        expands_to("(x+1)*(x+2)", "x^2 + 3*x + 2");
        expands_to("(x+1) + (x+2)", "2*x + 3");
    }

    #[test]
    fn expand_high_power() {
        // (x+1)^10: leading coeff 1, constant 1, binomial coefficients
        let s = expand_str("(x+1)^10");
        assert!(s.starts_with("x^10"));
        assert!(s.contains("45*x^8"));
        assert!(s.ends_with("+ 1"));
    }

    #[test]
    fn expand_nested_products() {
        expands_to("(x+1)*(x+2)*(x+3)", "x^3 + 6*x^2 + 11*x + 6");
    }

    #[test]
    fn expand_scalar_multiplication() {
        expands_to("3*(x + 4)", "3*x + 12");
        expands_to("2*(x+1)^2", "2*x^2 + 4*x + 2");
    }

    #[test]
    fn expand_division_by_constant() {
        expands_to("(2*x + 4)/2", "x + 2");
        expands_to("(x+1)/2*(x+3)", "0.5*x^2 + 2*x + 1.5");
    }

    #[test]
    fn expand_negative_coefficients() {
        expands_to("(x - 1)*(x + 1)", "x^2 - 1");
        expands_to("0 - (x+1)^2", "0 - x^2 - 2*x - 1");
        // Prefix unary minus applies to the complete power expression.
        expands_to("-(x + 1)^2", "-x^2 - 2*x - 1");
    }

    #[test]
    fn expand_keeps_functions_opaque() {
        // sin(x)*(x + 1): the function is not polynomial, so the product
        // stays, but the polynomial side is already collected.
        let e = expand(&Parser::parse("sin(x)*(x + 1)").unwrap());
        let Expr::Mul(_, b) = &e else { panic!("expected product") };
        assert!(b.equals(&Parser::parse("x + 1").unwrap()));

        // Arguments of functions are expanded.
        assert!(expand_str("sin((x+1)^2)").contains("x^2"));
    }

    #[test]
    fn expand_keeps_symbolic_powers() {
        // Non-numeric exponents stay symbolic.
        let s = expand_str("(x+1)^y");
        assert!(s.contains("(x + 1)^y") || s.contains("(x+1)^y") || s.contains("^y"), "got: {}", s);
    }

    #[test]
    fn expand_keeps_variable_denominators() {
        let s = expand_str("(x+1)/x");
        assert!(s.contains("/ x") || s.contains("/x"), "got: {}", s);
    }

    #[test]
    fn expand_too_large_power_stays_symbolic() {
        // Exponent beyond the cap: keep Pow form, don't explode.
        let s = expand_str("(x+1)^1000000");
        assert!(s.contains("1000000"));
        assert!(s.matches('+').count() < 100);
    }

    #[test]
    fn expand_multivariate_ordering() {
        let s = expand_str("(x + y + 1)^2");
        // Highest degree first
        assert!(s.starts_with("x^2"), "got: {}", s);
    }

    #[test]
    fn poly_round_trip_term_count() {
        // (x+2)(x+3) = x^2 + 5x + 6 → 3 collected terms
        let p = to_poly(&Parser::parse("(x+2)*(x+3)").unwrap()).unwrap();
        assert_eq!(p.len(), 3);
    }

    #[test]
    fn expand_no_op_on_monomial() {
        assert!(expand(&Parser::parse("x").unwrap()).equals(&Parser::parse("x").unwrap()));
        assert!(expand(&Parser::parse("5").unwrap()).equals(&Parser::parse("5").unwrap()));
    }

    // ===== Additional coverage =====

    #[test]
    fn to_poly_returns_none_for_function() {
        // sin(x) is not a polynomial
        let e = Parser::parse("sin(x)").unwrap();
        assert!(to_poly(&e).is_none());
    }

    #[test]
    fn to_poly_returns_none_for_variable_in_denominator() {
        // 1/x is not a polynomial
        let e = Parser::parse("1/x").unwrap();
        assert!(to_poly(&e).is_none());
    }

    #[test]
    fn to_poly_returns_none_for_symbolic_power() {
        // x^y has a symbolic exponent
        let e = Parser::parse("x^y").unwrap();
        assert!(to_poly(&e).is_none());
    }

    #[test]
    fn to_poly_returns_none_for_negative_exponent() {
        // x^(-1) is not a polynomial
        let e = Parser::parse("x^(-1)").unwrap();
        assert!(to_poly(&e).is_none());
    }

    #[test]
    fn to_poly_returns_none_for_fractional_exponent() {
        // x^(1/2) is not a polynomial
        let e = Parser::parse("x^(1/2)").unwrap();
        assert!(to_poly(&e).is_none());
    }

    #[test]
    fn poly_to_expr_round_trip() {
        // expand → to_poly → poly_to_expr → to_poly should be stable
        let e = Parser::parse("(x+1)*(x+2)*(x+3)").unwrap();
        let expanded = expand(&e);
        let p = to_poly(&expanded).unwrap();
        let back = poly_to_expr(&p);
        let p2 = to_poly(&back).unwrap();
        assert_eq!(p.len(), p2.len(), "term count should match: {:?} vs {:?}", p, p2);
    }

    #[test]
    fn poly_to_expr_empty_poly() {
        // Empty poly → should produce 0 or some representation
        let empty: Poly = vec![];
        let e = poly_to_expr(&empty);
        // Should evaluate to 0
        let ctx = crate::eval::Context::standard();
        let v = crate::eval::eval(&e, &ctx).unwrap();
        assert_eq!(v, 0.0);
    }

    #[test]
    fn expand_float_collecting() {
        // 0.1*x + 0.2*x should collect to 0.3*x (within float tolerance)
        let e = Parser::parse("0.1*x + 0.2*x").unwrap();
        let expanded = expand(&e);
        let p = to_poly(&expanded).unwrap();
        // Should be a single term with coefficient ~0.3
        assert_eq!(p.len(), 1, "should collect to one term: {:?}", p);
        assert!((p[0].coeff - 0.3).abs() < 1e-9, "coeff should be 0.3: {}", p[0].coeff);
    }

    #[test]
    fn expand_trinomial_cube() {
        // (x + y + z)^2 should have 6 terms
        let s = expand_str("(x + y + z)^2");
        let e = expand(&Parser::parse("(x + y + z)^2").unwrap());
        let p = to_poly(&e).unwrap();
        assert_eq!(p.len(), 6, "should have 6 terms: {} -> {:?}", s, p);
    }

    #[test]
    fn expand_negative_coefficient_term() {
        // (x - 2)^3 = x^3 - 6x^2 + 12x - 8
        expands_to("(x - 2)^3", "x^3 - 6*x^2 + 12*x - 8");
    }

    #[test]
    fn expand_constant_expression() {
        // (2 + 3)^2 = 25
        expands_to("(2 + 3)^2", "25");
    }

    #[test]
    fn expand_product_of_three_binomials() {
        // (x+1)(x+1)(x+1) = x^3 + 3x^2 + 3x + 1
        expands_to("(x+1)*(x+1)*(x+1)", "x^3 + 3*x^2 + 3*x + 1");
    }

    #[test]
    fn expand_with_subtraction() {
        // (x^2 - 1)*(x + 1) = x^3 + x^2 - x - 1
        expands_to("(x^2 - 1)*(x + 1)", "x^3 + x^2 - x - 1");
    }

    #[test]
    fn to_poly_simple_linear() {
        let e = Parser::parse("2*x + 3").unwrap();
        let p = to_poly(&e).unwrap();
        assert_eq!(p.len(), 2);
    }

    #[test]
    fn expand_preserves_zero_coefficient_cancellation() {
        // (x+1)*(x-1) = x^2 - 1 (no x term)
        let e = expand(&Parser::parse("(x+1)*(x-1)").unwrap());
        let p = to_poly(&e).unwrap();
        assert_eq!(p.len(), 2, "should have 2 terms (x^2 and -1): {:?}", p);
    }
}
