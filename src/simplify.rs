//! Algebraic simplification of [`Expr`] values.
//!
//! This implements a small but useful set of rewrite rules:
//! - constant folding
//! - identity flattening (`x + 0 = x`, `x * 1 = x`, `x * 0 = 0`, etc.)
//! - **coefficient folding** in `Mul`: numeric factors collapse into a single
//!   coefficient (`2 * (2 * x) -> 4 * x`, `x * 2 -> 2 * x` after sorting)
//! - **like-term collection** in `Add`: identical variable parts combine
//!   their coefficients (`2*x + 3*x -> 5*x`, `x - x -> 0`)
//! - double-negation collapse (`--x -> x`)
//! - divisor flattening (`0 / x -> 0`, `x / 1 -> x`)
//! - power identities (`x^0 -> 1`, `x^1 -> x`, `0^x -> 0`, `1^x -> 1`)
//!
//! Like-term collection handles multivariate products of bare variables only
//! (e.g. `2*x*y + 3*x*y -> 5*x*y`). Terms whose "variable part" contains a
//! function call or a non-integer power are treated as opaque and not
//! collected, so `2*sin(x) + 3*sin(x)` is left intact (would require matching
//! on the function structure).

use crate::expr::Expr;

/// Recursively simplify an expression.
pub fn simplify(e: &Expr) -> Expr {
    let e = preflatten(e);
    let e = collect(&e);
    fold_constants(&e)
}

/// Recursively flatten identities (`+0`, `*1`, `*0`, `^0`, `^1`, `--x`,
/// `/1`) on every node. After this pass the expression has no trivial
/// constants hiding in multiplicative positions; numeric coefficients that
/// multiply non-trivial terms still appear as `Num(n)` factors because the
/// like-term / coefficient-collection pass that comes next needs them.
fn preflatten(e: &Expr) -> Expr {
    match e {
        Expr::Num(_) | Expr::Var(_) => e.clone(),
        Expr::Neg(inner) => {
            let inner = preflatten(inner);
            match &inner {
                Expr::Num(n) => Expr::num(-n),
                Expr::Neg(x) => *x.clone(), // --x -> x
                _ => Expr::neg(inner),
            }
        }
        Expr::Add(a, b) => {
            let a = preflatten(a);
            let b = preflatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) => b,
                (_, Expr::Num(0.0)) => a,
                _ => Expr::add(a, b),
            }
        }
        Expr::Sub(a, b) => {
            let a = preflatten(a);
            let b = preflatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) => Expr::neg(b),
                (_, Expr::Num(0.0)) => a,
                _ => Expr::sub(a, b),
            }
        }
        Expr::Mul(a, b) => {
            let a = preflatten(a);
            let b = preflatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) | (_, Expr::Num(0.0)) => Expr::num(0.0),
                (Expr::Num(1.0), _) => b,
                (_, Expr::Num(1.0)) => a,
                (Expr::Num(-1.0), x) => Expr::neg(x.clone()),
                (x, Expr::Num(-1.0)) => Expr::neg(x.clone()),
                _ => Expr::mul(a, b),
            }
        }
        Expr::Div(a, b) => {
            let a = preflatten(a);
            let b = preflatten(b);
            match (&a, &b) {
                (Expr::Num(0.0), _) if !matches!(b, Expr::Num(0.0)) => Expr::num(0.0),
                (_, Expr::Num(1.0)) => a,
                _ => Expr::div(a, b),
            }
        }
        Expr::Pow(base, exp) => {
            let base = preflatten(base);
            let exp = preflatten(exp);
            match (&base, &exp) {
                (_, Expr::Num(0.0)) => Expr::num(1.0),
                (x, Expr::Num(1.0)) => x.clone(),
                (Expr::Num(1.0), _) => Expr::num(1.0),
                (Expr::Num(0.0), _) => Expr::num(0.0),
                _ => Expr::pow(base, exp),
            }
        }
        Expr::Func(name, args) => {
            Expr::func(name.clone(), args.iter().map(preflatten).collect())
        }
    }
}

/// Drive the simplifier over every leaf before any folding pass can interact
/// with partially-simplified subtrees.
fn collect(e: &Expr) -> Expr {
    match e {
        Expr::Num(_) | Expr::Var(_) => e.clone(),
        Expr::Neg(inner) => {
            let inner = collect(inner);
            // flip the coefficient on whatever multiplicative coefficient
            // the inner carries
            let (c, rest) = split_coefficient(&inner);
            let nc = -c;
            // -Neg(x) -> x (double negation handled in preflatten already,
            // but a collect pass may have produced Neg(Num(_) * rest))
            rebuild_term(nc, rest)
        }
        Expr::Add(a, b) => {
            let a = collect(a);
            let b = collect(b);
            collect_add(&[a, b])
        }
        Expr::Sub(a, b) => {
            let a = collect(a);
            let b = collect(b);
            // a - b  ==  a + (-b). Negation handled by rebuilding the term.
            let mut all = flatten_add(&a);
            let negated_b = {
                let (c, rest) = split_coefficient(&b);
                rebuild_term(-c, rest)
            };
            all.extend(flatten_add(&negated_b));
            collect_add(&all)
        }
        Expr::Mul(a, b) => {
            let a = collect(a);
            let b = collect(b);
            let (ac, ar) = split_coefficient(&a);
            let (bc, br) = split_coefficient(&b);
            let c = ac * bc;
            let rest = if ar == Expr::num(1.0) {
                br
            } else if br == Expr::num(1.0) {
                ar
            } else {
                Expr::mul(ar, br)
            };
            rebuild_term(c, rest)
        }
        Expr::Div(a, b) => {
            let a = collect(a);
            let b = collect(b);
            // Pull a numeric coefficient out of the numerator/denominator
            // (only safe when both are non-zero).
            let (ac, ar) = split_coefficient(&a);
            let (bc, br) = split_coefficient(&b);
            if ac == 0.0 {
                return Expr::num(0.0);
            }
            if bc == 0.0 {
                // 0/0 or x/0 — leave to the evaluator / domain error.
                return Expr::div(
                    Expr::mul(Expr::num(ac), ar),
                    Expr::mul(Expr::num(bc), br),
                );
            }
            let coe = ac / bc;
            if br == Expr::num(1.0) {
                // c * ar / 1
                return rebuild_term(coe, ar);
            }
            if ar == Expr::num(1.0) {
                // c / br
                return Expr::div(Expr::num(coe), br);
            }
            // (c * ar) / br — can't be simplified further at this layer
            // without trying to push the divisor's coefficient up too
            // (e.g. 6*x / (2*z) -> 3*x/z), but doing so requires deciding
            // what to do with cross-term denominators and we keep the
            // current behaviour conservative.
            Expr::div(rebuild_term(coe, ar), br)
        }
        Expr::Pow(base, exp) => {
            let base = collect(base);
            let exp = collect(exp);
            Expr::pow(base, exp)
        }
        Expr::Func(name, args) => {
            Expr::func(name.clone(), args.iter().map(collect).collect())
        }
    }
}

/// Flatten a `+`/`-` expression into its additive terms, with each `Sub`
/// contributing a negated right-hand side and each `Neg` pulling a minus
/// onto every inner term.
fn flatten_add(e: &Expr) -> Vec<Expr> {
    match e {
        Expr::Add(a, b) => {
            let mut v = flatten_add(a);
            v.extend(flatten_add(b));
            v
        }
        Expr::Sub(a, b) => {
            let mut v = flatten_add(a);
            v.extend(flatten_add(b).into_iter().map(|t| {
                let (c, r) = split_coefficient(&t);
                rebuild_term(-c, r)
            }));
            v
        }
        Expr::Neg(inner) => flatten_add(inner)
            .into_iter()
            .map(|t| {
                let (c, r) = split_coefficient(&t);
                rebuild_term(-c, r)
            })
            .collect(),
        _ => vec![e.clone()],
    }
}

/// Decompose an expression into a numeric coefficient `c` and a "rest"
/// such that `rebuild_term(c, rest)` matches the original (up to canonical
/// form). Numeric-only terms yield `(n, 1)`; pure-variable / function terms
/// yield `(1, e)`.
fn split_coefficient(e: &Expr) -> (f64, Expr) {
    match e {
        Expr::Num(n) => (*n, Expr::num(1.0)),
        Expr::Neg(inner) => {
            let (c, r) = split_coefficient(inner);
            (-c, r)
        }
        Expr::Mul(a, b) => {
            let (ac, ar) = split_coefficient(a);
            let (bc, br) = split_coefficient(b);
            let c = ac * bc;
            let r = if ar == Expr::num(1.0) {
                br
            } else if br == Expr::num(1.0) {
                ar
            } else {
                Expr::mul(ar, br)
            };
            (c, r)
        }
        _ => (1.0, e.clone()),
    }
}

/// Reconstruct `c * rest` in a form that respects the natural identity
/// collapse: `0 * r = 0`, `1 * r = r`, `-1 * r = -r`.
fn rebuild_term(c: f64, rest: Expr) -> Expr {
    if c == 0.0 {
        return Expr::num(0.0);
    }
    if rest == Expr::num(1.0) {
        return Expr::num(c);
    }
    if c == 1.0 {
        return rest;
    }
    if c == -1.0 {
        return Expr::neg(rest);
    }
    Expr::mul(Expr::num(c), rest)
}

/// Build an `Add` expression from a list of already-simplified terms that
/// may contain like terms. Combines same-`rest` terms by summing their
/// coefficients and keeps the result sorted by the rendered string of
/// the variable part so display is deterministic.
fn collect_add(terms: &[Expr]) -> Expr {
    let flat: Vec<Expr> = terms.iter().flat_map(flatten_add).collect();
    if flat.is_empty() {
        return Expr::num(0.0);
    }
    if flat.len() == 1 {
        return flat[0].clone();
    }
    // Group by canonical rest. Two terms "look the same" if their
    // coefficient-stripped canonical form matches; using the canonicalized
    // string keeps the matcher pure (no symbolic identity required).
    let mut groups: std::collections::BTreeMap<String, (f64, Expr)> =
        std::collections::BTreeMap::new();
    for t in &flat {
        let (c, rest) = split_coefficient(t);
        if rest == Expr::num(1.0) {
            // Pure-constant term — accumulate by keying on its numeric value.
            // Use the rendered numeric form so 2.5 and 3.5 collapse via
            // displayed-string equality rather than by canonical equivalence
            // (which can be lossy for floats).
            let key = format!("const:{c}");
            let entry = groups.entry(key).or_insert((0.0, rest));
            entry.0 += c;
        } else {
            let key = rest.canonicalize().to_string();
            let entry = groups.entry(key).or_insert((0.0, rest));
            entry.0 += c;
        }
    }
    // Now rebuild the Add as a left fold, jumping over zero-sum terms.
    // BTreeMap gives sorted iteration so the output's display is stable.
    let mut kept: Vec<(String, Expr)> = Vec::new();
    for (key, (c, rest)) in groups {
        let term = rebuild_term(c, rest);
        if matches!(&term, Expr::Num(n) if *n == 0.0) {
            continue;
        }
        kept.push((key, term));
    }
    if kept.is_empty() {
        return Expr::num(0.0);
    }
    let mut acc = kept[0].1.clone();
    for (_, t) in kept.into_iter().skip(1) {
        acc = match (&acc, &t) {
            (Expr::Num(0.0), _) => t,
            (_, Expr::Num(0.0)) => acc,
            _ => Expr::add(acc, t),
        };
    }
    acc
}

/// Fold sub-expressions that consist entirely of constants.
fn fold_constants(e: &Expr) -> Expr {
    if let Some(v) = try_eval_const(e) {
        Expr::num(v)
    } else {
        e.clone()
    }
}

fn try_eval_const(e: &Expr) -> Option<f64> {
    use Expr::*;
    match e {
        Num(n) => Some(*n),
        Var(_) => None,
        Neg(e) => try_eval_const(e).map(|x| -x),
        Add(a, b) => Some(try_eval_const(a)? + try_eval_const(b)?),
        Sub(a, b) => Some(try_eval_const(a)? - try_eval_const(b)?),
        Mul(a, b) => Some(try_eval_const(a)? * try_eval_const(b)?),
        Div(a, b) => {
            let bv = try_eval_const(b)?;
            if bv == 0.0 {
                None
            } else {
                Some(try_eval_const(a)? / bv)
            }
        }
        Pow(a, b) => Some(try_eval_const(a)?.powf(try_eval_const(b)?)),
        Func(name, args) => {
            let _ = name;
            let vs: Option<Vec<f64>> = args.iter().map(try_eval_const).collect();
            vs.and_then(|_| crate::eval::eval(e, standard_ctx()).ok())
        }
    }
}

use std::sync::OnceLock;

static STANDARD_CTX: OnceLock<crate::eval::Context> = OnceLock::new();

fn standard_ctx() -> &'static crate::eval::Context {
    STANDARD_CTX.get_or_init(crate::eval::Context::standard)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn s(src: &str) -> Expr {
        simplify(&Parser::parse(src).unwrap())
    }

    #[test]
    fn identities() {
        assert_eq!(s("x + 0"), Expr::var("x"));
        assert_eq!(s("0 + x"), Expr::var("x"));
        assert_eq!(s("x * 1"), Expr::var("x"));
        assert_eq!(s("x * 0"), Expr::num(0.0));
        assert_eq!(s("x^1"), Expr::var("x"));
        assert_eq!(s("x^0"), Expr::num(1.0));
    }

    #[test]
    fn constant_folding() {
        assert_eq!(s("2 + 3"), Expr::num(5.0));
        assert_eq!(s("2 * 3 + 4"), Expr::num(10.0));
    }

    #[test]
    fn signed_one_distribution() {
        // -1 * x should become -x
        let got = s("-1 * x");
        let want = Expr::neg(Expr::var("x"));
        assert_eq!(got, want);
    }

    #[test]
    fn function_with_constants() {
        // sin(0) -> 0
        assert_eq!(s("sin(0)"), Expr::num(0.0));
        // cos(0) -> 1
        assert_eq!(s("cos(0)"), Expr::num(1.0));
    }

    // ---- Like-term collection --------------------------------------------

    #[test]
    fn like_term_simple() {
        // 2*x + 3*x -> 5*x
        assert_eq!(s("2*x + 3*x"), Expr::mul(Expr::num(5.0), Expr::var("x")));
    }

    #[test]
    fn like_term_in_reverse_order() {
        // 3*x + 2*x -> 5*x (commutative)
        assert_eq!(s("3*x + 2*x"), Expr::mul(Expr::num(5.0), Expr::var("x")));
    }

    #[test]
    fn like_term_three_addends() {
        // 1*x + 2*x + 3*x -> 6*x
        assert_eq!(s("x + 2*x + 3*x"), Expr::mul(Expr::num(6.0), Expr::var("x")));
    }

    #[test]
    fn like_term_with_negative() {
        // 5*x - 3*x -> 2*x
        assert_eq!(s("5*x - 3*x"), Expr::mul(Expr::num(2.0), Expr::var("x")));
    }

    #[test]
    fn like_term_cancels_to_zero() {
        // x - x -> 0
        assert_eq!(s("x - x"), Expr::num(0.0));
    }

    #[test]
    fn like_term_three_var() {
        // 2*x*y + 3*x*y -> 5*x*y
        assert_eq!(
            s("2*x*y + 3*x*y"),
            Expr::mul(
                Expr::num(5.0),
                Expr::mul(Expr::var("x"), Expr::var("y"))
            )
        );
    }

    #[test]
    fn unlike_terms_not_collected() {
        // 2*x + 3*y has different rests and should stay distinct (only the
        // canonical-form Add survives; ensure we did not lose factors).
        let got = s("2*x + 3*y");
        assert!(matches!(got, Expr::Add(_, _)));
        // restore canonical sort so the assertion is stable
        let canonical = got.canonicalize();
        assert_eq!(
            canonical,
            Expr::add(
                Expr::mul(Expr::num(2.0), Expr::var("x")),
                Expr::mul(Expr::num(3.0), Expr::var("y"))
            )
        );
    }

    #[test]
    fn like_term_with_constant_mixed_in() {
        // 2*x + 3 + 5*x -> 7*x + 3 (coefficients + constant folded separately)
        let got = s("2*x + 3 + 5*x").canonicalize();
        assert_eq!(
            got,
            Expr::add(
                Expr::num(3.0),
                Expr::mul(Expr::num(7.0), Expr::var("x"))
            )
        );
    }

    // ---- Coefficient collapse in Mul -------------------------------------

    #[test]
    fn mul_two_numerics_collapses() {
        // 2 * 3 -> 6
        assert_eq!(s("2 * 3"), Expr::num(6.0));
    }

    #[test]
    fn mul_numeric_through_parens_collapses() {
        // 2 * (2 * x) -> 4*x
        assert_eq!(s("2 * (2 * x)"), Expr::mul(Expr::num(4.0), Expr::var("x")));
    }

    #[test]
    fn mul_chained_collapses() {
        // 2*3*x -> 6*x
        assert_eq!(s("2*3*x"), Expr::mul(Expr::num(6.0), Expr::var("x")));
    }

    #[test]
    fn mul_negative_coefficient() {
        // -2 * 3 * x -> -6*x (the minus rides on the numeric coefficient)
        assert_eq!(
            s("-2 * 3 * x"),
            Expr::mul(Expr::num(-6.0), Expr::var("x"))
        );
    }

    #[test]
    fn mul_signed_one_stripped() {
        // -1 * (2 * x) -> -(2*x) — the -1 fold into the inner coefficient
        // (-2*x).
        assert_eq!(
            s("-1 * (2 * x)"),
            Expr::mul(Expr::num(-2.0), Expr::var("x"))
        );
    }

    // ---- Neg / Sub arithmetic --------------------------------------------

    #[test]
    fn double_negation_collapses() {
        // -(-x) -> x
        assert_eq!(s("-(-x)"), Expr::var("x"));
    }

    #[test]
    fn neg_negative_constant_collapses() {
        // -(-5) -> 5
        assert_eq!(s("-(-5)"), Expr::num(5.0));
    }

    #[test]
    fn sub_flips_to_neg_in_add() {
        // 2*x - 3*x -> -x (Neg form is the canonical restructure of -1*x)
        assert_eq!(s("2*x - 3*x"), Expr::neg(Expr::var("x")));
    }

    #[test]
    fn nested_sub_cancel() {
        // x - y - (x - y) -> 0
        assert_eq!(s("x - y - (x - y)"), Expr::num(0.0));
    }

    // ---- Mul × Add interaction -------------------------------------------

    #[test]
    fn mul_distributing_does_not_collect_across_mul() {
        // a is still opaque, but the inner (2+3) folds to 5
        // (2 + 3) -> 5 (no x present)
        assert_eq!(s("2 + 3"), Expr::num(5.0));
        // x*(2+3)*x — only the constant Add inside folds; x*5*x stays as
        // a product of two variables after coefficient collection.
        let got = s("x*(2+3)*x");
        assert!(matches!(got, Expr::Mul(_, _) | Expr::Num(_)));
    }

    // ---- Edge cases ------------------------------------------------------

    #[test]
    fn like_term_factor_preserved() {
        // 2*x*y - x*y = x*y (coefficient 1 collapsed by rebuild_term)
        assert_eq!(
            s("2*x*y - x*y"),
            Expr::mul(Expr::var("x"), Expr::var("y"))
        );
    }

    #[test]
    fn nested_like_terms_collect() {
        // (2*x + 3*x) + x -> 6*x
        assert_eq!(
            s("(2*x + 3*x) + x"),
            Expr::mul(Expr::num(6.0), Expr::var("x"))
        );
    }

    #[test]
    fn polynomial_like_terms_collect() {
        // 2*x^2 + 3*x^2 + 5*x + 7*x -> 5*x^2 + 12*x
        let got = s("2*x^2 + 3*x^2 + 5*x + 7*x").canonicalize();
        // The exact ordering is canonical-form dependent; just make sure
        // the coefficients are right by checking equality after
        // canonicalization.
        let want = Expr::add(
            Expr::mul(Expr::num(5.0), Expr::pow(Expr::var("x"), Expr::num(2.0))),
            Expr::mul(Expr::num(12.0), Expr::var("x")),
        );
        let want_c = want.canonicalize();
        assert!(got.equals(&want_c), "got {got}, want {want_c}");
    }
}
