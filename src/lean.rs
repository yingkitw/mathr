//! Lean 4 export for [`Expr`] — renders expressions as Lean terms over
//! `ℝ` (requires **Mathlib** to compile).
//!
//! Mapping conventions:
//! - numeric literals → `({n} : ℝ)` (integers as `({n} : ℝ)`, decimals via
//!   Lean's `OfScientific` — the exported literal is *exactly* the f64's
//!   decimal value, so it is sound to paste into a proof)
//! - `pi`/`e`/`tau` parse-folded to f64 bits are exported as those decimal
//!   literals, **not** as `Real.pi`/`Real.exp 1` — the f64 bits are a
//!   rational approximation and upgrading them silently would be unsound
//! - `x^2` with an integral exponent → `(x ^ (2 : ℕ))` (monoid pow — correct
//!   for negative bases); non-integral exponents → rpow `(x ^ (y : ℝ))`
//! - elementary functions map to Mathlib names (`sin` → `Real.sin`,
//!   `ln` → `Real.log`, `log(x, b)` → `Real.logb b x` — note the argument
//!   swap, `log2`/`log10` → `Real.logb (2:ℝ) x` / `… (10:ℝ) x`,
//!   `floor`/`ceil` → coerced `↑⌊x⌋`/`↑⌈x⌉` since Mathlib returns `ℤ`)
//! - unmapped functions (e.g. `gamma`, `erf`, `bessel_j0`) are emitted as
//!   plain applications `(gamma x)` — define them in Lean before compiling
//!
//! Variable names are emitted as-is: the mathr parser only accepts
//! `[a-zA-Z_][a-zA-Z0-9_]*`, which is a subset of valid Lean identifiers.
//!
//! `NaN`/`±inf` have no `ℝ` counterpart and are a [`MathError::Domain`].
//!
//! # Example
//!
//! ```
//! use mathr::lean::to_lean;
//! use mathr::parser::Parser;
//!
//! let e = Parser::parse("x^2 + sin(y)").unwrap();
//! assert_eq!(to_lean(&e).unwrap(), "((x ^ (2 : ℕ)) + (Real.sin y))");
//! ```

use crate::error::{MathError, Result};
use crate::expr::Expr;

/// Export `e` as a Lean 4 term (requires Mathlib in scope to compile).
pub fn to_lean(e: &Expr) -> Result<String> {
    emit(e)
}

/// Export `e` as a self-contained Lean file: `import Mathlib` plus a
/// `noncomputable def` with every free variable bound as `(name : ℝ)`.
pub fn to_lean_doc(e: &Expr) -> Result<String> {
    let term = emit(e)?;
    let vars = e.variables();
    let binders = if vars.is_empty() {
        String::new()
    } else {
        format!(" ({} : ℝ)", vars.join(" "))
    };
    Ok(format!(
        "import Mathlib\n\n/-- Exported from mathr. -/\nnoncomputable def mathr_export{binders} : ℝ :=\n  {term}\n"
    ))
}

fn emit(e: &Expr) -> Result<String> {
    match e {
        Expr::Num(n) => {
            if !n.is_finite() {
                return Err(MathError::Domain(format!(
                    "cannot export {n} to Lean: ℝ has no NaN/inf"
                )));
            }
            let v = if *n == 0.0 { 0.0 } else { *n }; // normalise -0.0
            if v.fract() == 0.0 && v.abs() <= 9.0e15 {
                Ok(format!("({} : ℝ)", v as i64))
            } else {
                Ok(format!("({v} : ℝ)"))
            }
        }
        Expr::Var(v) => Ok(v.clone()),
        Expr::Neg(a) => Ok(format!("(-{})", emit(a)?)),
        Expr::Add(a, b) => Ok(format!("({} + {})", emit(a)?, emit(b)?)),
        Expr::Sub(a, b) => Ok(format!("({} - {})", emit(a)?, emit(b)?)),
        Expr::Mul(a, b) => Ok(format!("({} * {})", emit(a)?, emit(b)?)),
        Expr::Div(a, b) => Ok(format!("({} / {})", emit(a)?, emit(b)?)),
        Expr::Pow(a, b) => emit_pow(a, b),
        Expr::Func(name, args) => emit_func(name, args),
    }
}

/// `x^n` — integral non-negative exponents use monoid `^` (ℕ), which is
/// total and correct for negative bases; anything else is rpow (ℝ).
fn emit_pow(a: &Expr, b: &Expr) -> Result<String> {
    let base = emit(a)?;
    if let Expr::Num(e) = b {
        if e.fract() == 0.0 && *e >= 0.0 && *e <= 1.0e9 {
            return Ok(format!("({base} ^ ({} : ℕ))", *e as u64));
        }
    }
    Ok(format!("({base} ^ ({} : ℝ))", emit(b)?))
}

fn emit_func(name: &str, args: &[Expr]) -> Result<String> {
    // binary min/max fold (mathr is variadic, Lean is binary)
    if name == "min" || name == "max" {
        match args.len() {
            0 => return Err(MathError::Eval(format!("{name} needs at least one arg"))),
            // mathr's eval folds a single-arg min/max to the arg itself;
            // emitting `(min x)` would be invalid Lean (min is binary).
            1 => return emit(&args[0]),
            _ => {
                let mut acc = emit(&args[args.len() - 1])?;
                for a in args[..args.len() - 1].iter().rev() {
                    acc = format!("({name} {} {acc})", emit(a)?);
                }
                return Ok(acc);
            }
        }
    }
    if name == "pow" && args.len() == 2 {
        return emit_pow(&args[0], &args[1]);
    }
    if name == "log" && args.len() == 2 {
        // mathr log(x, b); Mathlib Real.logb b x — arguments swap sides.
        return Ok(format!("(Real.logb {} {})", emit(&args[1])?, emit(&args[0])?));
    }
    if args.len() == 1 {
        let a = emit(&args[0])?;
        let out = match name {
            "sin" => format!("(Real.sin {a})"),
            "cos" => format!("(Real.cos {a})"),
            "tan" => format!("(Real.tan {a})"),
            "asin" => format!("(Real.arcsin {a})"),
            "acos" => format!("(Real.arccos {a})"),
            "atan" => format!("(Real.arctan {a})"),
            "sinh" => format!("(Real.sinh {a})"),
            "cosh" => format!("(Real.cosh {a})"),
            "tanh" => format!("(Real.tanh {a})"),
            "exp" => format!("(Real.exp {a})"),
            "ln" => format!("(Real.log {a})"),
            "log2" => format!("(Real.logb (2 : ℝ) {a})"),
            "log10" => format!("(Real.logb (10 : ℝ) {a})"),
            "sqrt" => format!("(Real.sqrt {a})"),
            "abs" => format!("(abs {a})"),
            "round" => format!("(round {a})"),
            "fract" => format!("(Int.fract {a})"),
            // Mathlib's floor/ceil land in ℤ; coerce back to ℝ.
            "floor" => format!("(↑⌊{a}⌋ : ℝ)"),
            "ceil" => format!("(↑⌈{a}⌉ : ℝ)"),
            // Unmapped (gamma, erf, bessel_*, zeta, …): emit a plain
            // application; define `name : ℝ → ℝ` in Lean before compiling.
            other => format!("({other} {a})"),
        };
        return Ok(out);
    }
    // General fallback: plain application of however many args.
    let mut out = format!("({name}");
    for a in args {
        out.push(' ');
        out.push_str(&emit(a)?);
    }
    out.push(')');
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn l(src: &str) -> String {
        let e = crate::parser::Parser::parse(src).unwrap();
        to_lean(&e).unwrap()
    }

    #[test]
    fn num_integer() {
        assert_eq!(l("2"), "(2 : ℝ)");
    }

    #[test]
    fn num_decimal_exact() {
        assert_eq!(l("0.5"), "(0.5 : ℝ)");
    }

    #[test]
    fn num_negative_is_neg_node() {
        // parser: -x binds as Neg; -2 stays Neg(Num(2))
        assert_eq!(l("-2"), "(-(2 : ℝ))");
    }

    #[test]
    fn var_passthrough() {
        assert_eq!(l("x"), "x");
        assert_eq!(l("foo_1"), "foo_1");
    }

    #[test]
    fn add() {
        assert_eq!(l("x + 1"), "(x + (1 : ℝ))");
    }

    #[test]
    fn mul_with_coefficient() {
        assert_eq!(l("2*x"), "((2 : ℝ) * x)");
    }

    #[test]
    fn sub_and_div() {
        assert_eq!(l("x - y"), "(x - y)");
        assert_eq!(l("x/y"), "(x / y)");
    }

    #[test]
    fn pow_integral_uses_nat() {
        assert_eq!(l("x^2"), "(x ^ (2 : ℕ))");
    }

    #[test]
    fn pow_symbolic_uses_real_rpow() {
        assert_eq!(l("x^y"), "(x ^ (y : ℝ))");
    }

    #[test]
    fn pow_negative_base_integral_exponent() {
        // (-2)^3 (parser binds the minus into the base) must NOT be rpow —
        // rpow has junk on negative bases, monoid pow is exact.
        assert_eq!(l("(-2)^3"), "((-(2 : ℝ)) ^ (3 : ℕ))");
    }

    #[test]
    fn elementary_functions() {
        assert_eq!(l("sin(x)"), "(Real.sin x)");
        assert_eq!(l("cos(y)"), "(Real.cos y)");
        assert_eq!(l("ln(x)"), "(Real.log x)");
        assert_eq!(l("exp(x)"), "(Real.exp x)");
        assert_eq!(l("sqrt(x)"), "(Real.sqrt x)");
        assert_eq!(l("atan(x)"), "(Real.arctan x)");
    }

    #[test]
    fn log_argument_swap() {
        // mathr log(x, b) = Mathlib Real.logb b x
        assert_eq!(l("log(x, 10)"), "(Real.logb (10 : ℝ) x)");
    }

    #[test]
    fn log2_log10() {
        assert_eq!(l("log2(x)"), "(Real.logb (2 : ℝ) x)");
        assert_eq!(l("log10(x)"), "(Real.logb (10 : ℝ) x)");
    }

    #[test]
    fn floor_ceil_coerced_to_real() {
        assert_eq!(l("floor(x)"), "(↑⌊x⌋ : ℝ)");
        assert_eq!(l("ceil(x)"), "(↑⌈x⌉ : ℝ)");
    }

    #[test]
    fn abs_round_fract() {
        assert_eq!(l("abs(x)"), "(abs x)");
        assert_eq!(l("round(x)"), "(round x)");
        assert_eq!(l("fract(x)"), "(Int.fract x)");
    }

    #[test]
    fn min_max_folding() {
        assert_eq!(l("min(1, 2, 3)"), "(min (1 : ℝ) (min (2 : ℝ) (3 : ℝ)))");
        assert_eq!(l("max(x, y)"), "(max x y)");
    }

    #[test]
    fn unmapped_function_plain_application() {
        assert_eq!(l("gamma(x)"), "(gamma x)");
        assert_eq!(l("erf(0.5)"), "(erf (0.5 : ℝ))");
    }

    #[test]
    fn pow_func_matches_pow_operator() {
        assert_eq!(l("pow(x, 2)"), "(x ^ (2 : ℕ))");
    }

    #[test]
    fn combined_expression() {
        assert_eq!(l("x^2 + sin(y)"), "((x ^ (2 : ℕ)) + (Real.sin y))");
    }

    #[test]
    fn nan_rejected() {
        let e = Expr::num(f64::NAN);
        assert!(to_lean(&e).is_err());
        let e = Expr::num(f64::INFINITY);
        assert!(to_lean(&e).is_err());
    }

    #[test]
    fn doc_binds_variables() {
        let e = crate::parser::Parser::parse("2*x + y").unwrap();
        let doc = to_lean_doc(&e).unwrap();
        assert!(doc.starts_with("import Mathlib\n"), "got: {doc}");
        assert!(
            doc.contains("noncomputable def mathr_export (x y : ℝ) : ℝ :="),
            "got: {doc}"
        );
        assert!(doc.contains("((2 : ℝ) * x) + y"), "got: {doc}");
    }

    #[test]
    fn doc_no_variables() {
        let e = crate::parser::Parser::parse("2 + 3").unwrap();
        let doc = to_lean_doc(&e).unwrap();
        assert!(
            doc.contains("noncomputable def mathr_export : ℝ :="),
            "got: {doc}"
        );
    }

    // ---- Extended capability validation ----------------------------------

    /// Structural invariants every exported Lean term must satisfy
    /// (Lean isn't installed in CI, so this is the compile-proxy check).
    fn assert_wellformed(term: &str) {
        // balanced parentheses, never negative depth
        let mut depth: i32 = 0;
        for c in term.chars() {
            match c {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    assert!(depth >= 0, "unbalanced ')' in {term}");
                }
                _ => {}
            }
        }
        assert_eq!(depth, 0, "unbalanced '(' in {term}");
        // ℝ has no NaN/inf and Lean idents can't start with a digit
        assert!(!term.contains("NaN"), "NaN leaked into {term}");
        assert!(!term.contains("inf"), "inf leaked into {term}");
        // every numeral must be type-ascribed (bare numerals default to ℕ);
        // ℕ ascriptions are only legal as pow exponents, ℝ everywhere else
        let bytes = term.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            // numeral only where an expression can start: after '(', ' ', '-'
            // — a digit after an ident char (e.g. `bessel_j0`) is a name, not a numeral
            let numeral_start = bytes[i].is_ascii_digit()
                && (i == 0
                    || !(bytes[i - 1].is_ascii_digit()
                        || bytes[i - 1] == b'.'
                        || bytes[i - 1].is_ascii_alphanumeric()
                        || bytes[i - 1] == b'_'));
            if numeral_start {
                let mut j = i;
                while j < bytes.len() && (bytes[j].is_ascii_digit() || bytes[j] == b'.') {
                    j += 1;
                }
                let tail = &term[j..];
                assert!(
                    tail.starts_with(" : ℝ)") || tail.starts_with(" : ℕ)"),
                    "unascribed numeral in {term}"
                );
                i = j;
            } else {
                i += 1;
            }
        }
    }

    /// Golden input→output pairs; the workhorse capability table.
    #[test]
    fn golden_corpus() {
        let cases: [(&str, &str); 21] = [
            // associativity / precedence is preserved by explicit parens
            ("x - y - z", "((x - y) - z)"),
            ("x - (y - z)", "(x - (y - z))"),
            ("x/y/z", "((x / y) / z)"),
            ("x/(y + z)", "(x / (y + z))"),
            ("x/(y*z)", "(x / (y * z))"),
            ("(x + y)^2", "((x + y) ^ (2 : ℕ))"),
            ("(x + y)*(x - y)", "((x + y) * (x - y))"),
            // right-associative power tower: inner ℕ, outer ℝ (coerced)
            ("2^3^2", "((2 : ℝ) ^ (((3 : ℝ) ^ (2 : ℕ)) : ℝ))"),
            // nested functions
            ("sin(cos(ln(x)))", "(Real.sin (Real.cos (Real.log x)))"),
            // min/max folding edge cases
            ("min(x)", "x"),
            ("max(1, 2, 3, 4)", "(max (1 : ℝ) (max (2 : ℝ) (max (3 : ℝ) (4 : ℝ))))"),
            // div by zero literal: Lean field junk, but it compiles
            ("x/0", "(x / (0 : ℝ))"),
            // the `i` literal (ceval's imaginary unit) exports as an ident —
            // complex expressions belong to cval, documented limitation
            ("i + x", "(i + x)"),
            // large integers beyond the i64 fast path still print as digits
            ("10000000000000000", "(10000000000000000 : ℝ)"),
            // tiny decimals print in plain notation (no Rust exponent)
            ("0.00001", "(0.00001 : ℝ)"),
            // negative decimal rides a Neg node
            ("-2.5", "(-(2.5 : ℝ))"),
            // two-variable mixed expression
            ("x*y + y*x", "((x * y) + (y * x))"),
            // power of a function
            ("sin(x)^2", "((Real.sin x) ^ (2 : ℕ))"),
            // subtraction inside a denominator keeps parens (re-parse safe)
            ("1/(x - 1)", "((1 : ℝ) / (x - (1 : ℝ)))"),
            // unary minus on a product
            ("-(x*y)", "(-(x * y))"),
            // deep nesting does not lose structure
            ("((x + 1)*(x + 2))/((x + 3)*(x + 4))",
             "(((x + (1 : ℝ)) * (x + (2 : ℝ))) / ((x + (3 : ℝ)) * (x + (4 : ℝ))))"),
        ];
        for (src, want) in cases {
            let got = l(src);
            assert_eq!(got, want, "mismatch for {src}");
            assert_wellformed(&got);
        }
    }

    /// π / e / τ are parse-folded to f64 bits; export must keep the EXACT
    /// decimal (sound), not upgrade to Real.pi (which is a different real).
    #[test]
    fn constants_export_as_exact_decimals() {
        assert_eq!(l("pi"), "(3.141592653589793 : ℝ)");
        assert_eq!(l("e"), "(2.718281828459045 : ℝ)");
        assert_eq!(l("tau"), "(6.283185307179586 : ℝ)");
        // `e^(i*pi)` — i exports as a plain identifier (complex belongs to cval)
        assert_eq!(
            l("e^(i*pi) + 1"),
            "(((2.718281828459045 : ℝ) ^ ((i * (3.141592653589793 : ℝ)) : ℝ)) + (1 : ℝ))"
        );
    }

    /// -0.0 normalises to `(0 : ℝ)`; huge/small magnitudes stay literal.
    #[test]
    fn numeric_edge_cases() {
        assert_eq!(to_lean(&Expr::num(-0.0)).unwrap(), "(0 : ℝ)");
        assert_eq!(to_lean(&Expr::num(1.0e16)).unwrap(), "(10000000000000000 : ℝ)");
        assert_eq!(to_lean(&Expr::num(2.5e-7)).unwrap(), "(0.00000025 : ℝ)");
        // i64::MAX rounds to 2^63 in f64; Rust's shortest repr keeps the
        // rounded digits — still an exact, valid Lean numeral
        assert_eq!(to_lean(&Expr::num(i64::MAX as f64)).unwrap(), "(9223372036854776000 : ℝ)");
    }

    /// min/max arity validation.
    #[test]
    fn min_max_arity() {
        assert_eq!(l("min(x)"), "x");
        assert_eq!(l("max(y)"), "y");
        let e = Expr::func("min", vec![]);
        assert!(to_lean(&e).is_err());
    }

    /// Wellformedness must hold for every mapped elementary function.
    #[test]
    fn all_mapped_functions_wellformed() {
        for src in [
            "sin(x)", "cos(x)", "tan(x)", "asin(x)", "acos(x)", "atan(x)",
            "sinh(x)", "cosh(x)", "tanh(x)", "exp(x)", "ln(x)", "log2(x)",
            "log10(x)", "sqrt(x)", "abs(x)", "round(x)", "fract(x)",
            "floor(x)", "ceil(x)", "log(x, 3)", "min(x, y)", "max(1, 2, 3)",
            "pow(x, 3)", "pow(x, y)", "gamma(x)", "erf(x)", "bessel_j0(x)",
            "zeta(x)", "mod(x, 2)", "sign(x)", "cbrt(x)",
        ] {
            let got = l(src);
            assert_wellformed(&got);
        }
    }

    /// to_lean_doc output must also satisfy the invariants and bind
    /// variables in first-seen order.
    #[test]
    fn doc_wellformed_and_ordered() {
        let e = crate::parser::Parser::parse("y^2 + x + sin(y)").unwrap();
        let doc = to_lean_doc(&e).unwrap();
        assert_wellformed(doc.lines().last().unwrap());
        // y first (first seen), x second
        assert!(
            doc.contains("noncomputable def mathr_export (y x : ℝ) : ℝ :="),
            "got: {doc}"
        );
    }
}
