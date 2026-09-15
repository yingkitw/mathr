//! Complex-number evaluation of `Expr` ASTs.
//!
//! Zero changes to the parser or evaluator: `i` parses as `Var("i")` (implicit
//! multiplication makes `1 + 2i` parse naturally) and is interpreted here as
//! the imaginary unit; every other variable is looked up in the [`Context`]
//! and treated as a real value. Constants `pi`/`e`/`tau` are folded to `f64`
//! at parse time — unlike `bigdec.rs` no mapping-back is needed because `f64`
//! is this evaluator's precision ceiling.
//!
//! Branch conventions: principal branch throughout (`ln(-1) = i*pi`,
//! `sqrt(-1) = i`, `arg` in `(-pi, pi]`). Integer powers are evaluated by
//! repeated multiplication (exact for small exponents); general powers use
//! `exp(w * ln z)`.

use crate::complex::Complex;
use crate::error::{MathError, Result};
use crate::eval::Context;
use crate::expr::Expr;

const I: Complex<f64> = Complex::new(0.0, 1.0);
const TWO_I: Complex<f64> = Complex::new(0.0, 2.0);

/// Evaluate an expression over the complex numbers.
pub fn eval_complex(expr: &Expr, ctx: &Context) -> Result<Complex<f64>> {
    let z = eval_complex_inner(expr, ctx)?;
    // IEEE signed zero: `Neg` yields -0.0 components, and atan2(-0.0, x < 0)
    // returns -pi — turning sqrt(-1) into -i. Canonicalize to +0.0.
    Ok(Complex::new(z.re + 0.0, z.im + 0.0))
}

fn eval_complex_inner(expr: &Expr, ctx: &Context) -> Result<Complex<f64>> {
    match expr {
        Expr::Num(x) => Ok(Complex::new(*x, 0.0)),
        Expr::Var(v) => match v.as_str() {
            "i" => Ok(I),
            name => ctx
                .vars
                .get(name)
                .map(|&r| Complex::new(r, 0.0))
                .ok_or_else(|| MathError::UnknownVariable(name.to_string())),
        },
        Expr::Neg(a) => Ok(-eval_complex(a, ctx)?),
        Expr::Add(a, b) => Ok(eval_complex(a, ctx)? + eval_complex(b, ctx)?),
        Expr::Sub(a, b) => Ok(eval_complex(a, ctx)? - eval_complex(b, ctx)?),
        Expr::Mul(a, b) => Ok(eval_complex(a, ctx)? * eval_complex(b, ctx)?),
        Expr::Div(a, b) => {
            let d = eval_complex(b, ctx)?;
            if d == Complex::ZERO {
                return Err(MathError::Domain("division by zero".into()));
            }
            Ok(eval_complex(a, ctx)? / d)
        }
        Expr::Pow(a, b) => c_pow(&eval_complex(a, ctx)?, &eval_complex(b, ctx)?),
        Expr::Func(name, args) => eval_cfunc(name, args, ctx),
    }
}

/// Parse + evaluate an expression string over the complex numbers,
/// with optional real variable bindings (`with x=2, y=-1`).
pub fn eval_complex_str(input: &str) -> Result<Complex<f64>> {
    let (expr_src, assigns) = match input.rfind(" with ") {
        Some(pos) => (&input[..pos], &input[pos + 6..]), // 6 = " with ".len()
        None => (input, ""),
    };
    let e = crate::parser::Parser::parse(expr_src.trim())?;
    let mut ctx = Context::standard();
    for a in assigns.split(',') {
        let a = a.trim();
        if a.is_empty() {
            continue;
        }
        let (name, val_src) = a
            .split_once('=')
            .ok_or_else(|| MathError::Eval(format!("bad assignment (want var=value): {a}")))?;
        let val = crate::eval::eval(&crate::parser::Parser::parse(val_src.trim())?, &ctx)?;
        ctx.set(name.trim(), val);
    }
    eval_complex(&e, &ctx)
}

fn c_ln(z: &Complex<f64>) -> Result<Complex<f64>> {
    if *z == Complex::ZERO {
        return Err(MathError::Domain("ln(0) is undefined".into()));
    }
    Ok(Complex::new(z.abs().ln(), z.arg()))
}

fn c_sqrt(z: &Complex<f64>) -> Complex<f64> {
    if *z == Complex::ZERO {
        return Complex::ZERO;
    }
    Complex::from_polar(z.abs().sqrt(), z.arg() / 2.0)
}

fn c_pow(z: &Complex<f64>, w: &Complex<f64>) -> Result<Complex<f64>> {
    // integer exponents: exact repeated multiplication (|n| <= 4096)
    if w.im == 0.0 && w.re.fract() == 0.0 && w.re.abs() <= 4096.0 {
        let n = w.re as i64;
        if n == 0 {
            return Ok(Complex::ONE); // 0^0 = 1 convention (matches bigdec)
        }
        if *z == Complex::ZERO {
            if n > 0 {
                return Ok(Complex::ZERO);
            }
            return Err(MathError::Domain("0 to a negative power".into()));
        }
        let mut acc = Complex::ONE;
        for _ in 0..n.abs() {
            acc = acc * *z;
        }
        if n < 0 {
            acc = Complex::ONE / acc;
        }
        return Ok(acc);
    }
    if *z == Complex::ZERO {
        return Err(MathError::Domain("0^w for non-integer w".into()));
    }
    // general principal-branch power: z^w = exp(w * ln z)
    Ok((*w * c_ln(z)?).exp())
}

fn need_1_arg<'a>(name: &str, args: &'a [Expr]) -> Result<&'a Expr> {
    if args.len() != 1 {
        return Err(MathError::InvalidArgument(format!(
            "{name}() expects exactly 1 argument, got {}",
            args.len()
        )));
    }
    Ok(&args[0])
}

fn eval_cfunc(name: &str, args: &[Expr], ctx: &Context) -> Result<Complex<f64>> {
    let a = eval_complex(need_1_arg(name, args)?, ctx)?;
    match name {
        "sqrt" => Ok(c_sqrt(&a)),
        "exp" => Ok(a.exp()),
        "ln" => c_ln(&a),
        "log" | "log10" => Ok(c_ln(&a)? / Complex::new(std::f64::consts::LN_10, 0.0)),
        "log2" => Ok(c_ln(&a)? / Complex::new(std::f64::consts::LN_2, 0.0)),
        "sin" => Ok(c_sin_cos(&a).0),
        "cos" => Ok(c_sin_cos(&a).1),
        "tan" => {
            let (s, c) = c_sin_cos(&a);
            if c == Complex::ZERO {
                return Err(MathError::Domain("tan: cos(z) = 0".into()));
            }
            Ok(s / c)
        }
        "asin" => c_asin(&a),
        "acos" => Ok(Complex::new(std::f64::consts::FRAC_PI_2, 0.0) - c_asin(&a)?),
        "atan" => {
            let t = (Complex::ONE + I * a) / (Complex::ONE - I * a);
            Ok(Complex::new(-0.5, 0.0) * I * c_ln(&t)?)
        }
        "sinh" => Ok((a.exp() - (-a).exp()) / Complex::new(2.0, 0.0)),
        "cosh" => Ok((a.exp() + (-a).exp()) / Complex::new(2.0, 0.0)),
        "tanh" => {
            let s = (a.exp() - (-a).exp()) / Complex::new(2.0, 0.0);
            let c = (a.exp() + (-a).exp()) / Complex::new(2.0, 0.0);
            if c == Complex::ZERO {
                return Err(MathError::Domain("tanh: cosh(z) = 0".into()));
            }
            Ok(s / c)
        }
        "abs" => Ok(Complex::new(a.abs(), 0.0)),
        other => Err(MathError::Eval(format!(
            "function '{other}' is not supported in complex evaluation (supported: sqrt, exp, ln, log, log10, log2, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, abs)"
        ))),
    }
}

/// Principal-branch arcsine: asin(z) = -i * ln(i*z + sqrt(1 - z^2)).
fn c_asin(z: &Complex<f64>) -> Result<Complex<f64>> {
    let t = I * *z + c_sqrt(&(Complex::ONE - *z * *z));
    Ok(-I * c_ln(&t)?)
}

fn c_sin_cos(z: &Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    let s = ((I * *z).exp() - (I * -*z).exp()) / TWO_I;
    let c = ((I * *z).exp() + (I * -*z).exp()) / Complex::new(2.0, 0.0);
    (s, c)
}

/// Pretty-print a complex number, snapping tiny components to zero (so
/// `exp(i*pi)` renders `-1`, not `-1 + 1.2e-16i`). Matches the REPL's
/// real-value formatting style per component.
pub fn format_complex(z: Complex<f64>) -> String {
    let eps = 1e-12;
    let (re, im) = (z.re + 0.0, z.im + 0.0);
    let (re, im) = if im.abs() < eps { (re, 0.0) } else { (re, im) };
    let comp = |v: f64| -> String {
        if v == v.trunc() && v.abs() < 1e15 {
            format!("{}", v as i64)
        } else {
            format!("{:.10}", v).trim_end_matches('0').trim_end_matches('.').to_string()
        }
    };
    if im == 0.0 {
        comp(re)
    } else if re.abs() < eps {
        match im {
            x if x == 1.0 => "i".into(),
            x if x == -1.0 => "-i".into(),
            x => format!("{}i", comp(x)),
        }
    } else if im < 0.0 {
        format!("{} - {}i", comp(re), comp(-im))
    } else {
        format!("{} + {}i", comp(re), comp(im))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::Parser;

    fn cval(s: &str) -> Complex<f64> {
        eval_complex_str(s).unwrap()
    }

    fn close(z: Complex<f64>, re: f64, im: f64) -> bool {
        (z.re - re).abs() < 1e-9 && (z.im - im).abs() < 1e-9
    }

    #[test]
    fn arithmetic() {
        // (1 + 2i) * (3 - i) = 5 + 5i
        assert!(close(cval("(1 + 2i) * (3 - i)"), 5.0, 5.0));
        // (1 + i) / (1 - i) = i
        assert!(close(cval("(1 + i) / (1 - i)"), 0.0, 1.0));
        assert!(close(cval("-3i"), 0.0, -3.0));
    }

    #[test]
    fn integer_powers() {
        assert!(close(cval("i^2"), -1.0, 0.0));
        assert!(close(cval("i^4"), 1.0, 0.0));
        assert!(close(cval("(1 + i)^2"), 0.0, 2.0));
        assert!(close(cval("i^(-1)"), 0.0, -1.0));
        assert!(close(cval("0^3"), 0.0, 0.0));
        assert!(close(cval("0^0"), 1.0, 0.0));
    }

    #[test]
    fn euler_identity() {
        // exp(i*pi) = -1
        let z = cval("exp(i*pi)");
        assert!(close(z, -1.0, 0.0), "got {z}");
        // e^(i*pi) via power operator
        let z2 = cval("e^(i*pi)");
        assert!(close(z2, -1.0, 0.0), "got {z2}");
    }

    #[test]
    fn principal_branches() {
        assert!(close(cval("sqrt(-1)"), 0.0, 1.0));
        assert!(close(cval("sqrt(-4)"), 0.0, 2.0));
        assert!(close(cval("ln(-1)"), 0.0, std::f64::consts::PI));
        // ln(e^2) = 2
        assert!(close(cval("ln(e^2)"), 2.0, 0.0));
        // 4^0.5 = 2 via the general power path
        assert!(close(cval("4^0.5"), 2.0, 0.0));
    }

    #[test]
    fn trig_and_hyperbolic() {
        // sin(2i) = i*sinh(2)
        let sinh2 = 2.0_f64.sinh();
        assert!(close(cval("sin(2i)"), 0.0, sinh2));
        // cos(i*x) = cosh(x)
        assert!(close(cval("cos(i*1.5)"), 1.5_f64.cosh(), 0.0));
        // tan(0) = 0
        assert!(close(cval("tan(0)"), 0.0, 0.0));
        // asin(2) = pi/2 - i*ln(2 + sqrt(3)) (principal)
        let expected_im = -(2.0 + 3.0_f64.sqrt()).ln();
        assert!(close(cval("asin(2)"), std::f64::consts::FRAC_PI_2, expected_im));
        // atan(1) = pi/4
        assert!(close(cval("atan(1)"), std::f64::consts::FRAC_PI_4, 0.0));
        // tanh(ln(2)) = 3/5
        assert!(close(cval("tanh(ln(2))"), 0.6, 0.0));
    }

    #[test]
    fn variables_from_context() {
        // " with x=2, y=-1": x*i - y = 1 + 2i
        assert!(close(cval("x*i - y with x=2, y=-1"), 1.0, 2.0));
    }

    #[test]
    fn abs_is_real() {
        assert!(close(cval("abs(3 + 4i)"), 5.0, 0.0));
    }

    #[test]
    fn error_paths() {
        assert!(eval_complex_str("1/0 + i").is_err());
        assert!(eval_complex_str("ln(0)").is_err());
        assert!(eval_complex_str("0^(-2)").is_err());
        assert!(eval_complex_str("i + nosuchvar").is_err());
        assert!(eval_complex_str("floor(1 + i)").is_err());
        // parse errors propagate
        assert!(eval_complex_str("1 + * 2").is_err());
    }

    #[test]
    fn eval_complex_api_matches_str() {
        let e = Parser::parse("i^2 + 1").unwrap();
        let ctx = Context::standard();
        let z = eval_complex(&e, &ctx).unwrap();
        assert!(close(z, 0.0, 0.0));
    }
}
