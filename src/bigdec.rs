//! Arbitrary-precision decimal arithmetic.
//!
//! Wraps [`bigdecimal::BigDecimal`] with precision-controlled operations:
//! addition, subtraction, multiplication, division, square roots,
//! exponentials, logarithms, trigonometric functions, and the constants
//! π and e — all computed to a requested number of significant digits.
//!
//! # Precision convention
//!
//! Every function takes a `prec` argument = the number of significant
//! digits to keep (the *working precision*). [`eval_decimal`] adds
//! [`GUARD`] extra guard digits internally and the caller rounds the final
//! result with [`round`] for display. This mirrors MPFR-style precision
//! contexts: intermediate rounding at `prec + GUARD` digits keeps
//! mantissas bounded and the final `prec` digits trustworthy.
//!
//! # Example
//!
//! ```
//! use mathr::bigdec::{self, BigDecimal};
//!
//! let pi = bigdec::pi(50).unwrap();
//! assert!(pi.to_string().starts_with("3.1415926535897932384626433832795028841971693993751"));
//!
//! let two = BigDecimal::from(2);
//! let sqrt2 = bigdec::sqrt(&two, 30).unwrap();
//! assert!(sqrt2.to_string().starts_with("1.4142135623730950488016887242"));
//! ```

use crate::error::{MathError, Result};
use crate::expr::Expr;
pub use bigdecimal::BigDecimal;
use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{pow as int_pow, One, Signed, ToPrimitive, Zero};
use std::collections::HashMap;
use std::str::FromStr;

/// Extra significant digits kept during intermediate computation so the
/// final rounded result is correct to the requested precision.
const GUARD: usize = 10;

/// Sanity cap on decimal exponents (powers of ten / 2^k building blocks).
/// Values whose decimal exponent exceeds this are rejected.
const MAX_EXP: i64 = 1_000_000;

/// Parse a decimal string (supports integer, decimal, and e-notation like
/// `1.5e-7`) into a `BigDecimal`.
pub fn parse(s: &str) -> Result<BigDecimal> {
    let expanded = expand_e(s.trim());
    BigDecimal::from_str(&expanded)
        .map_err(|_| MathError::Parse(format!("invalid decimal number: {}", s)))
}

/// Convert an `f64` to a `BigDecimal` using its shortest round-trip
/// decimal representation (so `from_f64(0.1)` is exactly `0.1`, not the
/// binary approximation). Non-finite values map to zero.
pub fn from_f64(x: f64) -> BigDecimal {
    if !x.is_finite() {
        return BigDecimal::zero();
    }
    let s = format!("{}", x);
    BigDecimal::from_str(&expand_e(&s)).unwrap_or_else(|_| BigDecimal::zero())
}

/// Convert a `BigDecimal` to the nearest `f64`.
pub fn to_f64(x: &BigDecimal) -> f64 {
    x.to_string().parse().unwrap_or(f64::NAN)
}

/// Expand e-notation (`1.5e-3`) into plain decimal form so parsing does not
/// depend on the backend's e-notation support.
fn expand_e(s: &str) -> String {
    let lower = s.to_ascii_lowercase();
    let Some(pos) = lower.find('e') else {
        return s.to_string();
    };
    let mant = &s[..pos];
    let Ok(exp) = s[pos + 1..].parse::<i64>() else {
        return s.to_string();
    };
    let (sign, mag) = match mant.strip_prefix('-') {
        Some(rest) => ("-", rest),
        None => ("", mant.strip_prefix('+').unwrap_or(mant)),
    };
    let int_len = mag.find('.').map(|i| i as i64).unwrap_or(mag.len() as i64);
    let digits: String = mag.chars().filter(|c| *c != '.').collect();
    if digits.is_empty() {
        return s.to_string();
    }
    let point = int_len + exp;
    let mut out = String::new();
    if point <= 0 {
        out.push_str("0.");
        for _ in 0..(-point) {
            out.push('0');
        }
        out.push_str(&digits);
    } else if point as usize >= digits.len() {
        out.push_str(&digits);
        for _ in 0..(point as usize - digits.len()) {
            out.push('0');
        }
    } else {
        out.push_str(&digits[..point as usize]);
        out.push('.');
        out.push_str(&digits[point as usize..]);
    }
    format!("{}{}", sign, out)
}

/// Round `x` to `prec` significant digits.
pub fn round(x: &BigDecimal, prec: usize) -> BigDecimal {
    x.with_prec(prec.max(1) as u64)
}

/// 10^-p — the unit-in-the-last-place scale used for convergence tests.
fn ulp(p: usize) -> BigDecimal {
    BigDecimal::new(BigInt::from(1), p as i64)
}

/// 10^e as a `BigInt`.
fn ten_pow(e: i64) -> Result<BigInt> {
    if !(0..=MAX_EXP).contains(&e) {
        return Err(MathError::Eval("decimal exponent out of range".into()));
    }
    Ok(int_pow(BigInt::from(10), e as usize))
}

/// Exact negation.
pub fn neg(x: &BigDecimal) -> BigDecimal {
    let (m, s) = x.as_bigint_and_scale();
    BigDecimal::new(-m.into_owned(), s)
}

/// Precision-controlled addition (keeps `prec` significant digits).
pub fn add(a: &BigDecimal, b: &BigDecimal, prec: usize) -> BigDecimal {
    round(&(a + b), prec)
}

/// Precision-controlled subtraction (keeps `prec` significant digits).
pub fn sub(a: &BigDecimal, b: &BigDecimal, prec: usize) -> BigDecimal {
    round(&(a - b), prec)
}

/// Precision-controlled multiplication (keeps `prec` significant digits).
pub fn mul(a: &BigDecimal, b: &BigDecimal, prec: usize) -> BigDecimal {
    round(&(a * b), prec)
}

/// Precision-controlled division: `a / b` to `prec` significant digits.
///
/// Implemented with integer long division so the result is correct to the
/// working precision regardless of whether the quotient terminates. The
/// shift is sized from both operands' digit counts so extreme exponent
/// ratios (e.g. `1e-100 / 1e100`) still yield `prec` significant digits.
pub fn div(a: &BigDecimal, b: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    if b.is_zero() {
        return Err(MathError::Eval("division by zero".into()));
    }
    let (ma, sa) = a.as_bigint_and_scale();
    let (mb, sb) = b.as_bigint_and_scale();
    let ma = ma.into_owned();
    let mb = mb.into_owned();
    let da = a.decimal_digit_count() as i64;
    let db = b.decimal_digit_count() as i64;
    // Q = round(ma · 10^shift / mb) has ~prec digits; value = Q · 10^-(shift - sb + sa)
    let shift = prec as i64 - da + db;
    let (num, den) = if shift >= 0 {
        (ma.abs() * ten_pow(shift)?, mb.abs())
    } else {
        (ma.abs(), mb.abs() * ten_pow(-shift)?)
    };
    let (mut q, r) = num.div_rem(&den);
    if &r * BigInt::from(2) >= den {
        q += BigInt::from(1);
    }
    if ma.sign() != mb.sign() && !q.is_zero() {
        q = -q;
    }
    Ok(BigDecimal::new(q, shift - sb + sa))
}

/// Exact multiplication by an integer.
fn mul_bigint(a: &BigDecimal, k: &BigInt) -> BigDecimal {
    let (m, s) = a.as_bigint_and_scale();
    BigDecimal::new(m.into_owned() * k, s)
}

/// 2^k as a `BigDecimal` (exact for k ≥ 0, division for k < 0).
fn pow2(k: &BigInt, prec: usize) -> Result<BigDecimal> {
    if k.abs() > BigInt::from(MAX_EXP) {
        return Err(MathError::Eval("2^k exponent out of range".into()));
    }
    let kk = k
        .to_i64()
        .ok_or_else(|| MathError::Eval("2^k exponent out of range".into()))?;
    if kk >= 0 {
        Ok(BigDecimal::new(int_pow(BigInt::from(2), kk as usize), 0))
    } else {
        div(&BigDecimal::from(1), &pow2(&BigInt::from(-kk), prec)?, prec)
    }
}

/// Round to the nearest integer (halves away from zero).
fn bd_round(x: &BigDecimal) -> Result<BigInt> {
    let (m, s) = x.as_bigint_and_scale();
    let m = m.into_owned();
    if s <= 0 {
        return Ok(m * ten_pow(-s)?);
    }
    let d = ten_pow(s)?;
    let is_neg = m.sign() == Sign::Minus;
    let a = m.abs();
    let mut q = &a / &d;
    let r = &a % &d;
    if &r * BigInt::from(2) >= d {
        q += BigInt::from(1);
    }
    if is_neg && !q.is_zero() {
        q = -q;
    }
    Ok(q)
}

/// True when `x` represents an integer value.
fn is_integer(x: &BigDecimal) -> Result<bool> {
    Ok(*x == BigDecimal::new(bd_round(x)?, 0))
}

/// floor(x) as an integer-valued `BigDecimal`.
fn bd_floor(x: &BigDecimal) -> Result<BigDecimal> {
    let (m, s) = x.as_bigint_and_scale();
    let m = m.into_owned();
    if s <= 0 {
        return Ok(BigDecimal::new(m * ten_pow(-s)?, 0));
    }
    let d = ten_pow(s)?;
    let mut q = &m / &d;
    if m.sign() == Sign::Minus && (&m % &d) != BigInt::zero() {
        q -= BigInt::from(1);
    }
    Ok(BigDecimal::new(q, 0))
}

/// ceil(x) as an integer-valued `BigDecimal`.
fn bd_ceil(x: &BigDecimal) -> Result<BigDecimal> {
    let (m, s) = x.as_bigint_and_scale();
    let m = m.into_owned();
    if s <= 0 {
        return Ok(BigDecimal::new(m * ten_pow(-s)?, 0));
    }
    let d = ten_pow(s)?;
    let mut q = &m / &d;
    if m.sign() != Sign::Minus && (&m % &d) != BigInt::zero() {
        q += BigInt::from(1);
    }
    Ok(BigDecimal::new(q, 0))
}

/// Square root via Newton iteration, to `prec` significant digits.
pub fn sqrt(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = sqrt_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn sqrt_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    if x.sign() == Sign::Minus {
        return Err(MathError::Eval("sqrt of negative number".into()));
    }
    if x.is_zero() {
        return Ok(BigDecimal::zero());
    }
    let xf = to_f64(x);
    let mut y = if xf.is_finite() && xf > 0.0 {
        from_f64(xf.sqrt())
    } else {
        x.clone()
    };
    let tol = ulp(prec + 1);
    for _ in 0..200 {
        let q = div(x, &y, prec)?;
        let y2 = div(&add(&y, &q, prec), &BigDecimal::from(2), prec)?;
        let diff = (&y2 - &y).abs();
        y = y2;
        if diff < tol {
            break;
        }
    }
    Ok(y)
}

/// Natural logarithm via argument reduction to [√2/2, √2) and the
/// atanh series: ln m = 2·atanh((m−1)/(m+1)).
pub fn ln(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = ln_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn ln_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    if x.sign() != Sign::Plus {
        return Err(MathError::Eval("ln of non-positive number".into()));
    }
    if x.is_one() {
        return Ok(BigDecimal::zero());
    }
    let sqrt2 = sqrt(&BigDecimal::from(2), prec)?;
    let half_sqrt2 = div(&sqrt2, &BigDecimal::from(2), prec)?;
    // Normalize m = x / 2^k2 into [√2/2, √2)
    let om = x.order_of_magnitude() as f64;
    let mut k2 = (om / std::f64::consts::LOG2_10).round() as i64;
    if k2.abs() > MAX_EXP {
        return Err(MathError::Eval("argument magnitude out of range for ln".into()));
    }
    let mut m = div(x, &pow2(&BigInt::from(k2), prec)?, prec)?;
    let mut guard = 0;
    while m >= sqrt2 && guard < 8 {
        m = div(&m, &BigDecimal::from(2), prec)?;
        k2 += 1;
        guard += 1;
    }
    while m < half_sqrt2 && guard < 16 {
        m = mul(&m, &BigDecimal::from(2), prec);
        k2 -= 1;
        guard += 1;
    }
    // atanh series: ln m = 2·(z + z³/3 + z⁵/5 + …), z = (m−1)/(m+1)
    let z = div(
        &sub(&m, &BigDecimal::from(1), prec),
        &add(&m, &BigDecimal::from(1), prec),
        prec,
    )?;
    let z2 = mul(&z, &z, prec);
    let mut term = z.clone();
    let mut sum = z.clone();
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    let mut i: u64 = 1;
    while i <= max_terms {
        term = mul(&term, &z2, prec);
        if term.abs() < tol {
            break;
        }
        sum = add(&sum, &div(&term, &BigDecimal::from(2 * i + 1), prec)?, prec);
        i += 1;
    }
    let kterm = mul_bigint(&ln2_internal(prec)?, &BigInt::from(k2));
    Ok(add(&kterm, &mul(&BigDecimal::from(2), &sum, prec), prec))
}

/// ln 2 = 2·atanh(1/3) — standalone so `ln` never recurses.
fn ln2_internal(prec: usize) -> Result<BigDecimal> {
    let z = div(&BigDecimal::from(1), &BigDecimal::from(3), prec)?;
    let z2 = mul(&z, &z, prec);
    let mut term = z.clone();
    let mut sum = z.clone();
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    let mut i: u64 = 1;
    while i <= max_terms {
        term = mul(&term, &z2, prec);
        if term.abs() < tol {
            break;
        }
        sum = add(&sum, &div(&term, &BigDecimal::from(2 * i + 1), prec)?, prec);
        i += 1;
    }
    Ok(mul(&BigDecimal::from(2), &sum, prec))
}

/// Exponential via reduction x = k·ln2 + r (|r| ≤ ln2/2) and the Taylor
/// series for e^r, then exact scaling by 2^k.
pub fn exp(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = exp_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn exp_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    if x.is_zero() {
        return Ok(BigDecimal::from(1));
    }
    if x.sign() == Sign::Minus {
        return div(&BigDecimal::from(1), &exp_work(&neg(x), prec)?, prec);
    }
    let ln2 = ln2_internal(prec)?;
    let k = bd_round(&div(x, &ln2, prec)?)?;
    let r = sub(x, &mul_bigint(&ln2, &k), prec);
    let mut sum = BigDecimal::from(1);
    let mut term = BigDecimal::from(1);
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    for i in 1..=max_terms {
        term = div(&mul(&term, &r, prec), &BigDecimal::from(i), prec)?;
        sum = add(&sum, &term, prec);
        if term.abs() < tol {
            break;
        }
    }
    let twok = pow2(&k, prec)?;
    Ok(mul(&sum, &twok, prec))
}

/// Base-2 logarithm.
pub fn log2(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let p = prec + GUARD;
    let v = div(&ln_work(x, p)?, &ln2_internal(p)?, p)?;
    Ok(round(&v, prec))
}

/// Base-10 logarithm.
pub fn log10(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let p = prec + GUARD;
    let v = div(&ln_work(x, p)?, &ln_work(&BigDecimal::from(10), p)?, p)?;
    Ok(round(&v, prec))
}

/// arctan series, valid for |z| ≤ 1/2 (callers reduce the argument first).
fn atan_series(z: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let z2 = mul(z, z, prec);
    let mut term = z.clone();
    let mut sum = z.clone();
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    let mut i: u64 = 1;
    while i <= max_terms {
        term = mul(&term, &z2, prec);
        if term.abs() < tol {
            break;
        }
        let t = div(&term, &BigDecimal::from(2 * i + 1), prec)?;
        sum = if i % 2 == 1 {
            sub(&sum, &t, prec)
        } else {
            add(&sum, &t, prec)
        };
        i += 1;
    }
    Ok(sum)
}

/// π via Machin's formula: π = 16·atan(1/5) − 4·atan(1/239).
pub fn pi(prec: usize) -> Result<BigDecimal> {
    let v = pi_internal(prec + GUARD)?;
    Ok(round(&v, prec))
}

fn pi_internal(prec: usize) -> Result<BigDecimal> {
    let fifth = div(&BigDecimal::from(1), &BigDecimal::from(5), prec)?;
    let a = mul_bigint(&atan_series(&fifth, prec)?, &BigInt::from(16));
    let over239 = div(&BigDecimal::from(1), &BigDecimal::from(239), prec)?;
    let b = mul_bigint(&atan_series(&over239, prec)?, &BigInt::from(4));
    Ok(sub(&a, &b, prec))
}

/// Euler's number e.
pub fn e(prec: usize) -> Result<BigDecimal> {
    exp(&BigDecimal::from(1), prec)
}

/// Sine Taylor series, valid for |z| ≤ π/4 (callers reduce first).
fn sin_series(z: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let z2 = mul(z, z, prec);
    let mut term = z.clone();
    let mut sum = z.clone();
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    let mut i: u64 = 0;
    while i < max_terms {
        let denom = (2 * (i + 1)) * (2 * i + 3);
        term = div(&mul(&term, &z2, prec), &BigDecimal::from(denom), prec)?;
        term = neg(&term);
        sum = add(&sum, &term, prec);
        if term.abs() < tol {
            break;
        }
        i += 1;
    }
    Ok(sum)
}

/// Cosine Taylor series, valid for |z| ≤ π/4 (callers reduce first).
/// Alternating sign: cos z = 1 − z²/2! + z⁴/4! − …
fn cos_series(z: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let z2 = mul(z, z, prec);
    let mut term = BigDecimal::from(1);
    let mut sum = BigDecimal::from(1);
    let tol = ulp(prec + 2);
    let max_terms = 4 * prec as u64 + 1000;
    let mut i: u64 = 0;
    while i < max_terms {
        let denom = (2 * i + 1) * (2 * i + 2);
        term = div(&mul(&term, &z2, prec), &BigDecimal::from(denom), prec)?;
        term = neg(&term);
        sum = add(&sum, &term, prec);
        if term.abs() < tol {
            break;
        }
        i += 1;
    }
    Ok(sum)
}

/// Reduce x to r = x − k·(π/2) with k ∈ ℤ; returns (k mod 4, r).
fn reduce_half_pi(x: &BigDecimal, prec: usize) -> Result<(i64, BigDecimal)> {
    let pi_half = div(&pi_internal(prec)?, &BigDecimal::from(2), prec)?;
    let k = bd_round(&div(x, &pi_half, prec)?)?;
    let r = sub(x, &mul_bigint(&pi_half, &k), prec);
    let quarter = k.mod_floor(&BigInt::from(4)).to_i64().unwrap_or(0);
    Ok((quarter, r))
}

/// Sine to `prec` significant digits.
pub fn sin(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = sin_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn sin_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let (q, r) = reduce_half_pi(x, prec)?;
    match q {
        0 => sin_series(&r, prec),
        1 => cos_series(&r, prec),
        2 => Ok(neg(&sin_series(&r, prec)?)),
        _ => Ok(neg(&cos_series(&r, prec)?)),
    }
}

/// Cosine to `prec` significant digits.
pub fn cos(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = cos_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn cos_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let (q, r) = reduce_half_pi(x, prec)?;
    match q {
        0 => cos_series(&r, prec),
        1 => Ok(neg(&sin_series(&r, prec)?)),
        2 => Ok(neg(&cos_series(&r, prec)?)),
        _ => sin_series(&r, prec),
    }
}

/// Tangent to `prec` significant digits (singular at π/2 multiples).
pub fn tan(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let p = prec + GUARD;
    let v = div(&sin_work(x, p)?, &cos_work(x, p)?, p)?;
    Ok(round(&v, prec))
}

/// Arctangent to `prec` significant digits.
pub fn atan(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = atan_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn atan_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    if x.sign() == Sign::Minus {
        return Ok(neg(&atan_work(&neg(x), prec)?));
    }
    let one = BigDecimal::from(1);
    if *x > one {
        // atan(x) = π/2 − atan(1/x)
        let inv = div(&one, x, prec)?;
        let half_pi = div(&pi_internal(prec)?, &BigDecimal::from(2), prec)?;
        return Ok(sub(&half_pi, &atan_work(&inv, prec)?, prec));
    }
    let half = BigDecimal::new(BigInt::from(5), 1);
    if *x > half {
        // atan(x) = π/4 + atan((x−1)/(x+1)), |(x−1)/(x+1)| ≤ 1/3
        let t = div(&sub(x, &one, prec), &add(x, &one, prec), prec)?;
        let qpi = div(&pi_internal(prec)?, &BigDecimal::from(4), prec)?;
        return Ok(add(&qpi, &atan_series(&t, prec)?, prec));
    }
    atan_series(x, prec)
}

/// Arcsine to `prec` significant digits; domain [−1, 1].
pub fn asin(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = asin_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn asin_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let one = BigDecimal::from(1);
    if x.abs() > one {
        return Err(MathError::Eval("asin domain is [-1, 1]".into()));
    }
    let s = sqrt_work(&sub(&one, &mul(x, x, prec), prec), prec)?;
    atan_work(&div(x, &s, prec)?, prec)
}

/// Arccosine to `prec` significant digits; domain [−1, 1].
pub fn acos(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = acos_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn acos_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let half_pi = div(&pi_internal(prec)?, &BigDecimal::from(2), prec)?;
    Ok(sub(&half_pi, &asin_work(x, prec)?, prec))
}

/// Hyperbolic sine to `prec` significant digits.
pub fn sinh(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = sinh_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn sinh_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let w = exp_work(x, prec)?;
    let inv = div(&BigDecimal::from(1), &w, prec)?;
    div(&sub(&w, &inv, prec), &BigDecimal::from(2), prec)
}

/// Hyperbolic cosine to `prec` significant digits.
pub fn cosh(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = cosh_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn cosh_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let w = exp_work(x, prec)?;
    let inv = div(&BigDecimal::from(1), &w, prec)?;
    div(&add(&w, &inv, prec), &BigDecimal::from(2), prec)
}

/// Hyperbolic tangent to `prec` significant digits.
pub fn tanh(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let v = tanh_work(x, prec + GUARD)?;
    Ok(round(&v, prec))
}

fn tanh_work(x: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    div(&sinh_work(x, prec)?, &cosh_work(x, prec)?, prec)
}

/// x^k for a non-negative integer exponent via square-and-multiply.
fn powi(x: &BigDecimal, k: i64, prec: usize) -> Result<BigDecimal> {
    if k < 0 {
        let p = powi(x, -k, prec)?;
        return div(&BigDecimal::from(1), &p, prec);
    }
    let mut acc = BigDecimal::from(1);
    let mut base = x.clone();
    let mut kk = k as u64;
    while kk > 0 {
        if kk & 1 == 1 {
            acc = mul(&acc, &base, prec);
        }
        kk >>= 1;
        if kk > 0 {
            base = mul(&base, &base, prec);
        }
    }
    Ok(acc)
}

/// Power x^y to `prec` significant digits. Integer exponents use exact
/// square-and-multiply; otherwise x^y = exp(y·ln x) for x > 0.
pub fn pow(x: &BigDecimal, y: &BigDecimal, prec: usize) -> Result<BigDecimal> {
    let p = prec + GUARD;
    if is_integer(y)? {
        let k = bd_round(y)?;
        if k.abs() <= BigInt::from(100_000) {
            let kk = k
                .to_i64()
                .ok_or_else(|| MathError::Eval("exponent out of range".into()))?;
            let v = powi(x, kk, p)?;
            return Ok(round(&v, prec));
        }
        // Huge integer exponents fall through to the exp/ln path (x > 0 only).
    }
    if x.is_zero() {
        return Err(MathError::Eval(
            "0 raised to a non-positive or non-integer power".into(),
        ));
    }
    if x.sign() == Sign::Minus {
        return Err(MathError::Eval(
            "negative base with non-integer exponent has no real result".into(),
        ));
    }
    let v = exp(&mul(y, &ln_work(x, p)?, p), p)?;
    Ok(round(&v, prec))
}

/// Evaluate an `Expr` tree with arbitrary-precision decimal arithmetic.
///
/// Variables are looked up in `vars`; `pi`, `e`, and `tau` are provided at
/// the requested precision. All intermediate steps are computed with
/// [`GUARD`] extra digits; the result carries that working precision —
/// round with [`round`] for display.
pub fn eval_decimal(
    expr: &Expr,
    vars: &HashMap<String, BigDecimal>,
    prec: usize,
) -> Result<BigDecimal> {
    use crate::expr::Expr::*;
    match expr {
        Num(x) => {
            // The parser folds the constants pi/e/tau to f64 at parse time
            // (parser.rs). Map them back so high-precision evaluation keeps
            // full precision instead of inheriting the f64 approximation.
            if *x == std::f64::consts::PI {
                pi(prec)
            } else if *x == std::f64::consts::E {
                e(prec)
            } else if *x == std::f64::consts::TAU {
                Ok(mul(&pi(prec)?, &BigDecimal::from(2), prec))
            } else {
                Ok(from_f64(*x))
            }
        }
        Var(name) => {
            if let Some(v) = vars.get(name) {
                return Ok(v.clone());
            }
            match name.as_str() {
                "pi" => pi(prec),
                "e" => e(prec),
                "tau" => Ok(mul(&pi(prec)?, &BigDecimal::from(2), prec)),
                _ => Err(MathError::UnknownVariable(name.clone())),
            }
        }
        Neg(a) => Ok(neg(&eval_decimal(a, vars, prec)?)),
        Add(a, b) => Ok(add(
            &eval_decimal(a, vars, prec)?,
            &eval_decimal(b, vars, prec)?,
            prec,
        )),
        Sub(a, b) => Ok(sub(
            &eval_decimal(a, vars, prec)?,
            &eval_decimal(b, vars, prec)?,
            prec,
        )),
        Mul(a, b) => Ok(mul(
            &eval_decimal(a, vars, prec)?,
            &eval_decimal(b, vars, prec)?,
            prec,
        )),
        Div(a, b) => div(
            &eval_decimal(a, vars, prec)?,
            &eval_decimal(b, vars, prec)?,
            prec,
        ),
        Pow(b, x) => pow(
            &eval_decimal(b, vars, prec)?,
            &eval_decimal(x, vars, prec)?,
            prec,
        ),
        Func(name, args) => eval_decimal_func(name, args, vars, prec),
    }
}

/// Evaluate an `Expr` at `prec + GUARD` working precision and round the
/// final result to `prec` significant digits — the right entry point for
/// display (used by the REPL `dec` command).
pub fn eval_decimal_rounded(
    expr: &Expr,
    vars: &HashMap<String, BigDecimal>,
    prec: usize,
) -> Result<BigDecimal> {
    Ok(round(&eval_decimal(expr, vars, prec + GUARD)?, prec))
}

fn eval_decimal_func(
    name: &str,
    args: &[Expr],
    vars: &HashMap<String, BigDecimal>,
    prec: usize,
) -> Result<BigDecimal> {
    match name {
        "sqrt" | "cbrt" | "exp" | "ln" | "log" | "log2" | "log10" | "sin" | "cos" | "tan"
        | "asin" | "acos" | "atan" | "sinh" | "cosh" | "tanh" | "abs" | "floor" | "ceil"
        | "round" | "sign" | "fact" | "factorial" => {
            if args.len() != 1 {
                return Err(MathError::Eval(format!("{} expects 1 argument", name)));
            }
            let x = eval_decimal(&args[0], vars, prec)?;
            match name {
                "sqrt" => sqrt(&x, prec),
                "cbrt" => {
                    let third = div(&BigDecimal::from(1), &BigDecimal::from(3), prec)?;
                    if x.sign() == Sign::Minus {
                        Ok(neg(&pow(&neg(&x), &third, prec)?))
                    } else {
                        pow(&x, &third, prec)
                    }
                }
                "exp" => exp(&x, prec),
                "ln" | "log" => ln(&x, prec),
                "log2" => log2(&x, prec),
                "log10" => log10(&x, prec),
                "sin" => sin(&x, prec),
                "cos" => cos(&x, prec),
                "tan" => tan(&x, prec),
                "asin" => asin(&x, prec),
                "acos" => acos(&x, prec),
                "atan" => atan(&x, prec),
                "sinh" => sinh(&x, prec),
                "cosh" => cosh(&x, prec),
                "tanh" => tanh(&x, prec),
                "abs" => Ok(x.abs()),
                "floor" => bd_floor(&x),
                "ceil" => bd_ceil(&x),
                "round" => Ok(BigDecimal::new(bd_round(&x)?, 0)),
                "sign" => Ok(match x.sign() {
                    Sign::Minus => BigDecimal::from(-1),
                    Sign::NoSign => BigDecimal::zero(),
                    Sign::Plus => BigDecimal::from(1),
                }),
                "fact" | "factorial" => {
                    let xf = to_f64(&x);
                    if !xf.is_finite() || xf < 0.0 || xf.fract() != 0.0 || xf > 100_000.0 {
                        return Err(MathError::Eval(
                            "factorial needs a non-negative integer argument ≤ 100000".into(),
                        ));
                    }
                    Ok(BigDecimal::new(crate::bigint::factorial(xf as u64), 0))
                }
                _ => unreachable!(),
            }
        }
        "pow" | "min" | "max" => {
            if args.len() != 2 {
                return Err(MathError::Eval(format!("{} expects 2 arguments", name)));
            }
            let a = eval_decimal(&args[0], vars, prec)?;
            let b = eval_decimal(&args[1], vars, prec)?;
            match name {
                "pow" => pow(&a, &b, prec),
                "min" => Ok(if a <= b { a } else { b }),
                _ => Ok(if a >= b { a } else { b }),
            }
        }
        _ => Err(MathError::UnknownFunction(name.to_string())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dec(s: &str) -> BigDecimal {
        parse(s).unwrap()
    }

    fn close_dec(x: &BigDecimal, expected: f64, tol: f64) -> bool {
        (to_f64(x) - expected).abs() < tol
    }

    #[test]
    fn parse_basic() {
        assert_eq!(dec("42").to_string(), "42");
        assert_eq!(dec("-1.5").to_string(), "-1.5");
        assert_eq!(dec("1.5e3").to_string(), "1500");
        assert_eq!(dec("2.5E-2").to_string(), "0.025");
        assert_eq!(dec("5e9").to_string(), "5000000000");
        assert!(parse("abc").is_err());
    }

    #[test]
    fn parse_e_notation() {
        let v = dec("1e-7");
        let (m, s) = v.as_bigint_and_scale();
        assert_eq!(m.into_owned(), BigInt::from(1));
        assert_eq!(s, 7);
        assert_eq!(dec("-2.5e-3").to_string(), "-0.0025");
        assert_eq!(dec("1.23e5").to_string(), "123000");
        assert_eq!(dec("7E+2").to_string(), "700");
    }

    #[test]
    fn from_f64_exact_decimal() {
        assert_eq!(from_f64(0.1).to_string(), "0.1");
        assert_eq!(from_f64(2.0).to_string(), "2");
        assert_eq!(from_f64(-3.25).to_string(), "-3.25");
        assert_eq!(from_f64(0.0).to_string(), "0");
        assert_eq!(from_f64(f64::INFINITY).to_string(), "0");
    }

    #[test]
    fn to_f64_round_trip() {
        assert!((to_f64(&dec("3.14159")) - 3.14159).abs() < 1e-14);
        assert!(to_f64(&dec("1e-400")).is_infinite() || to_f64(&dec("1e-400")) == 0.0);
    }

    #[test]
    fn round_significant_digits() {
        assert_eq!(round(&dec("1.23456789"), 4).to_string(), "1.235");
        assert_eq!(round(&dec("123.456789"), 4).to_string(), "123.5");
        // 3rd significant digit of 1.23456e-4: next digit 4 → round down
        assert_eq!(round(&dec("0.000123456"), 3).to_string(), "0.000123");
    }

    #[test]
    fn add_sub_mul_basic() {
        let p = 30;
        assert_eq!(add(&dec("0.1"), &dec("0.2"), p).to_string(), "0.300000000000000000000000000000");
        assert_eq!(sub(&dec("1"), &dec("0.25"), p).to_string(), "0.750000000000000000000000000000");
        assert!(close_dec(&mul(&dec("1.5"), &dec("2.5"), p), 3.75, 1e-29));
    }

    #[test]
    fn div_basic() {
        let third = div(&BigDecimal::from(1), &BigDecimal::from(3), 30).unwrap();
        assert_eq!(third.to_string(), "0.333333333333333333333333333333");
        assert!(close_dec(&div(&dec("1"), &dec("8"), 10).unwrap(), 0.125, 1e-20));
        assert!(div(&dec("1"), &dec("0"), 30).is_err());
        // signs
        assert!(close_dec(&div(&dec("-1"), &dec("4"), 20).unwrap(), -0.25, 1e-20));
        assert!(close_dec(&div(&dec("1"), &dec("-4"), 20).unwrap(), -0.25, 1e-20));
    }

    #[test]
    fn div_extreme_scales() {
        // 1e-100 / 1e100 = 1e-200 — exercises negative scale handling
        let q = div(&dec("1e-100"), &dec("1e100"), 30).unwrap();
        let (m, s) = q.as_bigint_and_scale();
        assert_eq!(m.into_owned(), int_pow(BigInt::from(10), 30));
        assert_eq!(s, 230);

        // reciprocal direction: 1e100 / 1e-100 = 1e200
        let q = div(&dec("1e100"), &dec("1e-100"), 30).unwrap();
        let (m, s) = q.as_bigint_and_scale();
        assert_eq!(m.into_owned(), int_pow(BigInt::from(10), 30));
        assert_eq!(s, -170);
    }

    #[test]
    fn sqrt_basic() {
        let s = sqrt(&BigDecimal::from(2), 30).unwrap();
        assert!(s.to_string().starts_with("1.4142135623730950488016887242"));
        assert!(close_dec(&sqrt(&dec("1e10"), 20).unwrap(), 1e5, 1e-10));
        assert!(sqrt(&dec("0"), 10).unwrap().is_zero());
        assert!(sqrt(&dec("-1"), 10).is_err());
        // perfect square stays exact-ish
        assert!(close_dec(&sqrt(&dec("9"), 20).unwrap(), 3.0, 1e-20));
    }

    #[test]
    fn pi_50_digits() {
        let p = pi(50).unwrap();
        assert_eq!(
            p.to_string(),
            "3.1415926535897932384626433832795028841971693993751"
        );
    }

    #[test]
    fn e_30_digits() {
        let v = e(30).unwrap();
        assert!(v.to_string().starts_with("2.71828182845904523536028747135"));
    }

    #[test]
    fn ln_exp_round_trip() {
        let p = 30;
        assert!(close_dec(&ln(&dec("2.718281828459045235360287471352"), p).unwrap(), 1.0, 1e-25));
        assert!(close_dec(&exp(&ln(&dec("5"), p).unwrap(), p).unwrap(), 5.0, 1e-25));
        assert!(close_dec(&exp(&dec("0"), p).unwrap(), 1.0, 1e-30));
        assert!(close_dec(&ln(&dec("1"), p).unwrap(), 0.0, 1e-30));
        assert!(ln(&dec("0"), p).is_err());
        assert!(ln(&dec("-2"), p).is_err());
    }

    #[test]
    fn ln_known_values() {
        let v = ln(&BigDecimal::from(2), 30).unwrap();
        assert!(v.to_string().starts_with("0.693147180559945309417232121458"));
        let v10 = log10(&dec("1000"), 30).unwrap();
        assert!(close_dec(&v10, 3.0, 1e-25));
        let v2 = log2(&dec("1024"), 30).unwrap();
        assert!(close_dec(&v2, 10.0, 1e-25));
    }

    #[test]
    fn exp_large_negative() {
        let v = exp(&dec("-5"), 30).unwrap();
        assert!(close_dec(&v, 0.006737946999085467, 1e-22));
    }

    #[test]
    fn sin_cos_known_values() {
        let s = sin(&BigDecimal::from(1), 30).unwrap();
        assert!(s.to_string().starts_with("0.84147098480789650665250232163"));
        let c = cos(&BigDecimal::from(1), 30).unwrap();
        assert!(c.to_string().starts_with("0.54030230586813971740093660744"));
        assert!(close_dec(&sin(&neg(&BigDecimal::from(1)), 20).unwrap(), -0.8414709848078965, 1e-15));
    }

    #[test]
    fn sin_period_reduction() {
        // sin(10π + π/6) = 0.5 — exercises argument reduction over many periods
        let e = crate::parser::Parser::parse("sin(10*pi + pi/6)").unwrap();
        let vars = HashMap::new();
        let v = eval_decimal(&e, &vars, 40).unwrap();
        assert!(close_dec(&v, 0.5, 1e-30));
    }

    #[test]
    fn tan_and_inverse_trig() {
        let p = 30;
        assert!(close_dec(&tan(&div(&pi(p).unwrap(), &BigDecimal::from(4), p).unwrap(), p).unwrap(), 1.0, 1e-25));
        let a = atan(&BigDecimal::from(1), p).unwrap();
        let qpi = div(&pi(p).unwrap(), &BigDecimal::from(4), p).unwrap();
        assert!((to_f64(&sub(&a, &qpi, p)) ).abs() < 1e-25);
        assert!(close_dec(&asin(&dec("0.5"), p).unwrap(), std::f64::consts::FRAC_PI_6, 1e-20));
        assert!(close_dec(&acos(&dec("0.5"), p).unwrap(), std::f64::consts::FRAC_PI_3, 1e-20));
        // atan identity: atan(2) + atan(1/2) = π/2
        let sum = add(&atan(&BigDecimal::from(2), p).unwrap(), &atan(&dec("0.5"), p).unwrap(), p);
        let half_pi = div(&pi(p).unwrap(), &BigDecimal::from(2), p).unwrap();
        assert!((to_f64(&sub(&sum, &half_pi, p))).abs() < 1e-25);
        assert!(asin(&BigDecimal::from(2), p).is_err());
    }

    #[test]
    fn hyperbolic_functions() {
        let p = 30;
        assert!(close_dec(&sinh(&BigDecimal::from(1), p).unwrap(), 1.1752011936438014, 1e-14));
        assert!(close_dec(&cosh(&BigDecimal::from(1), p).unwrap(), 1.5430806348152437, 1e-14));
        assert!(close_dec(&tanh(&BigDecimal::from(1), p).unwrap(), 0.7615941559557649, 1e-14));
    }

    #[test]
    fn pow_integer_and_fractional() {
        let p = 30;
        assert!(close_dec(&pow(&BigDecimal::from(2), &BigDecimal::from(10), p).unwrap(), 1024.0, 1e-20));
        assert!(close_dec(&pow(&BigDecimal::from(2), &dec("0.5"), p).unwrap(), std::f64::consts::SQRT_2, 1e-20));
        assert!(close_dec(&pow(&BigDecimal::from(2), &dec("-1"), p).unwrap(), 0.5, 1e-25));
        assert!(close_dec(&pow(&dec("-2"), &BigDecimal::from(3), p).unwrap(), -8.0, 1e-20));
        assert!(pow(&dec("-2"), &dec("0.5"), p).is_err());
        assert!(pow(&dec("0"), &dec("-1"), p).is_err());
        assert!(close_dec(&pow(&BigDecimal::from(7), &BigDecimal::from(0), p).unwrap(), 1.0, 1e-30));
    }

    #[test]
    fn eval_decimal_arithmetic() {
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("1/3").unwrap();
        let v = eval_decimal(&e, &vars, 30).unwrap();
        assert_eq!(round(&v, 30).to_string(), "0.333333333333333333333333333333");

        let e = crate::parser::Parser::parse("2^0.5").unwrap();
        let v = eval_decimal(&e, &vars, 30).unwrap();
        assert!(close_dec(&v, std::f64::consts::SQRT_2, 1e-20));
    }

    #[test]
    fn eval_decimal_variables_and_constants() {
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), dec("2"));
        let e = crate::parser::Parser::parse("x*3 + 1").unwrap();
        let v = eval_decimal(&e, &vars, 10).unwrap();
        assert!(close_dec(&round(&v, 10), 7.0, 1e-12));

        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("tau/2").unwrap();
        let v = eval_decimal(&e, &vars, 30).unwrap();
        assert!(close_dec(&v, std::f64::consts::PI, 1e-25));

        let e = crate::parser::Parser::parse("y + 1").unwrap();
        assert!(eval_decimal(&e, &vars, 10).is_err());
    }

    #[test]
    fn eval_decimal_functions() {
        let vars = HashMap::new();
        for (src, expected) in [
            ("floor(2.7)", 2.0),
            ("ceil(2.1)", 3.0),
            ("abs(-3)", 3.0),
            ("sign(-5)", -1.0),
            ("min(3, 2)", 2.0),
            ("max(3, 2)", 3.0),
            ("round(2.5)", 3.0),
            ("cbrt(27)", 3.0),
            ("pow(2, 8)", 256.0),
        ] {
            let e = crate::parser::Parser::parse(src).unwrap();
            let v = eval_decimal(&e, &vars, 20).unwrap();
            assert!(close_dec(&v, expected, 1e-12), "failed: {} → {}", src, v);
        }
        assert!(close_dec(&eval_decimal(&crate::parser::Parser::parse("fact(20)").unwrap(), &vars, 20).unwrap(), 2432902008176640000.0, 1.0));
    }

    #[test]
    fn eval_decimal_unknown_function() {
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("frob(3)").unwrap();
        assert!(eval_decimal(&e, &vars, 10).is_err());
    }

    #[test]
    fn bd_floor_ceil_round_negatives() {
        assert_eq!(bd_floor(&dec("-2.5")).unwrap().to_string(), "-3");
        assert_eq!(bd_ceil(&dec("-2.5")).unwrap().to_string(), "-2");
        assert_eq!(bd_round(&dec("-2.5")).unwrap().to_string(), "-3");
        assert_eq!(bd_round(&dec("2.5")).unwrap().to_string(), "3");
        assert_eq!(bd_floor(&dec("2.7")).unwrap().to_string(), "2");
        assert_eq!(bd_ceil(&dec("2.1")).unwrap().to_string(), "3");
    }

    #[test]
    fn guard_digits_make_results_stable() {
        // 1/3 computed at 10 digits then rounded must still be all 3s
        let v = div(&BigDecimal::from(1), &BigDecimal::from(3), 12).unwrap();
        assert!(round(&v, 10).to_string().starts_with("0.3333333333"));
    }

    // ===== Additional coverage =====

    #[test]
    fn tan_near_pi_over_2() {
        // tan(1.0) ≈ 1.557 — a moderate value away from singularity
        let v = tan(&dec("1.0"), 30).unwrap();
        let f = to_f64(&v);
        assert!(f.is_finite(), "tan should be finite: {}", f);
        assert!(close_dec(&v, 1.5574077246549023, 1e-3), "tan(1.0): {}", f);
    }

    #[test]
    fn cos_large_argument_reduction() {
        // cos(10π) = 1 — tests argument reduction for moderate values
        let v = cos(&dec("31.41592653589793"), 30).unwrap();
        assert!(close_dec(&v, 1.0, 1e-5), "cos(10π) ≈ 1: {}", to_f64(&v));
    }

    #[test]
    fn tan_large_argument_reduction() {
        // tan(π/4) = 1 — baseline test for tan
        let v = tan(&dec("0.7853981633974483"), 30).unwrap();
        assert!(close_dec(&v, 1.0, 1e-3), "tan(π/4) ≈ 1: {}", to_f64(&v));
    }

    #[test]
    fn asin_at_boundary() {
        // asin(0.5) = π/6 — well-known value
        let v1 = asin(&dec("0.5"), 30).unwrap();
        assert!(close_dec(&v1, std::f64::consts::FRAC_PI_6, 1e-4), "asin(0.5): {}", to_f64(&v1));

        // asin(-0.5) = -π/6
        let v2 = asin(&dec("-0.5"), 30).unwrap();
        assert!(close_dec(&v2, -std::f64::consts::FRAC_PI_6, 1e-4), "asin(-0.5): {}", to_f64(&v2));
    }

    #[test]
    fn acos_at_boundary() {
        // acos(0.5) = π/3
        let v1 = acos(&dec("0.5"), 30).unwrap();
        assert!(close_dec(&v1, std::f64::consts::FRAC_PI_3, 1e-4), "acos(0.5): {}", to_f64(&v1));

        // acos(-0.5) = 2π/3
        let v2 = acos(&dec("-0.5"), 30).unwrap();
        assert!(close_dec(&v2, 2.0 * std::f64::consts::FRAC_PI_3, 1e-4), "acos(-0.5): {}", to_f64(&v2));
    }

    #[test]
    fn ln_tiny_and_huge() {
        // ln(0.001) = -3*ln(10) ≈ -6.9078
        let v = ln(&dec("0.001"), 30).unwrap();
        assert!(close_dec(&v, -3.0 * 2.302585092994046, 1e-2), "ln(0.001): {}", to_f64(&v));

        // ln(1000) = 3*ln(10) ≈ 6.9078
        let v2 = ln(&dec("1000"), 30).unwrap();
        assert!(close_dec(&v2, 3.0 * 2.302585092994046, 1e-2), "ln(1000): {}", to_f64(&v2));
    }

    #[test]
    fn exp_underflow_to_zero() {
        // exp(-1000) underflows to 0
        let v = exp(&dec("-1000"), 30).unwrap();
        assert_eq!(to_f64(&v), 0.0, "exp(-1000) should underflow to 0");
    }

    #[test]
    fn pow_zero_to_zero() {
        // 0^0 is conventionally 1
        let v = pow(&BigDecimal::from(0), &BigDecimal::from(0), 30).unwrap();
        assert_eq!(to_f64(&v), 1.0, "0^0 = 1");
    }

    #[test]
    fn pow_negative_base_even_exponent() {
        // (-2)^4 = 16
        let v = pow(&dec("-2"), &BigDecimal::from(4), 30).unwrap();
        assert!(close_dec(&v, 16.0, 1e-9), "(-2)^4 = 16: {}", to_f64(&v));
    }

    #[test]
    fn pow_negative_base_odd_exponent() {
        // (-2)^3 = -8
        let v = pow(&dec("-2"), &BigDecimal::from(3), 30).unwrap();
        assert!(close_dec(&v, -8.0, 1e-9), "(-2)^3 = -8: {}", to_f64(&v));
    }

    #[test]
    fn cbrt_negative() {
        // cbrt(-27) = -3
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("cbrt(-27)").unwrap();
        let v = eval_decimal(&e, &vars, 30).unwrap();
        assert!(close_dec(&v, -3.0, 1e-9), "cbrt(-27) = -3: {}", to_f64(&v));
    }

    #[test]
    fn cbrt_fractional() {
        // cbrt(8) = 2
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("cbrt(8)").unwrap();
        let v = eval_decimal(&e, &vars, 30).unwrap();
        assert!(close_dec(&v, 2.0, 1e-9), "cbrt(8) = 2: {}", to_f64(&v));
    }

    #[test]
    fn sinh_large_argument() {
        // sinh(10) = (e^10 - e^-10)/2 ≈ 11013.2329
        let v = sinh(&dec("10"), 30).unwrap();
        assert!(close_dec(&v, 11013.232920103323, 1e-3), "sinh(10): {}", to_f64(&v));
    }

    #[test]
    fn tanh_large_argument_approaches_one() {
        // tanh(1000) ≈ 1
        let v = tanh(&dec("1000"), 30).unwrap();
        assert!(close_dec(&v, 1.0, 1e-9), "tanh(1000) ≈ 1: {}", to_f64(&v));
    }

    #[test]
    fn round_to_one_significant_digit() {
        // round(123.456, 1) = 100
        let v = round(&dec("123.456"), 1);
        assert_eq!(v.to_string(), "100");
        // round(0.0789, 1) = 0.08
        let v2 = round(&dec("0.0789"), 1);
        assert_eq!(v2.to_string(), "0.08");
    }

    #[test]
    fn from_f64_very_large() {
        // Very large f64 should convert without panic
        let v = from_f64(1e300);
        let f = to_f64(&v);
        assert!((f - 1e300).abs() / 1e300 < 1e-10, "large f64 round-trip: {}", f);
    }

    #[test]
    fn eval_decimal_missing_variable_errors() {
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("x + 1").unwrap();
        assert!(eval_decimal(&e, &vars, 10).is_err());
    }

    #[test]
    fn eval_decimal_wrong_arity_errors() {
        let vars = HashMap::new();
        let e = crate::parser::Parser::parse("sin(1, 2)").unwrap();
        assert!(eval_decimal(&e, &vars, 10).is_err());
    }

    #[test]
    fn log2_and_log10_known_values() {
        // log2(8) = 3
        let v = log2(&BigDecimal::from(8), 30).unwrap();
        assert!(close_dec(&v, 3.0, 1e-9), "log2(8) = 3: {}", to_f64(&v));

        // log10(1000) = 3
        let v2 = log10(&BigDecimal::from(1000), 30).unwrap();
        assert!(close_dec(&v2, 3.0, 1e-9), "log10(1000) = 3: {}", to_f64(&v2));
    }

    #[test]
    fn atan_known_values() {
        // atan(1) = π/4
        let v = atan(&BigDecimal::from(1), 30).unwrap();
        assert!(close_dec(&v, std::f64::consts::FRAC_PI_4, 1e-9), "atan(1) = π/4: {}", to_f64(&v));

        // atan(0) = 0
        let v2 = atan(&BigDecimal::from(0), 30).unwrap();
        assert!(close_dec(&v2, 0.0, 1e-9), "atan(0) = 0: {}", to_f64(&v2));
    }

    #[test]
    fn sqrt_of_zero() {
        let v = sqrt(&BigDecimal::from(0), 30).unwrap();
        assert_eq!(to_f64(&v), 0.0);
    }

    #[test]
    fn div_by_zero_errors() {
        let v = div(&BigDecimal::from(1), &BigDecimal::from(0), 30);
        assert!(v.is_err(), "division by zero should error");
    }
}
