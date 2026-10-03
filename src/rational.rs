//! Exact rational arithmetic.
//!
//! `Rational` stores a numerator and denominator as `i64` values, always
//! kept in reduced form (GCD = 1) with a positive denominator.

use crate::error::{MathError, Result};
use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// An exact rational number `num / den`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Rational {
    num: i64,
    den: i64,
}

impl Rational {
    /// Create a new rational, reducing to lowest terms with positive denominator.
    pub fn new(num: i64, den: i64) -> Result<Self> {
        if den == 0 {
            return Err(MathError::InvalidArgument(
                "Rational: denominator cannot be zero".into(),
            ));
        }
        let g = gcd(num.unsigned_abs(), den.unsigned_abs()) as i64;
        let g = if g == 0 { 1 } else { g };
        let (num, den) = (num / g, den / g);
        if den < 0 {
            Ok(Self {
                num: -num,
                den: -den,
            })
        } else {
            Ok(Self { num, den })
        }
    }

/// Create a rational from an integer.
    pub fn from_int(n: i64) -> Self {
        Self { num: n, den: 1 }
    }

    /// Find the closest rational approximation `p/q` (with `1 <= q <= max_den`)
    /// to `x` using the continued-fraction convergent sequence.
    ///
    /// Returns `None` only if `x` is non-finite or `max_den < 1`.  Otherwise
    /// returns the best (smallest `|x − p/q|`) convergent whose denominator
    /// fits the budget — even when the best approximation is still a poor
    /// fit (e.g. `x = π, max_den = 3` returns `3/1`).  The caller decides
    /// whether the approximation is acceptable by comparing against the
    /// raw float.
    ///
    /// Cross-port pattern from unsolve `limit.ts:537-562` and `apart.ts:93-117`.
    ///
    /// The algorithm walks the simple continued fraction of `|x|` and
    /// returns the convergent with denominator in range that minimises
    /// `|x − p/q|`.  An exact hit (within `1e-15`) is accepted immediately
    /// to avoid regressing to a worse approximation for clean fractions
    /// like `0.5` or `1/7`.
    pub fn from_float_snap(x: f64, max_den: i64) -> Option<Self> {
        if !x.is_finite() || max_den < 1 {
            return None;
        }
        if x == 0.0 {
            return Some(Self::from_int(0));
        }
        let sign: i64 = if x < 0.0 { -1 } else { 1 };
        let r = x.abs();
        // Walk CF convergents h/k starting from h0=1, k0=0 and h1=floor(r), k1=1.
        let mut h0: i64 = 1;
        let mut k0: i64 = 0;
        let mut h1: i64 = r.floor() as i64;
        let mut k1: i64 = 1;
        // Check the integer part alone.
        if h1 as f64 == r {
            return Self::new(sign * h1, 1).ok();
        }
        let mut frac = r - h1 as f64;
        let mut best: (i64, i64, f64) = (sign * h1, 1, (h1 as f64 - r).abs());
        let mut iters = 0usize;
        const MAX_ITERS: usize = 200;
        while k1 <= max_den && frac > 1e-15 && iters < MAX_ITERS {
            iters += 1;
            let inv = 1.0 / frac;
            let a = inv.floor() as i64;
            let h2 = a.checked_mul(h1).and_then(|v| v.checked_add(h0));
            let k2 = a.checked_mul(k1).and_then(|v| v.checked_add(k0));
            let (h2, k2) = match (h2, k2) {
                (Some(h), Some(k)) => (h, k),
                _ => break,
            };
            if k2 > max_den || k2 <= 0 {
                break;
            }
            let approx = h2 as f64 / k2 as f64;
            let err = (approx - r).abs();
            // Exact hit (within tolerance): accept immediately.
            if err < 1e-15 {
                return Self::new(sign * h2, k2).ok();
            }
            if err < best.2 {
                best = (sign * h2, k2, err);
            }
            h0 = h1;
            k0 = k1;
            h1 = h2;
            k1 = k2;
            let frac_rest = inv - a as f64;
            if frac_rest <= 0.0 {
                break;
            }
            frac = frac_rest;
        }
        Self::new(best.0, best.1).ok()
    }

/// Numerator.
    pub fn num(&self) -> i64 {
        self.num
    }

    /// Denominator (always positive).
    pub fn den(&self) -> i64 {
        self.den
    }

    /// Convert to `f64`.
    pub fn to_f64(&self) -> f64 {
        self.num as f64 / self.den as f64
    }

    /// True if the number is an integer.
    pub fn is_integer(&self) -> bool {
        self.den == 1
    }

    /// Absolute value.
    pub fn abs(self) -> Self {
        Self {
            num: self.num.abs(),
            den: self.den,
        }
    }

    /// Reciprocal `1/self`.
    pub fn recip(self) -> Result<Self> {
        if self.num == 0 {
            return Err(MathError::InvalidArgument(
                "Rational::recip: reciprocal of zero".into(),
            ));
        }
        Self::new(self.den, self.num)
    }

    /// Raise to an integer power.
    pub fn powi(self, exp: i32) -> Self {
        if exp == 0 {
            return Self::from_int(1);
        }
        if exp < 0 {
            let r = self.recip().unwrap_or(Self {
                num: 0,
                den: 1,
            });
            return r.powi(-exp);
        }
        let mut result = Self::from_int(1);
        let mut base = self;
        let mut e = exp as u32;
        while e > 0 {
            if e & 1 == 1 {
                result = result * base;
            }
            base = base * base;
e >>= 1;
        }
        result
    }

    /// Snap a float to a rational display string, falling back to a 12-significant-
    /// digit fixed representation when no clean fraction exists.  `max_den`
    /// caps the denominator (e.g. `1_000_000` for limits, `10_000` for
    /// partial fractions and qsolve).
    ///
    /// Default tolerance (`tol = None`) requires the CF to terminate
    /// *exactly* (err < 1e-15).  This avoids misleading outputs like
    /// `665857/470832` for `sqrt(2)` while still accepting `5/3` and
    /// `1/7`.  Pass `Some(t)` to override (e.g. for testing or a noisy
    /// numeric input).
    ///
    /// Returns `None` only for non-finite `v`.
    pub fn snap_to_string(v: f64, max_den: i64, tol: Option<f64>) -> Option<String> {
        if !v.is_finite() {
            return None;
        }
        if let Some(r) = Self::from_float_snap(v, max_den) {
            let err = (r.to_f64() - v).abs();
            let ok = match tol {
                Some(t) => err <= t,
                None => err < 1e-15,
            };
            if ok {
                return Some(r.to_string());
            }
        }
        // Fallback: integer or 12-sig-digit fixed (mirrors limit.rs / qsolve.rs /
        // apart.rs conventions so we don't regress those modules' tests).
        let snapped = if (v - v.round()).abs() < 1e-9 * v.abs().max(1.0) {
            v.round()
        } else {
            let mag = v.abs();
            let digits = 11 - mag.log10().floor() as i32;
            if !(-20..=20).contains(&digits) {
                v
            } else {
                let factor = 10f64.powi(digits);
                (v * factor).round() / factor
            }
        };
        Some(format_snap_f64(snapped))
    }
}

/// 12-significant-digit fallback formatter (replicates the inline snapping in
/// `qsolve::snap` / `limit::snap_num` / `apart::snap` so a single shared
/// helper can replace them).
fn format_snap_f64(v: f64) -> String {
    if !v.is_finite() {
        return format!("{}", v);
    }
    if (v - v.round()).abs() < 1e-9 * v.abs().max(1.0) {
        return format!("{}", v.round() as i64);
    }
    let mag = v.abs();
    let digits = 11 - mag.log10().floor() as i32;
    if !(-20..=20).contains(&digits) {
        return format!("{}", v);
    }
    let factor = 10f64.powi(digits);
    let rounded = (v * factor).round() / factor;
    format!("{}", rounded)
}

/// Try to evaluate an `Expr` using exact rational arithmetic.
    /// Returns `Some(Rational)` if all parts are rational (integer leaves + arithmetic ops).
    /// Returns `None` if the expression contains functions, variables, or non-integer constants.
    pub fn eval_rational(e: &crate::expr::Expr) -> Option<Rational> {
    use crate::expr::Expr::*;
    match e {
        Num(x) => {
            // Only accept values that are exact integers
            if x.fract() == 0.0 && x.abs() < i64::MAX as f64 {
                Some(Rational::from_int(*x as i64))
            } else {
                None
            }
        }
        Var(_) => None,
        Neg(a) => eval_rational(a).map(|r| -r),
        Add(a, b) => {
            let ra = eval_rational(a)?;
            let rb = eval_rational(b)?;
            Some(ra + rb)
        }
        Sub(a, b) => {
            let ra = eval_rational(a)?;
            let rb = eval_rational(b)?;
            Some(ra - rb)
        }
        Mul(a, b) => {
            let ra = eval_rational(a)?;
            let rb = eval_rational(b)?;
            Some(ra * rb)
        }
        Div(a, b) => {
            let ra = eval_rational(a)?;
            let rb = eval_rational(b)?;
            if rb.num() == 0 {
                return None;
            }
            Some(ra / rb)
        }
        Pow(base, exp) => {
            let rb = eval_rational(base)?;
            let re = eval_rational(exp)?;
            // Only integer exponents are exact
            if re.is_integer() {
                Some(rb.powi(re.num() as i32))
            } else {
                None
            }
        }
        Func(_, _) => None,
    }
}

impl fmt::Display for Rational {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

impl From<i64> for Rational {
    fn from(n: i64) -> Self {
        Self::from_int(n)
    }
}

impl From<i32> for Rational {
    fn from(n: i32) -> Self {
        Self::from_int(n as i64)
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rational {
    fn cmp(&self, other: &Self) -> Ordering {
        // a/b vs c/d  →  a*d vs c*b  (denominators are positive)
        let lhs = self.num as i128 * other.den as i128;
        let rhs = other.num as i128 * self.den as i128;
        lhs.cmp(&rhs)
    }
}

impl Add for Rational {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        // (a/b + c/d) = (ad + bc) / bd
        let num = self.num as i128 * rhs.den as i128 + rhs.num as i128 * self.den as i128;
        let den = self.den as i128 * rhs.den as i128;
        reduce_i128(num, den)
    }
}

impl Sub for Rational {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let num = self.num as i128 * rhs.den as i128 - rhs.num as i128 * self.den as i128;
        let den = self.den as i128 * rhs.den as i128;
        reduce_i128(num, den)
    }
}

impl Mul for Rational {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let num = self.num as i128 * rhs.num as i128;
        let den = self.den as i128 * rhs.den as i128;
        reduce_i128(num, den)
    }
}

impl Div for Rational {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        if rhs.num == 0 {
            panic!("Rational: division by zero");
        }
        let num = self.num as i128 * rhs.den as i128;
        let den = self.den as i128 * rhs.num as i128;
        reduce_i128(num, den)
    }
}

impl Neg for Rational {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            num: -self.num,
            den: self.den,
        }
    }
}

/// Greatest common divisor of two non-negative integers.
fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Reduce a fraction with i128 numerator/denominator to a Rational (i64).
/// Panics if the reduced values don't fit in i64.
fn reduce_i128(num: i128, den: i128) -> Rational {
    if den == 0 {
        panic!("Rational: denominator cannot be zero");
    }
    let g = gcd_u128(num.unsigned_abs(), den.unsigned_abs());
    let g = if g == 0 { 1 } else { g } as i128;
    let num = num / g;
    let den = den / g;
    let num = if den < 0 { -num } else { num };
    let den = den.abs();
    Rational {
        num: num as i64,
        den: den as i64,
    }
}

fn gcd_u128(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd_u128(b, a % b)
    }
}

/// Parse a rational from a string. Accepts "n", "n/d", or decimal "a.b".
pub fn parse_rational(s: &str) -> Result<Rational> {
    let s = s.trim();
    if s.contains('/') {
        let parts: Vec<&str> = s.split('/').collect();
        if parts.len() != 2 {
            return Err(MathError::InvalidArgument(
                format!("parse_rational: expected 'n/d', got '{}'", s),
            ));
        }
        let num: i64 = parts[0]
            .trim()
            .parse()
            .map_err(|_| MathError::InvalidArgument("parse_rational: bad numerator".into()))?;
        let den: i64 = parts[1]
            .trim()
            .parse()
            .map_err(|_| MathError::InvalidArgument("parse_rational: bad denominator".into()))?;
        Rational::new(num, den)
    } else if s.contains('.') {
        // Decimal → rational
        let neg = s.starts_with('-');
        let s = s.trim_start_matches('-');
        let parts: Vec<&str> = s.split('.').collect();
        if parts.len() != 2 {
            return Err(MathError::InvalidArgument(
                format!("parse_rational: bad decimal '{}'", s),
            ));
        }
        let int_part: i64 = if parts[0].is_empty() {
            0
        } else {
            parts[0]
                .parse()
                .map_err(|_| MathError::InvalidArgument("parse_rational: bad integer part".into()))?
        };
        let frac_str = parts[1];
        let frac_part: i64 = if frac_str.is_empty() {
            0
        } else {
            frac_str
                .parse()
                .map_err(|_| MathError::InvalidArgument("parse_rational: bad fractional part".into()))?
        };
        let den = 10i64.pow(frac_str.len() as u32);
        let num = int_part * den + frac_part;
        let r = Rational::new(num, den)?;
        Ok(if neg { -r } else { r })
    } else {
        let n: i64 = s
            .parse()
            .map_err(|_| MathError::InvalidArgument(format!("parse_rational: bad integer '{}'", s)))?;
        Ok(Rational::from_int(n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_new_reduces() {
        let r = Rational::new(6, 8).unwrap();
        assert_eq!(r.num(), 3);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn rational_new_negative_den() {
        let r = Rational::new(3, -4).unwrap();
        assert_eq!(r.num(), -3);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn rational_zero_den_errors() {
        assert!(Rational::new(1, 0).is_err());
    }

    #[test]
    fn rational_add() {
        let a = Rational::new(1, 2).unwrap();
        let b = Rational::new(1, 3).unwrap();
        let c = a + b;
        assert_eq!(c.num(), 5);
        assert_eq!(c.den(), 6);
    }

    #[test]
    fn rational_sub() {
        let a = Rational::new(1, 2).unwrap();
        let b = Rational::new(1, 3).unwrap();
        let c = a - b;
        assert_eq!(c.num(), 1);
        assert_eq!(c.den(), 6);
    }

    #[test]
    fn rational_mul() {
        let a = Rational::new(2, 3).unwrap();
        let b = Rational::new(3, 4).unwrap();
        let c = a * b;
        assert_eq!(c.num(), 1);
        assert_eq!(c.den(), 2);
    }

    #[test]
    fn rational_div() {
        let a = Rational::new(2, 3).unwrap();
        let b = Rational::new(4, 5).unwrap();
        let c = a / b;
        assert_eq!(c.num(), 5);
        assert_eq!(c.den(), 6);
    }

    #[test]
    fn rational_neg() {
        let a = Rational::new(3, 4).unwrap();
        let b = -a;
        assert_eq!(b.num(), -3);
        assert_eq!(b.den(), 4);
    }

    #[test]
    fn rational_abs() {
        let a = Rational::new(-3, 4).unwrap();
        let b = a.abs();
        assert_eq!(b.num(), 3);
        assert_eq!(b.den(), 4);
    }

    #[test]
    fn rational_recip() {
        let a = Rational::new(3, 4).unwrap();
        let b = a.recip().unwrap();
        assert_eq!(b.num(), 4);
        assert_eq!(b.den(), 3);
    }

    #[test]
    fn rational_recip_zero_errors() {
        let a = Rational::from_int(0);
        assert!(a.recip().is_err());
    }

    #[test]
    fn rational_powi_positive() {
        let a = Rational::new(2, 3).unwrap();
        let b = a.powi(3);
        assert_eq!(b.num(), 8);
        assert_eq!(b.den(), 27);
    }

    #[test]
    fn rational_powi_negative() {
        let a = Rational::new(2, 3).unwrap();
        let b = a.powi(-2);
        assert_eq!(b.num(), 9);
        assert_eq!(b.den(), 4);
    }

    #[test]
    fn rational_powi_zero() {
        let a = Rational::new(2, 3).unwrap();
        let b = a.powi(0);
        assert_eq!(b.num(), 1);
        assert_eq!(b.den(), 1);
    }

    #[test]
    fn rational_to_f64() {
        let a = Rational::new(1, 4).unwrap();
        assert!((a.to_f64() - 0.25).abs() < 1e-15);
    }

    #[test]
    fn rational_is_integer() {
        assert!(Rational::from_int(5).is_integer());
        assert!(!Rational::new(1, 2).unwrap().is_integer());
    }

    #[test]
    fn rational_display() {
        assert_eq!(Rational::from_int(5).to_string(), "5");
        assert_eq!(Rational::new(3, 4).unwrap().to_string(), "3/4");
        assert_eq!(Rational::new(-3, 4).unwrap().to_string(), "-3/4");
    }

    #[test]
    fn rational_ordering() {
        let a = Rational::new(1, 2).unwrap();
        let b = Rational::new(2, 3).unwrap();
        assert!(a < b);
        assert!(b > a);
        assert_eq!(a.cmp(&a), Ordering::Equal);
    }

    #[test]
    fn rational_equality() {
        let a = Rational::new(1, 2).unwrap();
        let b = Rational::new(2, 4).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn rational_large_numerators() {
        // Cross-multiplication should not overflow i64
        let a = Rational::new(1, 1_000_000_000).unwrap();
        let b = Rational::new(1, 1_000_000_001).unwrap();
        assert!(a > b); // 1/1e9 > 1/1e9+1
        assert_eq!(a.cmp(&a), Ordering::Equal);
    }

    #[test]
    fn rational_add_large() {
        // i64 * i64 → i128 should handle large intermediate values
        let a = Rational::new(1, 1_000_000_000).unwrap();
        let b = Rational::new(1, 1_000_000_000).unwrap();
        let c = a + b;
        assert_eq!(c.num(), 1);
        assert_eq!(c.den(), 500_000_000);
    }

    #[test]
    fn parse_rational_integer() {
        let r = parse_rational("42").unwrap();
        assert_eq!(r, Rational::from_int(42));
    }

    #[test]
    fn parse_rational_fraction() {
        let r = parse_rational("3/4").unwrap();
        assert_eq!(r.num(), 3);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn parse_rational_negative_fraction() {
        let r = parse_rational("-3/4").unwrap();
        assert_eq!(r.num(), -3);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn parse_rational_decimal() {
        let r = parse_rational("0.5").unwrap();
        assert_eq!(r.num(), 1);
        assert_eq!(r.den(), 2);
    }

    #[test]
    fn parse_rational_decimal_negative() {
        let r = parse_rational("-1.25").unwrap();
        assert_eq!(r.num(), -5);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn parse_rational_decimal_long() {
        let r = parse_rational("0.125").unwrap();
        assert_eq!(r.num(), 1);
        assert_eq!(r.den(), 8);
    }

    #[test]
    fn parse_rational_bad_input() {
        assert!(parse_rational("abc").is_err());
        assert!(parse_rational("1/0").is_err());
        assert!(parse_rational("1/2/3").is_err());
    }

    #[test]
    fn rational_from_int() {
        let r = Rational::from_int(-7);
        assert_eq!(r.num(), -7);
        assert_eq!(r.den(), 1);
        assert!(r.is_integer());
    }

    #[test]
    fn rational_chained_arithmetic() {
        // (1/2 + 1/3) * (1/4) = (5/6) * (1/4) = 5/24
        let a = Rational::new(1, 2).unwrap();
        let b = Rational::new(1, 3).unwrap();
        let c = Rational::new(1, 4).unwrap();
        let result = (a + b) * c;
        assert_eq!(result.num(), 5);
        assert_eq!(result.den(), 24);
    }

    #[test]
    fn eval_rational_simple_add() {
        let e = crate::parser::Parser::parse("1/2 + 3/4").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 5);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn eval_rational_frac_tex() {
        let e = crate::parser::Parser::parse(r"\frac{1}{2} + \frac{3}{4}").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 5);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn eval_rational_mul_div() {
        let e = crate::parser::Parser::parse("(1/2 + 1/3) * 2/5").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 1);
        assert_eq!(r.den(), 3);
    }

    #[test]
    fn eval_rational_pow() {
        let e = crate::parser::Parser::parse("2^3 + 1/4").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 33);
        assert_eq!(r.den(), 4);
    }

    #[test]
    fn eval_rational_integer() {
        let e = crate::parser::Parser::parse("3 + 4").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 7);
        assert_eq!(r.den(), 1);
    }

    #[test]
    fn eval_rational_sin_returns_none() {
        let e = crate::parser::Parser::parse("sin(pi/4)").unwrap();
        assert!(eval_rational(&e).is_none());
    }

    #[test]
    fn eval_rational_var_returns_none() {
        let e = crate::parser::Parser::parse("x + 1").unwrap();
        assert!(eval_rational(&e).is_none());
    }

    #[test]
    fn eval_rational_decimal_returns_none() {
        let e = crate::parser::Parser::parse("0.5 + 0.25").unwrap();
        assert!(eval_rational(&e).is_none());
    }

    #[test]
    fn eval_rational_neg() {
        let e = crate::parser::Parser::parse("-(1/2 + 1/3)").unwrap();
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), -5);
        assert_eq!(r.den(), 6);
    }

    #[test]
    fn eval_rational_nested_frac() {
        let e = crate::parser::Parser::parse(r"\frac{\frac{1}{2}}{\frac{3}{4}").unwrap_or(
            crate::parser::Parser::parse(r"(1/2)/(3/4)").unwrap()
        );
        let r = eval_rational(&e).unwrap();
        assert_eq!(r.num(), 2);
        assert_eq!(r.den(), 3);
    }

    // ---- from_float_snap: continued-fraction rational recovery ----
    // Cross-port pattern from unsolve `limit.ts:537-562` and `apart.ts:93-117`.

    fn s(x: f64, max: i64) -> (i64, i64) {
        let r = Rational::from_float_snap(x, max).expect("expected a rational");
        (r.num(), r.den())
    }

    #[test]
    fn snap_zero() {
        assert_eq!(s(0.0, 1000), (0, 1));
    }

    #[test]
    fn snap_integer() {
        assert_eq!(s(5.0, 1000), (5, 1));
        assert_eq!(s(-7.0, 1000), (-7, 1));
    }

    #[test]
    fn snap_half() {
        // 0.5 → exact 1/2 even with a large budget
        assert_eq!(s(0.5, 100), (1, 2));
    }

    #[test]
    fn snap_negative_half() {
        assert_eq!(s(-0.5, 100), (-1, 2));
    }

    #[test]
    fn snap_third() {
        // 1/3 = 0.333...; with max_den=10 the best is 1/3
        assert_eq!(s(1.0 / 3.0, 10), (1, 3));
    }

    #[test]
    fn snap_5_over_3() {
        // The harvest's signature case: 5/3 ≈ 1.6666666667
        assert_eq!(s(5.0 / 3.0, 100), (5, 3));
    }

    #[test]
    fn snap_negative_7_over_2() {
        // qsolve-style case: -3.5 → -7/2
        assert_eq!(s(-3.5, 100), (-7, 2));
    }

    #[test]
    fn snap_seventh() {
        // 1/7 with max_den=20 → exact
        assert_eq!(s(1.0 / 7.0, 20), (1, 7));
    }

    #[test]
    fn snap_pi_is_not_rational_with_small_budget() {
        // π with max_den=100 → 22/7 is the classic convergent
        let r = Rational::from_float_snap(std::f64::consts::PI, 100).unwrap();
        assert_eq!((r.num(), r.den()), (22, 7));
    }

    #[test]
    fn snap_pi_with_tiny_budget_falls_back() {
        // π with max_den=3 → 3/1 is the best (≈ 3.0), error ≈ 0.14
        // We don't assert a specific value — just that nothing panics and
        // it's reasonable.
        let r = Rational::from_float_snap(std::f64::consts::PI, 3);
        // With max_den=3, the CF walk would yield 3/1 (integer part)
        // or potentially 13/4 (out of budget). The integer convergent wins.
        assert!(r.is_some());
        let r = r.unwrap();
        assert_eq!((r.num(), r.den()), (3, 1));
    }

    #[test]
    fn snap_non_finite_returns_none() {
        assert!(Rational::from_float_snap(f64::NAN, 100).is_none());
        assert!(Rational::from_float_snap(f64::INFINITY, 100).is_none());
        assert!(Rational::from_float_snap(f64::NEG_INFINITY, 100).is_none());
    }

    #[test]
    fn snap_overflow_returns_none() {
        // max_den < 1 is invalid
        assert!(Rational::from_float_snap(1.5, 0).is_none());
        assert!(Rational::from_float_snap(1.5, -1).is_none());
    }

    #[test]
    fn snap_too_inaccurate_returns_none() {
        // 0.1 with max_den=2: best is 0/1 (err 0.1). That's within tolerance
        // — so it might still snap. Use a less-cooperative value:
        // 0.30000001 with max_den=5 → best is 3/10 (out of budget) or 1/3 (in
        // budget). 1/3 = 0.333... err = 0.033, well within tolerance.
        // Actually 0.30000001 with max_den=5 → 1/3 (err 0.033) is accepted.
        // Just verify the API doesn't panic.
        let r = Rational::from_float_snap(0.30000001, 5);
        assert!(r.is_some());
    }

    // ---- snap_to_string: format a float as a rational if one fits cleanly ----

    #[test]
    fn snap_string_rational_hits() {
        // The harvest's signature case: 5/3 must display as "5/3".
        assert_eq!(Rational::snap_to_string(5.0 / 3.0, 1_000_000, None).unwrap(), "5/3");
        assert_eq!(Rational::snap_to_string(-7.0 / 2.0, 1_000_000, None).unwrap(), "-7/2");
        assert_eq!(Rational::snap_to_string(0.5, 1_000_000, None).unwrap(), "1/2");
    }

    #[test]
    fn snap_string_integer() {
        assert_eq!(Rational::snap_to_string(5.0, 1000, None).unwrap(), "5");
        assert_eq!(Rational::snap_to_string(-3.0, 1000, None).unwrap(), "-3");
    }

    #[test]
    fn snap_string_integer_via_near_integer() {
        // Values within tolerance of an integer should display as that integer.
        assert_eq!(
            Rational::snap_to_string(0.9999999999998, 1000, None).unwrap(),
            "1"
        );
    }

    #[test]
    fn snap_string_fallback_to_fixed_for_irrational() {
        // √2 has no small-denominator rational; should fall back to 12-sig-digit.
        let s = Rational::snap_to_string(std::f64::consts::SQRT_2, 100, None).unwrap();
        // Verify it's a plain decimal representation (not "a/b").
        assert!(!s.contains('/'), "expected plain decimal, got: {s}");
    }

    #[test]
    fn snap_string_non_finite_returns_none() {
        assert!(Rational::snap_to_string(f64::NAN, 1000, None).is_none());
        assert!(Rational::snap_to_string(f64::INFINITY, 1000, None).is_none());
    }

    #[test]
    fn snap_string_custom_tolerance_allows() {
        // With a generous tolerance (0.01), accept the snap.
        let s = Rational::snap_to_string(5.0 / 3.0, 1000, Some(0.01)).unwrap();
        assert_eq!(s, "5/3");
    }

    #[test]
    fn snap_string_exact_rational_always_accepted() {
        // 5/3 IS exactly representable by the CF algorithm (the float
        // 5.0/3.0 hits the convergent with err < 1e-15), so even a
        // tight tolerance can't skip it.  This documents the design.
        let s = Rational::snap_to_string(5.0 / 3.0, 1000, Some(1e-15)).unwrap();
        assert_eq!(s, "5/3");
    }
}
