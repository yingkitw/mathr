//! Interval arithmetic for rigorous bounds on computations over the [`Expr`] AST.
//!
//! An [`Interval`] `[lo, hi]` represents the set of all real numbers `x` with
//! `lo ≤ x ≤ hi`. Arithmetic operations are extended to intervals so that the
//! result always contains every possible output value — giving **guaranteed
//! bounds** rather than point estimates.
//!
//! ## Example
//! ```
//! use mathr::interval::{Interval, eval_interval};
//! use mathr::parser::Parser;
//! use std::collections::HashMap;
//!
//! let e = Parser::parse("x^2 + 1").unwrap();
//! let mut vars = HashMap::new();
//! vars.insert("x".to_string(), Interval::new(-2.0, 3.0));
//! let bounds = eval_interval(&e, &vars).unwrap();
//! // x^2 ∈ [0, 9] on [-2, 3], so x^2 + 1 ∈ [1, 10]
//! assert!(bounds.contains(1.0));
//! assert!(bounds.contains(10.0));
//! ```
//!
//! ## Limitations
//!
//! - **No outward rounding**: bounds use plain `f64` arithmetic without
//!   IEEE 1788 directed rounding, so they may be slightly tight at the last
//!   bit. For safety-critical use, widen results by a small epsilon.
//! - **Dependency problem**: `x - x` over `[a, b]` yields `[a-b, b-a]`, not
//!   `[0, 0]`, because interval arithmetic treats each occurrence as
//!   independent. This is a fundamental limitation of the method.
//! - **Division by zero-containing intervals**: returns the whole real line
//!   `(-∞, +∞)` rather than a split interval.

use std::collections::HashMap;
use std::ops::{Add, Div, Mul, Neg, Sub};

use crate::error::{MathError, Result};
use crate::expr::Expr;

/// A closed real interval `[lo, hi]`, or the empty interval when `is_empty`
/// is set. The whole real line is `(-∞, +∞)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub lo: f64,
    pub hi: f64,
    pub is_empty: bool,
}

impl Interval {
    /// Create a closed interval `[lo, hi]`. If `lo > hi`, returns the empty
    /// interval.
    pub fn new(lo: f64, hi: f64) -> Self {
        if lo > hi {
            Self::empty()
        } else {
            Self { lo, hi, is_empty: false }
        }
    }

    /// A degenerate (point) interval `[x, x]`.
    pub fn point(x: f64) -> Self {
        Self { lo: x, hi: x, is_empty: false }
    }

    /// The empty interval ∅.
    pub fn empty() -> Self {
        Self { lo: f64::NAN, hi: f64::NAN, is_empty: true }
    }

    /// The whole real line `(-∞, +∞)`.
    pub fn whole() -> Self {
        Self { lo: f64::NEG_INFINITY, hi: f64::INFINITY, is_empty: false }
    }

    /// Width of the interval: `hi - lo`. Empty → NaN.
    pub fn width(&self) -> f64 {
        if self.is_empty {
            f64::NAN
        } else {
            self.hi - self.lo
        }
    }

    /// Midpoint `(lo + hi) / 2`. Empty → NaN.
    pub fn midpoint(&self) -> f64 {
        if self.is_empty {
            f64::NAN
        } else {
            (self.lo + self.hi) / 2.0
        }
    }

    /// Does this interval contain `x`?
    pub fn contains(&self, x: f64) -> bool {
        !self.is_empty && x >= self.lo && x <= self.hi
    }

    /// Do two intervals overlap (have non-empty intersection)?
    pub fn overlaps(&self, other: &Self) -> bool {
        if self.is_empty || other.is_empty {
            return false;
        }
        self.lo <= other.hi && other.lo <= self.hi
    }

    /// Intersection of two intervals. Returns empty if disjoint.
    pub fn intersect(&self, other: &Self) -> Self {
        if self.is_empty || other.is_empty {
            return Self::empty();
        }
        let lo = self.lo.max(other.lo);
        let hi = self.hi.min(other.hi);
        if lo <= hi {
            Self::new(lo, hi)
        } else {
            Self::empty()
        }
    }

    /// Convex hull (union) of two intervals — the smallest interval containing
    /// both.
    pub fn hull(&self, other: &Self) -> Self {
        if self.is_empty {
            return *other;
        }
        if other.is_empty {
            return *self;
        }
        Self::new(self.lo.min(other.lo), self.hi.max(other.hi))
    }

    /// Square: `x²` over the interval. Tighter than `self * self` because it
    /// knows both operands are the same value.
    pub fn sqr(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        let lo_sq = self.lo * self.lo;
        let hi_sq = self.hi * self.hi;
        if self.contains(0.0) {
            Self::new(0.0, lo_sq.max(hi_sq))
        } else {
            Self::new(lo_sq.min(hi_sq), lo_sq.max(hi_sq))
        }
    }

    /// Integer power `x^n` for `n ≥ 0`.
    pub fn powi(&self, n: u32) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if n == 0 {
            return Self::point(1.0);
        }
        if n == 1 {
            return *self;
        }
        let lo_p = self.lo.powi(n as i32);
        let hi_p = self.hi.powi(n as i32);
        if n.is_multiple_of(2) {
            // Even power: non-negative. If 0 ∈ [lo, hi], min is 0.
            if self.contains(0.0) {
                Self::new(0.0, lo_p.max(hi_p))
            } else {
                Self::new(lo_p.min(hi_p), lo_p.max(hi_p))
            }
        } else {
            // Odd power: monotonic increasing.
            Self::new(lo_p, hi_p)
        }
    }

    /// Absolute value: `|x|` over the interval.
    pub fn abs(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.lo >= 0.0 {
            *self
        } else if self.hi <= 0.0 {
            Self::new(-self.hi, -self.lo)
        } else {
            // Straddles zero
            Self::new(0.0, self.lo.abs().max(self.hi.abs()))
        }
    }

    /// `exp(x)` — monotonically increasing.
    pub fn exp(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(self.lo.exp(), self.hi.exp())
    }

    /// `ln(x)` — monotonically increasing. Domain: `(0, ∞)`.
    pub fn ln(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.hi <= 0.0 {
            return Self::empty(); // entirely outside domain
        }
        let lo = if self.lo > 0.0 { self.lo.ln() } else { f64::NEG_INFINITY };
        Self::new(lo, self.hi.ln())
    }

    /// `log10(x)` — monotonically increasing. Domain: `(0, ∞)`.
    pub fn log10(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.hi <= 0.0 {
            return Self::empty();
        }
        let lo = if self.lo > 0.0 { self.lo.log10() } else { f64::NEG_INFINITY };
        Self::new(lo, self.hi.log10())
    }

    /// `log2(x)` — monotonically increasing. Domain: `(0, ∞)`.
    pub fn log2(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.hi <= 0.0 {
            return Self::empty();
        }
        let lo = if self.lo > 0.0 { self.lo.log2() } else { f64::NEG_INFINITY };
        Self::new(lo, self.hi.log2())
    }

    /// `sqrt(x)` — monotonically increasing. Domain: `[0, ∞)`.
    pub fn sqrt(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.hi < 0.0 {
            return Self::empty();
        }
        let lo = if self.lo > 0.0 { self.lo.sqrt() } else { 0.0 };
        Self::new(lo, self.hi.sqrt())
    }

    /// `sin(x)` over `[lo, hi]`. Tracks the global extrema ±1 when the
    /// interval spans a peak or trough.
    pub fn sin(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.width().is_infinite() {
            return Self::new(-1.0, 1.0);
        }
        // Sinusoid extrema: max at π/2 + 2kπ, min at -π/2 + 2kπ (i.e. 3π/2 + 2kπ).
        let has_max = contains_extremum(self.lo, self.hi, std::f64::consts::FRAC_PI_2);
        let has_min = contains_extremum(self.lo, self.hi, -std::f64::consts::FRAC_PI_2);
        let lo_s = self.lo.sin();
        let hi_s = self.hi.sin();
        let lo = if has_min { -1.0 } else { lo_s.min(hi_s) };
        let hi = if has_max { 1.0 } else { lo_s.max(hi_s) };
        Self::new(lo, hi)
    }

    /// `cos(x)` over `[lo, hi]`. Tracks the global extrema ±1.
    pub fn cos(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.width().is_infinite() {
            return Self::new(-1.0, 1.0);
        }
        // Cosine maxima at 2kπ, minima at π + 2kπ.
        let has_max = contains_extremum(self.lo, self.hi, 0.0);
        let has_min = contains_extremum(self.lo, self.hi, std::f64::consts::PI);
        let lo_c = self.lo.cos();
        let hi_c = self.hi.cos();
        let lo = if has_min { -1.0 } else { lo_c.min(hi_c) };
        let hi = if has_max { 1.0 } else { lo_c.max(hi_c) };
        Self::new(lo, hi)
    }

    /// `tan(x)` — has singularities at π/2 + kπ. Returns whole if a
    /// singularity lies in the interval.
    pub fn tan(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if interval_contains_tan_pole(self) {
            return Self::whole();
        }
        Self::new(self.lo.tan(), self.hi.tan())
    }

    /// `atan(x)` — monotonically increasing, bounded by `(-π/2, π/2)`.
    pub fn atan(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(self.lo.atan(), self.hi.atan())
    }

    /// `asin(x)` — monotonically increasing. Domain `[-1, 1]`.
    pub fn asin(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        let clamped = self.intersect(&Self::new(-1.0, 1.0));
        if clamped.is_empty {
            return Self::empty();
        }
        Self::new(clamped.lo.asin(), clamped.hi.asin())
    }

    /// `acos(x)` — monotonically **decreasing**. Domain `[-1, 1]`.
    pub fn acos(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        let clamped = self.intersect(&Self::new(-1.0, 1.0));
        if clamped.is_empty {
            return Self::empty();
        }
        // Decreasing: swap bounds.
        Self::new(clamped.hi.acos(), clamped.lo.acos())
    }

    /// `sinh(x)` — monotonically increasing.
    pub fn sinh(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(self.lo.sinh(), self.hi.sinh())
    }

    /// `cosh(x)` — even function, minimum at 0.
    pub fn cosh(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        if self.contains(0.0) {
            Self::new(1.0, self.lo.cosh().max(self.hi.cosh()))
        } else {
            Self::new(self.lo.cosh().min(self.hi.cosh()), self.lo.cosh().max(self.hi.cosh()))
        }
    }

    /// `tanh(x)` — monotonically increasing, bounded by `(-1, 1)`.
    pub fn tanh(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(self.lo.tanh(), self.hi.tanh())
    }

    /// `cbrt(x)` — cube root, monotonically increasing over all reals.
    pub fn cbrt(&self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(self.lo.cbrt(), self.hi.cbrt())
    }
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_empty {
            return write!(f, "∅");
        }
        write!(f, "[{}, {}]", fmt_bound(self.lo), fmt_bound(self.hi))
    }
}

fn fmt_bound(x: f64) -> String {
    if x.is_infinite() {
        if x > 0.0 { "∞".to_string() } else { "-∞".to_string() }
    } else if x == 0.0 {
        "0".to_string()
    } else if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else {
        format!("{}", x)
    }
}

// --- Operator overloads ---

impl Add for Interval {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        if self.is_empty || rhs.is_empty {
            return Self::empty();
        }
        Self::new(self.lo + rhs.lo, self.hi + rhs.hi)
    }
}

impl Sub for Interval {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        if self.is_empty || rhs.is_empty {
            return Self::empty();
        }
        Self::new(self.lo - rhs.hi, self.hi - rhs.lo)
    }
}

impl Mul for Interval {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        if self.is_empty || rhs.is_empty {
            return Self::empty();
        }
        let (a, b) = (self.lo, self.hi);
        let (c, d) = (rhs.lo, rhs.hi);
        let p1 = a * c;
        let p2 = a * d;
        let p3 = b * c;
        let p4 = b * d;
        let lo = p1.min(p2).min(p3).min(p4);
        let hi = p1.max(p2).max(p3).max(p4);
        Self::new(lo, hi)
    }
}

impl Div for Interval {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        if self.is_empty || rhs.is_empty {
            return Self::empty();
        }
        // If the divisor contains 0, the result is unbounded.
        if rhs.contains(0.0) {
            if rhs.lo == 0.0 && rhs.hi == 0.0 {
                // Division by exactly [0, 0] — undefined.
                return Self::empty();
            }
            return Self::whole();
        }
        let (a, b) = (self.lo, self.hi);
        let (c, d) = (rhs.lo, rhs.hi);
        let p1 = a / c;
        let p2 = a / d;
        let p3 = b / c;
        let p4 = b / d;
        let lo = p1.min(p2).min(p3).min(p4);
        let hi = p1.max(p2).max(p3).max(p4);
        Self::new(lo, hi)
    }
}

impl Neg for Interval {
    type Output = Self;
    fn neg(self) -> Self {
        if self.is_empty {
            return Self::empty();
        }
        Self::new(-self.hi, -self.lo)
    }
}

// --- Extrema helpers ---

/// Does `[lo, hi]` contain any value congruent to `target` modulo `2π`?
fn contains_extremum(lo: f64, hi: f64, target: f64) -> bool {
    let two_pi = 2.0 * std::f64::consts::PI;
    // Find the smallest k such that target + 2π·k ≥ lo.
    let k_start = ((lo - target) / two_pi).ceil() as i64 - 1;
    let k_end = ((hi - target) / two_pi).floor() as i64 + 1;
    for k in k_start..=k_end {
        let pt = target + two_pi * k as f64;
        if pt >= lo && pt <= hi {
            return true;
        }
    }
    false
}

/// Does the interval contain a pole of tan at `π/2 + kπ`?
fn interval_contains_tan_pole(iv: &Interval) -> bool {
    let pi = std::f64::consts::PI;
    let half = std::f64::consts::FRAC_PI_2;
    let period = pi;
    let k_start = ((iv.lo - half) / period).ceil() as i64 - 1;
    let k_end = ((iv.hi - half) / period).floor() as i64 + 1;
    for k in k_start..=k_end {
        let pole = half + period * k as f64;
        if iv.contains(pole) {
            return true;
        }
    }
    false
}

// =========================================================================
// Expr evaluation with intervals
// =========================================================================

/// Evaluate an [`Expr`] over interval-valued variables, returning rigorous
/// bounds on the result.
///
/// `vars` maps variable names to their intervals. Constants `pi`, `e`, `tau`
/// are provided as thin point intervals if not overridden. Unknown variables
/// produce an error. Functions are handled for the standard elementary set;
/// unknown functions produce an error.
pub fn eval_interval(expr: &Expr, vars: &HashMap<String, Interval>) -> Result<Interval> {
    match expr {
        Expr::Num(n) => Ok(Interval::point(*n)),
        Expr::Var(name) => {
            if let Some(iv) = vars.get(name) {
                Ok(*iv)
            } else {
                match name.as_str() {
                    "pi" => Ok(Interval::point(std::f64::consts::PI)),
                    "e" => Ok(Interval::point(std::f64::consts::E)),
                    "tau" => Ok(Interval::point(std::f64::consts::TAU)),
                    "inf" => Ok(Interval::point(f64::INFINITY)),
                    _ => Err(MathError::UnknownVariable(name.clone())),
                }
            }
        }
        Expr::Neg(a) => Ok(-eval_interval(a, vars)?),
        Expr::Add(a, b) => Ok(eval_interval(a, vars)? + eval_interval(b, vars)?),
        Expr::Sub(a, b) => Ok(eval_interval(a, vars)? - eval_interval(b, vars)?),
        Expr::Mul(a, b) => Ok(eval_interval(a, vars)? * eval_interval(b, vars)?),
        Expr::Div(a, b) => Ok(eval_interval(a, vars)? / eval_interval(b, vars)?),
        Expr::Pow(base, exp) => {
            let bv = eval_interval(base, vars)?;
            let ev = eval_interval(exp, vars)?;
            // If exponent is a non-negative integer point, use powi for tight bounds.
            if let Interval { lo, hi, is_empty: false } = ev {
                if lo == hi && lo >= 0.0 && lo.fract() == 0.0 && lo < 1e15 {
                    return Ok(bv.powi(lo as u32));
                }
            }
            // General case: exp(b·ln(a)) via intervals. Only valid for a > 0.
            if bv.lo <= 0.0 {
                return Err(MathError::Eval(format!(
                    "interval pow: base {} must be positive for non-integer exponent",
                    bv
                )));
            }
            Ok((bv.ln() * ev).exp())
        }
        Expr::Func(name, args) => {
            let ivs: Result<Vec<Interval>> =
                args.iter().map(|a| eval_interval(a, vars)).collect();
            let ivs = ivs?;
            eval_interval_func(name, &ivs)
        }
    }
}

fn eval_interval_func(name: &str, args: &[Interval]) -> Result<Interval> {
    let one = |args: &[Interval]| -> Result<Interval> {
        args.first()
            .copied()
            .ok_or_else(|| MathError::InvalidArgument(format!("{} needs 1 arg", name)))
    };
    match name {
        "sin" => Ok(one(args)?.sin()),
        "cos" => Ok(one(args)?.cos()),
        "tan" => Ok(one(args)?.tan()),
        "asin" => Ok(one(args)?.asin()),
        "acos" => Ok(one(args)?.acos()),
        "atan" => Ok(one(args)?.atan()),
        "sinh" => Ok(one(args)?.sinh()),
        "cosh" => Ok(one(args)?.cosh()),
        "tanh" => Ok(one(args)?.tanh()),
        "exp" => Ok(one(args)?.exp()),
        "ln" => Ok(one(args)?.ln()),
        "log" | "log10" => Ok(one(args)?.log10()),
        "log2" => Ok(one(args)?.log2()),
        "sqrt" => Ok(one(args)?.sqrt()),
        "abs" => Ok(one(args)?.abs()),
        "sqr" => Ok(one(args)?.sqr()),
        "cbrt" => {
            // cbrt is monotonically increasing over all reals.
            let a = one(args)?;
            if a.is_empty {
                Ok(Interval::empty())
            } else {
                Ok(Interval::new(a.lo.cbrt(), a.hi.cbrt()))
            }
        }
        "floor" | "ceil" | "round" | "sign" | "fract" => {
            // These are discontinuous; return the hull of applying to both endpoints.
            let a = one(args)?;
            if a.is_empty {
                return Ok(Interval::empty());
            }
            let f: fn(f64) -> f64 = match name {
                "floor" => f64::floor,
                "ceil" => f64::ceil,
                "round" => f64::round,
                "sign" => f64::signum,
                "fract" => f64::fract,
                _ => unreachable!(),
            };
            // For fract, the range is [0, 1) regardless. Be conservative.
            if name == "fract" {
                return Ok(Interval::new(0.0, 1.0));
            }
            if name == "sign" {
                // signum ∈ {-1, 0, 1}; hull over the interval.
                if a.contains(0.0) {
                    return Ok(Interval::new(-1.0, 1.0));
                }
                return Ok(Interval::point(if a.hi < 0.0 { -1.0 } else { 1.0 }));
            }
            let lo = f(a.lo);
            let hi = f(a.hi);
            Ok(Interval::new(lo.min(hi), lo.max(hi)))
        }
        "min" => {
            if args.is_empty() {
                return Err(MathError::InvalidArgument("min needs args".into()));
            }
            // min over intervals: take the hull of lower bounds → lower, and
            // the minimum of upper bounds → upper.
            let lo = args.iter().map(|a| a.lo).fold(f64::INFINITY, f64::min);
            let hi = args.iter().map(|a| a.hi).fold(f64::INFINITY, f64::min);
            Ok(Interval::new(lo, hi))
        }
        "max" => {
            if args.is_empty() {
                return Err(MathError::InvalidArgument("max needs args".into()));
            }
            let lo = args.iter().map(|a| a.lo).fold(f64::NEG_INFINITY, f64::max);
            let hi = args.iter().map(|a| a.hi).fold(f64::NEG_INFINITY, f64::max);
            Ok(Interval::new(lo, hi))
        }
        "pow" => {
            // pow(base, exp) — same logic as Expr::Pow.
            let (b, e) = two(args, name)?;
            if let Interval { lo, hi, is_empty: false } = e {
                if lo == hi && lo >= 0.0 && lo.fract() == 0.0 && lo < 1e15 {
                    return Ok(b.powi(lo as u32));
                }
            }
            if b.lo <= 0.0 {
                return Err(MathError::Eval(format!(
                    "interval pow: base {} must be positive for non-integer exponent",
                    b
                )));
            }
            Ok((b.ln() * e).exp())
        }
        "mod" => {
            // mod(a, b) — conservative hull over [0, |b|).
            let (_, b) = two(args, name)?;
            if b.is_empty || b.hi <= 0.0 {
                return Ok(Interval::empty());
            }
            let upper = b.hi.abs().max(b.lo.abs());
            Ok(Interval::new(0.0, upper))
        }
        _ => Err(MathError::Eval(format!(
            "interval eval: unsupported function `{}`",
            name
        ))),
    }
}

fn two(args: &[Interval], name: &str) -> Result<(Interval, Interval)> {
    if args.len() != 2 {
        return Err(MathError::InvalidArgument(format!("{} needs 2 args", name)));
    }
    Ok((args[0], args[1]))
}

// =========================================================================
// Tests
// =========================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn iv(lo: f64, hi: f64) -> Interval {
        Interval::new(lo, hi)
    }

    // --- Basic properties ---

    #[test]
    fn construct_and_contains() {
        let a = iv(1.0, 3.0);
        assert!(a.contains(1.0));
        assert!(a.contains(2.5));
        assert!(a.contains(3.0));
        assert!(!a.contains(0.999));
        assert!(!a.contains(3.001));
    }

    #[test]
    fn empty_interval() {
        let e = Interval::empty();
        assert!(e.is_empty);
        assert!(!e.contains(0.0));
        assert_eq!(e.to_string(), "∅");
    }

    #[test]
    fn inverted_bounds_become_empty() {
        let a = Interval::new(5.0, 1.0);
        assert!(a.is_empty);
    }

    #[test]
    fn width_and_midpoint() {
        let a = iv(2.0, 8.0);
        assert_abs_diff_eq!(a.width(), 6.0);
        assert_abs_diff_eq!(a.midpoint(), 5.0);
    }

    // --- Arithmetic ---

    #[test]
    fn add_sub() {
        let a = iv(1.0, 2.0);
        let b = iv(3.0, 4.0);
        assert_eq!(a + b, iv(4.0, 6.0));
        assert_eq!(a - b, iv(-3.0, -1.0));
    }

    #[test]
    fn mul_positive() {
        assert_eq!(iv(1.0, 2.0) * iv(3.0, 4.0), iv(3.0, 8.0));
    }

    #[test]
    fn mul_straddles_zero() {
        // [-1, 2] * [1, 3] → products: -1*1=-1, -1*3=-3, 2*1=2, 2*3=6 → [-3, 6]
        assert_eq!(iv(-1.0, 2.0) * iv(1.0, 3.0), iv(-3.0, 6.0));
    }

    #[test]
    fn div_no_zero() {
        assert_eq!(iv(1.0, 2.0) / iv(3.0, 4.0), iv(0.25, 2.0 / 3.0));
    }

    #[test]
    fn div_by_zero_interval_is_whole() {
        assert_eq!(iv(1.0, 2.0) / iv(-1.0, 1.0), Interval::whole());
    }

    #[test]
    fn neg_swaps_bounds() {
        assert_eq!(-iv(1.0, 3.0), iv(-3.0, -1.0));
    }

    // --- sqr and powi ---

    #[test]
    fn sqr_straddles_zero() {
        // [-2, 3]² → [0, 9]
        assert_eq!(iv(-2.0, 3.0).sqr(), iv(0.0, 9.0));
    }

    #[test]
    fn sqr_positive() {
        assert_eq!(iv(2.0, 3.0).sqr(), iv(4.0, 9.0));
    }

    #[test]
    fn sqr_negative() {
        assert_eq!(iv(-3.0, -1.0).sqr(), iv(1.0, 9.0));
    }

    #[test]
    fn powi_even_zero_crossing() {
        // [-2, 3]^4 → [0, 81]
        assert_eq!(iv(-2.0, 3.0).powi(4), iv(0.0, 81.0));
    }

    #[test]
    fn powi_odd_monotonic() {
        // [-2, 3]^3 → [-8, 27]
        assert_eq!(iv(-2.0, 3.0).powi(3), iv(-8.0, 27.0));
    }

    #[test]
    fn powi_zero_is_one() {
        assert_eq!(iv(-5.0, 5.0).powi(0), Interval::point(1.0));
    }

    // --- Elementary functions ---

    #[test]
    fn exp_monotonic() {
        let r = iv(0.0, 1.0).exp();
        assert_abs_diff_eq!(r.lo, 1.0);
        assert_abs_diff_eq!(r.hi, std::f64::consts::E);
    }

    #[test]
    fn ln_domain() {
        let r = iv(1.0, std::f64::consts::E).ln();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 1.0);
        // ln of negative interval → empty
        assert!(iv(-2.0, -1.0).ln().is_empty);
    }

    #[test]
    fn sqrt_domain() {
        let r = iv(0.0, 4.0).sqrt();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 2.0);
        assert!(iv(-4.0, -1.0).sqrt().is_empty);
    }

    #[test]
    fn abs_straddles_zero() {
        assert_eq!(iv(-3.0, 1.0).abs(), iv(0.0, 3.0));
    }

    #[test]
    fn sin_reaches_extrema() {
        // [0, 2π] → sin covers full [-1, 1]
        assert_eq!(iv(0.0, 2.0 * std::f64::consts::PI).sin(), iv(-1.0, 1.0));
        // [0, π/2] → [0, 1]
        let r = iv(0.0, std::f64::consts::FRAC_PI_2).sin();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 1.0);
    }

    #[test]
    fn cos_reaches_extrema() {
        // [0, π] → [-1, 1]
        assert_eq!(iv(0.0, std::f64::consts::PI).cos(), iv(-1.0, 1.0));
        // [-π/2, π/2] → [0, 1]
        let r = iv(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2).cos();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 1.0);
    }

    #[test]
    fn tan_pole_returns_whole() {
        // [0, π] contains π/2 → whole
        assert_eq!(iv(0.0, std::f64::consts::PI).tan(), Interval::whole());
        // [0, 1] (no pole) → [tan(0), tan(1)]
        let r = iv(0.0, 1.0).tan();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 1.0_f64.tan());
    }

    #[test]
    fn atan_bounded() {
        let r = iv(-1e10, 1e10).atan();
        assert!(r.lo < -1.5);
        assert!(r.hi > 1.5);
    }

    #[test]
    fn cosh_minimum_at_zero() {
        let r = iv(-1.0, 2.0).cosh();
        assert_abs_diff_eq!(r.lo, 1.0);
        assert_abs_diff_eq!(r.hi, 2.0_f64.cosh());
    }

    // --- Set operations ---

    #[test]
    fn intersect_and_hull() {
        let a = iv(1.0, 3.0);
        let b = iv(2.0, 5.0);
        assert_eq!(a.intersect(&b), iv(2.0, 3.0));
        assert_eq!(a.hull(&b), iv(1.0, 5.0));
    }

    #[test]
    fn intersect_disjoint_is_empty() {
        let a = iv(1.0, 2.0);
        let b = iv(3.0, 4.0);
        assert!(a.intersect(&b).is_empty);
    }

    // --- Expr evaluation ---

    #[test]
    fn eval_simple_poly() {
        use crate::parser::Parser;
        let e = Parser::parse("x^2 + 1").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-2.0, 3.0));
        let r = eval_interval(&e, &vars).unwrap();
        // x^2 ∈ [0, 9], +1 → [1, 10]
        assert_eq!(r, iv(1.0, 10.0));
    }

    #[test]
    fn eval_dependency_problem() {
        use crate::parser::Parser;
        // x - x over [1, 2] → [1-2, 2-1] = [-1, 1] (NOT [0,0] — dependency problem)
        let e = Parser::parse("x - x").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 2.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-1.0, 1.0));
    }

    #[test]
    fn eval_sin_range() {
        use crate::parser::Parser;
        let e = Parser::parse("sin(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(0.0, 2.0 * std::f64::consts::PI));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-1.0, 1.0));
    }

    #[test]
    fn eval_exp_and_const() {
        use crate::parser::Parser;
        let e = Parser::parse("exp(x) + pi").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(0.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, 1.0 + std::f64::consts::PI, epsilon = 1e-12);
        assert_abs_diff_eq!(r.hi, std::f64::consts::E + std::f64::consts::PI, epsilon = 1e-12);
    }

    #[test]
    fn eval_div_by_zero_interval() {
        use crate::parser::Parser;
        let e = Parser::parse("1 / x").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-1.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, Interval::whole());
    }

    #[test]
    fn eval_sqrt_of_negative_empty() {
        use crate::parser::Parser;
        let e = Parser::parse("sqrt(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-4.0, -1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert!(r.is_empty);
    }

    #[test]
    fn eval_integer_pow_tight() {
        use crate::parser::Parser;
        // x^3 over [-2, 1] → [-8, 1] (odd power, monotonic)
        let e = Parser::parse("x^3").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-2.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-8.0, 1.0));
    }

    #[test]
    fn eval_unknown_var_errors() {
        use crate::parser::Parser;
        let e = Parser::parse("y + 1").unwrap();
        let vars = HashMap::new();
        assert!(eval_interval(&e, &vars).is_err());
    }

    #[test]
    fn eval_nested_function() {
        use crate::parser::Parser;
        // sin(x^2) over [0, 2] → x^2 ∈ [0, 4], sin over [0, 4] (4 > π so hits max 1)
        let e = Parser::parse("sin(x^2)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(0.0, 2.0));
        let r = eval_interval(&e, &vars).unwrap();
        // sin on [0, 4]: 4 > π/2 so max=1; 4 < 3π/2 so no -1; min = min(sin(0), sin(4)) = sin(4) < 0
        assert_abs_diff_eq!(r.hi, 1.0);
        assert!(r.lo < 0.0);
    }

    #[test]
    fn display_formatting() {
        assert_eq!(iv(1.0, 2.0).to_string(), "[1, 2]");
        assert_eq!(Interval::whole().to_string(), "[-∞, ∞]");
        assert_eq!(iv(-1.5, 0.5).to_string(), "[-1.5, 0.5]");
    }

    // =================================================================
    // Additional coverage: arithmetic edge cases, domain handling,
    // function-specific behavior, Expr eval compositions, set ops,
    // and the fundamental containment property.
    // =================================================================

    // --- Arithmetic edge cases ---

    #[test]
    fn mul_both_negative() {
        // [-3, -1] * [-4, -2] → products: 12, 6, 4, 2 → [2, 12]
        assert_eq!(iv(-3.0, -1.0) * iv(-4.0, -2.0), iv(2.0, 12.0));
    }

    #[test]
    fn mul_both_straddle_zero() {
        // [-1, 1] * [-2, 2] → products: 2, -2, -2, 2 → [-2, 2]
        assert_eq!(iv(-1.0, 1.0) * iv(-2.0, 2.0), iv(-2.0, 2.0));
    }

    #[test]
    fn sub_negative_intervals() {
        // [-3, -1] - [-4, -2] = [-3 - (-2), -1 - (-4)] = [-1, 3]
        assert_eq!(iv(-3.0, -1.0) - iv(-4.0, -2.0), iv(-1.0, 3.0));
    }

    #[test]
    fn div_negative_divisor() {
        // [1, 2] / [-4, -3] → [1/-3, 2/-4] = [-1/3, -1/2]... but -1/3 > -1/2
        // products: 1/-4=-0.25, 1/-3=-0.333, 2/-4=-0.5, 2/-3=-0.666
        // lo = -0.666, hi = -0.25
        let r = iv(1.0, 2.0) / iv(-4.0, -3.0);
        assert_abs_diff_eq!(r.lo, 2.0 / -3.0);
        assert_abs_diff_eq!(r.hi, 1.0 / -4.0);
    }

    #[test]
    fn div_by_exact_zero_is_empty() {
        let r = iv(1.0, 2.0) / Interval::point(0.0);
        assert!(r.is_empty);
    }

    #[test]
    fn empty_propagates_through_arithmetic() {
        let e = Interval::empty();
        assert!((e + iv(1.0, 2.0)).is_empty);
        assert!((e - iv(1.0, 2.0)).is_empty);
        assert!((e * iv(1.0, 2.0)).is_empty);
        assert!((e / iv(1.0, 2.0)).is_empty);
        assert!((-e).is_empty);
    }

    #[test]
    fn point_interval_arithmetic() {
        let p = Interval::point(3.0);
        assert_eq!(p + Interval::point(4.0), Interval::point(7.0));
        assert_eq!(p * Interval::point(2.0), Interval::point(6.0));
        assert_eq!(p - Interval::point(1.0), Interval::point(2.0));
    }

    #[test]
    fn whole_interval_arithmetic() {
        let w = Interval::whole();
        // whole + anything = whole
        assert_eq!(w + iv(1.0, 2.0), Interval::whole());
        // whole * non-point = whole (conservatively)
        assert_eq!(w * iv(1.0, 2.0), Interval::whole());
    }

    // --- powi additional cases ---

    #[test]
    fn powi_two_no_zero_crossing() {
        // [2, 3]^2 → [4, 9] (even, no zero crossing)
        assert_eq!(iv(2.0, 3.0).powi(2), iv(4.0, 9.0));
    }

    #[test]
    fn powi_one_is_identity() {
        assert_eq!(iv(-5.0, 3.0).powi(1), iv(-5.0, 3.0));
    }

    #[test]
    fn powi_even_negative_only() {
        // [-3, -1]^2 → [1, 9] (even, no zero crossing, negative interval)
        assert_eq!(iv(-3.0, -1.0).powi(2), iv(1.0, 9.0));
    }

    #[test]
    fn powi_large_even() {
        // [-1, 2]^6 → [0, 64]
        assert_eq!(iv(-1.0, 2.0).powi(6), iv(0.0, 64.0));
    }

    // --- Function domain handling ---

    #[test]
    fn log10_domain() {
        let r = iv(1.0, 100.0).log10();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 2.0);
        assert!(iv(-1.0, 0.0).log10().is_empty);
    }

    #[test]
    fn log2_domain() {
        let r = iv(1.0, 8.0).log2();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, 3.0);
        assert!(iv(-1.0, 0.0).log2().is_empty);
    }

    #[test]
    fn ln_straddling_zero() {
        // ln([-1, 2]) → ln of [0, 2] (negative part excluded) → [-inf, ln(2)]
        let r = iv(-1.0, 2.0).ln();
        assert!(r.lo.is_infinite() && r.lo < 0.0);
        assert_abs_diff_eq!(r.hi, 2.0_f64.ln());
    }

    #[test]
    fn asin_in_domain() {
        let r = iv(-0.5, 0.5).asin();
        assert_abs_diff_eq!(r.lo, (-0.5_f64).asin());
        assert_abs_diff_eq!(r.hi, 0.5_f64.asin());
    }

    #[test]
    fn asin_outside_domain() {
        // [2, 3] is entirely outside [-1, 1] → empty
        assert!(iv(2.0, 3.0).asin().is_empty);
    }

    #[test]
    fn asin_straddling_domain() {
        // [-2, 0.5] → clamped to [-1, 0.5] → [asin(-1), asin(0.5)]
        let r = iv(-2.0, 0.5).asin();
        assert_abs_diff_eq!(r.lo, (-1.0_f64).asin());
        assert_abs_diff_eq!(r.hi, 0.5_f64.asin());
    }

    #[test]
    fn acos_decreasing_swaps_bounds() {
        // acos is decreasing: [0, 1] → [acos(1), acos(0)] = [0, π/2]
        let r = iv(0.0, 1.0).acos();
        assert_abs_diff_eq!(r.lo, 0.0);
        assert_abs_diff_eq!(r.hi, std::f64::consts::FRAC_PI_2);
    }

    #[test]
    fn acos_outside_domain() {
        assert!(iv(2.0, 3.0).acos().is_empty);
    }

    #[test]
    fn tanh_bounded() {
        let r = iv(-100.0, 100.0).tanh();
        // tanh saturates to ±1 for large |x|; bounds are [-1, 1] inclusive
        assert!(r.lo >= -1.0);
        assert!(r.hi <= 1.0);
        assert!(r.lo < -0.99);
        assert!(r.hi > 0.99);
    }

    #[test]
    fn sinh_monotonic() {
        let r = iv(-1.0, 1.0).sinh();
        assert_abs_diff_eq!(r.lo, (-1.0_f64).sinh());
        assert_abs_diff_eq!(r.hi, 1.0_f64.sinh());
    }

    #[test]
    fn cosh_no_zero_crossing() {
        // [2, 3] — cosh is increasing for x > 0
        let r = iv(2.0, 3.0).cosh();
        assert_abs_diff_eq!(r.lo, 2.0_f64.cosh());
        assert_abs_diff_eq!(r.hi, 3.0_f64.cosh());
    }

    #[test]
    fn cbrt_negative_domain() {
        // cbrt is defined for all reals, monotonically increasing
        let r = iv(-8.0, -1.0);
        let c = r.cbrt();
        assert_abs_diff_eq!(c.lo, -2.0);
        assert_abs_diff_eq!(c.hi, -1.0);
    }

    #[test]
    fn cbrt_straddles_zero() {
        let r = iv(-8.0, 27.0).cbrt();
        assert_abs_diff_eq!(r.lo, -2.0);
        assert_abs_diff_eq!(r.hi, 3.0);
    }

    #[test]
    fn abs_positive_interval() {
        assert_eq!(iv(2.0, 5.0).abs(), iv(2.0, 5.0));
    }

    #[test]
    fn abs_negative_interval() {
        assert_eq!(iv(-5.0, -2.0).abs(), iv(2.0, 5.0));
    }

    // --- Set operations additional ---

    #[test]
    fn overlaps_basic() {
        assert!(iv(1.0, 3.0).overlaps(&iv(2.0, 5.0)));
        assert!(!iv(1.0, 2.0).overlaps(&iv(3.0, 4.0)));
        // Touching at a point counts as overlap
        assert!(iv(1.0, 2.0).overlaps(&iv(2.0, 3.0)));
    }

    #[test]
    fn overlaps_with_empty() {
        assert!(!Interval::empty().overlaps(&iv(1.0, 2.0)));
        assert!(!iv(1.0, 2.0).overlaps(&Interval::empty()));
    }

    #[test]
    fn hull_with_empty() {
        let a = iv(1.0, 3.0);
        assert_eq!(a.hull(&Interval::empty()), a);
        assert_eq!(Interval::empty().hull(&a), a);
    }

    #[test]
    fn width_midpoint_of_empty() {
        assert!(Interval::empty().width().is_nan());
        assert!(Interval::empty().midpoint().is_nan());
    }

    #[test]
    fn contains_nan_is_false() {
        let a = iv(1.0, 3.0);
        assert!(!a.contains(f64::NAN));
    }

    #[test]
    fn contains_infinity() {
        let w = Interval::whole();
        assert!(w.contains(f64::INFINITY));
        assert!(w.contains(f64::NEG_INFINITY));
    }

    // --- Display additional ---

    #[test]
    fn display_point_interval() {
        assert_eq!(Interval::point(5.0).to_string(), "[5, 5]");
    }

    #[test]
    fn display_half_infinite() {
        assert_eq!(iv(0.0, f64::INFINITY).to_string(), "[0, ∞]");
        assert_eq!(iv(f64::NEG_INFINITY, 0.0).to_string(), "[-∞, 0]");
    }

    // --- Expr eval: compositions and domain ---

    #[test]
    fn eval_multivar_complex() {
        use crate::parser::Parser;
        // (x + y) * (x - y) over x∈[1,3], y∈[2,4]
        // x+y ∈ [3, 7], x-y ∈ [-3, 1], product ∈ [-21, 7]
        let e = Parser::parse("(x + y) * (x - y)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 3.0));
        vars.insert("y".to_string(), iv(2.0, 4.0));
        let r = eval_interval(&e, &vars).unwrap();
        // Conservative: x+y ∈ [3,7], x-y ∈ [-3,1], so product ∈ [-21, 7]
        assert!(r.contains(-21.0));
        assert!(r.contains(7.0));
    }

    #[test]
    fn eval_constant_expression() {
        use crate::parser::Parser;
        // No variables — pure constant
        let e = Parser::parse("2 + 3 * 4").unwrap();
        let vars = HashMap::new();
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, Interval::point(14.0));
    }

    #[test]
    fn eval_neg_node() {
        use crate::parser::Parser;
        let e = Parser::parse("-x").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 3.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-3.0, -1.0));
    }

    #[test]
    fn eval_division_expression() {
        use crate::parser::Parser;
        let e = Parser::parse("1 / (x + 1)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 3.0));
        // x+1 ∈ [2, 4], 1/[2,4] = [1/4, 1/2]
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, 0.25);
        assert_abs_diff_eq!(r.hi, 0.5);
    }

    #[test]
    fn eval_constants_pi_e_tau() {
        use crate::parser::Parser;
        for (src, expected) in [
            ("pi", std::f64::consts::PI),
            ("e", std::f64::consts::E),
            ("tau", std::f64::consts::TAU),
        ] {
            let e = Parser::parse(src).unwrap();
            let vars = HashMap::new();
            let r = eval_interval(&e, &vars).unwrap();
            assert!(r.is_empty == false);
            assert_abs_diff_eq!(r.lo, expected);
            assert_abs_diff_eq!(r.hi, expected);
        }
    }

    #[test]
    fn eval_unsupported_function_errors() {
        use crate::parser::Parser;
        // bessel_j0 is not in the interval function set
        let e = Parser::parse("bessel_j0(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(0.0, 1.0));
        assert!(eval_interval(&e, &vars).is_err());
    }

    #[test]
    fn eval_pow_negative_base_integer_exp() {
        use crate::parser::Parser;
        // x^2 over [-3, 1] → [0, 9] (even power, zero crossing)
        let e = Parser::parse("x^2").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-3.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(0.0, 9.0));
    }

    #[test]
    fn eval_pow_negative_base_noninteger_errors() {
        use crate::parser::Parser;
        // x^0.5 over [-1, 1] — non-integer exponent with negative base → error
        let e = Parser::parse("x^0.5").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-1.0, 1.0));
        assert!(eval_interval(&e, &vars).is_err());
    }

    #[test]
    fn eval_pow_positive_base_noninteger() {
        use crate::parser::Parser;
        // x^0.5 over [1, 4] → [1, 2]
        let e = Parser::parse("x^0.5").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 4.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, 1.0);
        assert_abs_diff_eq!(r.hi, 2.0);
    }

    #[test]
    fn eval_func_min_max() {
        use crate::parser::Parser;
        let e = Parser::parse("min(x, y)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 5.0));
        vars.insert("y".to_string(), iv(2.0, 3.0));
        let r = eval_interval(&e, &vars).unwrap();
        // min of [1,5] and [2,3] → [min(1,2), min(5,3)] = [1, 3]
        assert_eq!(r, iv(1.0, 3.0));
    }

    #[test]
    fn eval_func_max() {
        use crate::parser::Parser;
        let e = Parser::parse("max(x, y)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.0, 5.0));
        vars.insert("y".to_string(), iv(2.0, 3.0));
        let r = eval_interval(&e, &vars).unwrap();
        // max of [1,5] and [2,3] → [max(1,2), max(5,3)] = [2, 5]
        assert_eq!(r, iv(2.0, 5.0));
    }

    #[test]
    fn eval_func_abs() {
        use crate::parser::Parser;
        let e = Parser::parse("abs(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-3.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(0.0, 3.0));
    }

    #[test]
    fn eval_func_floor_ceil() {
        use crate::parser::Parser;
        let e = Parser::parse("floor(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(1.5, 3.7));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(1.0, 3.0));

        let e = Parser::parse("ceil(x)").unwrap();
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(2.0, 4.0));
    }

    #[test]
    fn eval_func_sign() {
        use crate::parser::Parser;
        let e = Parser::parse("sign(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-1.0, 1.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-1.0, 1.0));

        let mut vars2 = HashMap::new();
        vars2.insert("x".to_string(), iv(2.0, 5.0));
        let r2 = eval_interval(&e, &vars2).unwrap();
        assert_eq!(r2, Interval::point(1.0));
    }

    #[test]
    fn eval_func_log10_log2() {
        use crate::parser::Parser;
        let e = Parser::parse("log10(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(10.0, 1000.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, 1.0);
        assert_abs_diff_eq!(r.hi, 3.0);

        let e = Parser::parse("log2(x)").unwrap();
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, 10.0_f64.log2());
        assert_abs_diff_eq!(r.hi, 1000.0_f64.log2());
    }

    #[test]
    fn eval_func_cbrt() {
        use crate::parser::Parser;
        let e = Parser::parse("cbrt(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-8.0, 27.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert_abs_diff_eq!(r.lo, -2.0);
        assert_abs_diff_eq!(r.hi, 3.0);
    }

    #[test]
    fn eval_func_tanh() {
        use crate::parser::Parser;
        let e = Parser::parse("tanh(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(-100.0, 100.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert!(r.lo >= -1.0 && r.hi <= 1.0);
    }

    #[test]
    fn eval_nested_arithmetic() {
        use crate::parser::Parser;
        // (a + b) * (c - d) over a∈[1,2], b∈[3,4], c∈[5,6], d∈[1,2]
        // a+b ∈ [4,6], c-d ∈ [3,5], product ∈ [12, 30]
        let e = Parser::parse("(a + b) * (c - d)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("a".to_string(), iv(1.0, 2.0));
        vars.insert("b".to_string(), iv(3.0, 4.0));
        vars.insert("c".to_string(), iv(5.0, 6.0));
        vars.insert("d".to_string(), iv(1.0, 2.0));
        let r = eval_interval(&e, &vars).unwrap();
        assert!(r.contains(12.0));
        assert!(r.contains(30.0));
    }

    #[test]
    fn eval_composition_sin_plus_cos() {
        use crate::parser::Parser;
        // sin(x) + cos(x) over [0, 2π] → each is [-1,1], sum is [-2, 2]
        // (true range is [-√2, √2] but interval arithmetic is conservative)
        let e = Parser::parse("sin(x) + cos(x)").unwrap();
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), iv(0.0, 2.0 * std::f64::consts::PI));
        let r = eval_interval(&e, &vars).unwrap();
        assert_eq!(r, iv(-2.0, 2.0));
    }

    // --- Containment property (fundamental correctness) ---

    #[test]
    fn containment_property_polynomial() {
        use crate::parser::Parser;
        // For any x in [1, 3], x^2 + 2*x must lie in the interval result
        let e = Parser::parse("x^2 + 2*x").unwrap();
        let x_iv = iv(1.0, 3.0);
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), x_iv);
        let bounds = eval_interval(&e, &vars).unwrap();
        // Sample many points in the interval and verify containment
        let n = 100;
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = x_iv.lo + t * x_iv.width();
            let mut v = HashMap::new();
            v.insert("x".to_string(), Interval::point(x));
            let val = eval_interval(&e, &v).unwrap();
            assert!(
                bounds.contains(val.lo) && bounds.contains(val.hi),
                "point x={} gave value {} not in bounds {}",
                x, val, bounds
            );
        }
    }

    #[test]
    fn containment_property_trig() {
        use crate::parser::Parser;
        let e = Parser::parse("sin(x) + cos(x)").unwrap();
        let x_iv = iv(0.0, std::f64::consts::PI);
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), x_iv);
        let bounds = eval_interval(&e, &vars).unwrap();
        let n = 100;
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = x_iv.lo + t * x_iv.width();
            let mut v = HashMap::new();
            v.insert("x".to_string(), Interval::point(x));
            let val = eval_interval(&e, &v).unwrap();
            assert!(
                bounds.contains(val.lo) && bounds.contains(val.hi),
                "point x={} gave value {} not in bounds {}",
                x, val, bounds
            );
        }
    }

    #[test]
    fn containment_property_exp_log() {
        use crate::parser::Parser;
        let e = Parser::parse("exp(x) + ln(x + 1)").unwrap();
        let x_iv = iv(0.5, 2.0);
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), x_iv);
        let bounds = eval_interval(&e, &vars).unwrap();
        let n = 50;
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = x_iv.lo + t * x_iv.width();
            let mut v = HashMap::new();
            v.insert("x".to_string(), Interval::point(x));
            let val = eval_interval(&e, &v).unwrap();
            assert!(
                bounds.contains(val.lo) && bounds.contains(val.hi),
                "point x={} gave value {} not in bounds {}",
                x, val, bounds
            );
        }
    }

    #[test]
    fn containment_property_rational() {
        use crate::parser::Parser;
        // 1 / (1 + x^2) over [-2, 2]
        let e = Parser::parse("1 / (1 + x^2)").unwrap();
        let x_iv = iv(-2.0, 2.0);
        let mut vars = HashMap::new();
        vars.insert("x".to_string(), x_iv);
        let bounds = eval_interval(&e, &vars).unwrap();
        let n = 100;
        for i in 0..=n {
            let t = i as f64 / n as f64;
            let x = x_iv.lo + t * x_iv.width();
            let mut v = HashMap::new();
            v.insert("x".to_string(), Interval::point(x));
            let val = eval_interval(&e, &v).unwrap();
            assert!(
                bounds.contains(val.lo) && bounds.contains(val.hi),
                "point x={} gave value {} not in bounds {}",
                x, val, bounds
            );
        }
    }
}
