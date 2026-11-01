//! Number-base conversion.
//!
//! Converts an integer to a representation in another base — binary,
//! octal, hexadecimal, or any base from 2 to 36.
//!
//! The conversion is **exact**: the value is carried in a [`BigInt`], so a
//! 40-digit decimal converts to binary without loss. Fractional parts are
//! not supported — this is an *integer* base conversion.
//!
//! Cross-port note: this fills a gap identified in the qalculate competitive
//! pass (number bases), and complements the `is-prime`/`factor` BigInt
//! auto-upgrade.

use crate::error::{MathError, Result};
use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive, Zero};

/// The digits used for bases above 10: `A`–`Z` (and lowercase is accepted
/// on input).
const DIGITS: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// The smallest supported base.
pub const MIN_BASE: u32 = 2;
/// The largest supported base.
pub const MAX_BASE: u32 = 36;

/// Parse a base name or number into a base.
///
/// Accepts the names `bin`/`binary` (2), `oct`/`octal` (8), `hex`/`hexadecimal`
/// (16), and any number from 2 to 36.
pub fn parse_base(name: &str) -> Result<u32> {
    let n = name.trim().to_ascii_lowercase();
    let base = match n.as_str() {
        "bin" | "binary" => 2,
        "oct" | "octal" => 8,
        "dec" | "decimal" => 10,
        "hex" | "hexadecimal" => 16,
        other => {
            let b: u32 = other.parse().map_err(|_| {
                MathError::Eval(format!(
                    "unknown base '{}' (use bin, oct, dec, hex, or a number from 2 to 36)",
                    name
                ))
            })?;
            if !(MIN_BASE..=MAX_BASE).contains(&b) {
                return Err(MathError::Eval(format!(
                    "base {} is out of range ({} to {})",
                    b, MIN_BASE, MAX_BASE
                )));
            }
            b
        }
    };
    Ok(base)
}

/// Convert `n` to a string in the given base.
///
/// The value is signed: a negative number gets a leading `-`.
pub fn to_base(n: &BigInt, base: u32) -> Result<String> {
    if !(MIN_BASE..=MAX_BASE).contains(&base) {
        return Err(MathError::Eval(format!(
            "base {} is out of range ({} to {})",
            base, MIN_BASE, MAX_BASE
        )));
    }
    if n.is_zero() {
        return Ok("0".into());
    }
    let negative = n < &BigInt::zero();
    let mut magnitude = n.abs();
    let b = BigInt::from(base);
    let mut digits: Vec<u8> = Vec::new();
    while !magnitude.is_zero() {
        let rem = (&magnitude % &b).to_u32().unwrap_or(0) as usize;
        digits.push(DIGITS[rem]);
        magnitude /= &b;
    }
    if negative {
        digits.push(b'-');
    }
    digits.reverse();
    Ok(String::from_utf8(digits).expect("digits are ASCII"))
}

/// Parse a string written in `base` into a [`BigInt`].
///
/// Letters are accepted in either case (`ff`, `FF`).
pub fn from_base(s: &str, base: u32) -> Result<BigInt> {
    if !(MIN_BASE..=MAX_BASE).contains(&base) {
        return Err(MathError::Eval(format!(
            "base {} is out of range ({} to {})",
            base, MIN_BASE, MAX_BASE
        )));
    }
    let t = s.trim();
    let (negative, digits) = match t.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, t.strip_prefix('+').unwrap_or(t)),
    };
    if digits.is_empty() {
        return Err(MathError::Eval(format!("no digits in '{}'", s)));
    }
    let mut value = BigInt::zero();
    let b = BigInt::from(base);
    for c in digits.chars() {
        if c == '_' {
            // Digit separators are ignored (as in most languages).
            continue;
        }
        let d = c
            .to_digit(36)
            .ok_or_else(|| MathError::Eval(format!("invalid digit '{}' in base {}", c, base)))?;
        if d >= base {
            return Err(MathError::Eval(format!(
                "digit '{}' is not valid in base {}",
                c, base
            )));
        }
        value = value * &b + BigInt::from(d);
    }
    if negative {
        value = -value;
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(s: &str) -> BigInt {
        s.parse().unwrap()
    }

    #[test]
    fn to_base_named_bases() {
        assert_eq!(to_base(&n("255"), 2).unwrap(), "11111111");
        assert_eq!(to_base(&n("255"), 8).unwrap(), "377");
        assert_eq!(to_base(&n("255"), 16).unwrap(), "FF");
        assert_eq!(to_base(&n("255"), 10).unwrap(), "255");
    }

    #[test]
    fn to_base_zero_and_one() {
        assert_eq!(to_base(&n("0"), 16).unwrap(), "0");
        assert_eq!(to_base(&n("1"), 16).unwrap(), "1");
    }

    #[test]
    fn to_base_negative() {
        assert_eq!(to_base(&n("-255"), 16).unwrap(), "-FF");
        assert_eq!(to_base(&n("-10"), 2).unwrap(), "-1010");
    }

    #[test]
    fn to_base_is_exact_for_big_values() {
        // 2^100 in binary.
        let v = BigInt::from(1u8) << 100usize;
        let s = to_base(&v, 2).unwrap();
        assert_eq!(s.len(), 101, "got {}", s);
        assert!(s.starts_with('1'), "got {}", s);
        // Round-trip.
        assert_eq!(from_base(&s, 2).unwrap(), v);
    }

    #[test]
    fn to_base_base36() {
        assert_eq!(to_base(&n("35"), 36).unwrap(), "Z");
        assert_eq!(to_base(&n("36"), 36).unwrap(), "10");
    }

    #[test]
    fn to_base_rejects_bad_base() {
        assert!(to_base(&n("10"), 1).is_err());
        assert!(to_base(&n("10"), 37).is_err());
        assert!(to_base(&n("10"), 0).is_err());
    }

    #[test]
    fn parse_base_names() {
        assert_eq!(parse_base("bin").unwrap(), 2);
        assert_eq!(parse_base("binary").unwrap(), 2);
        assert_eq!(parse_base("oct").unwrap(), 8);
        assert_eq!(parse_base("octal").unwrap(), 8);
        assert_eq!(parse_base("dec").unwrap(), 10);
        assert_eq!(parse_base("hex").unwrap(), 16);
        assert_eq!(parse_base("hexadecimal").unwrap(), 16);
        assert_eq!(parse_base("HEX").unwrap(), 16);
        assert_eq!(parse_base("7").unwrap(), 7);
        assert_eq!(parse_base("36").unwrap(), 36);
    }

    #[test]
    fn parse_base_rejects_bad_input() {
        assert!(parse_base("1").is_err());
        assert!(parse_base("37").is_err());
        assert!(parse_base("nonsense").is_err());
    }

    #[test]
    fn from_base_parses_digits() {
        assert_eq!(from_base("11111111", 2).unwrap(), n("255"));
        assert_eq!(from_base("377", 8).unwrap(), n("255"));
        assert_eq!(from_base("FF", 16).unwrap(), n("255"));
        assert_eq!(from_base("ff", 16).unwrap(), n("255"));
        assert_eq!(from_base("-10", 2).unwrap(), n("-2"));
    }

    #[test]
    fn from_base_rejects_invalid_digits() {
        // '2' is not a binary digit.
        assert!(from_base("12", 2).is_err());
        // '9' is not an octal digit.
        assert!(from_base("9", 8).is_err());
        assert!(from_base("", 16).is_err());
        assert!(from_base("g", 16).is_err());
    }

    #[test]
    fn round_trip_across_bases() {
        for v in [0i64, 1, 7, 255, 65535, 1_000_000_007] {
            for base in 2u32..=36 {
                let s = to_base(&n(&v.to_string()), base).unwrap();
                assert_eq!(from_base(&s, base).unwrap(), n(&v.to_string()), "v={} base={}", v, base);
            }
        }
    }
}
