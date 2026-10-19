//! Euclidean vector primitives (dot and 3-D cross product).

use crate::error::{MathError, Result};

/// Dot product of two equal-length vectors.
pub fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() {
        return Err(MathError::InvalidArgument(format!(
            "dot: length mismatch ({} vs {})",
            a.len(),
            b.len()
        )));
    }
    Ok(a.iter().zip(b.iter()).map(|(x, y)| x * y).sum())
}

/// Cross product of two 3-vectors.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Cross product; both slices must have length 3.
pub fn cross(a: &[f64], b: &[f64]) -> Result<[f64; 3]> {
    if a.len() != 3 || b.len() != 3 {
        return Err(MathError::InvalidArgument(
            "cross: requires two 3-vectors (use `a1 a2 a3 | b1 b2 b3`)".into(),
        ));
    }
    Ok(cross3([a[0], a[1], a[2]], [b[0], b[1], b[2]]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn dot_basic() {
        assert_relative_eq!(dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap(), 32.0);
    }

    #[test]
    fn cross_basic() {
        let c = cross(&[1.0, 0.0, 0.0], &[0.0, 1.0, 0.0]).unwrap();
        assert_relative_eq!(c[0], 0.0);
        assert_relative_eq!(c[1], 0.0);
        assert_relative_eq!(c[2], 1.0);
    }

    #[test]
    fn dot_length_mismatch() {
        assert!(dot(&[1.0], &[1.0, 2.0]).is_err());
    }
}
