//! Matrix operations from scratch.
//!
//! Provides a [`Matrix`] type with arithmetic (add, sub, mul, scalar),
//! transpose, determinant via Gaussian elimination, inverse, and
//! linear system solving (Ax = b).

use crate::error::{MathError, Result};
use std::fmt;
use std::ops::{Add, Mul, Sub};

/// A row-major `f64` matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    data: Vec<f64>,
}

impl Matrix {
    /// Create a matrix from a flat row-major vector.
    pub fn from_row_major(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self> {
        if data.len() != rows * cols {
            return Err(MathError::InvalidArgument(format!(
                "from_row_major: expected {} elements, got {}",
                rows * cols,
                data.len()
            )));
        }
        Ok(Self { rows, cols, data })
    }

    /// Create a matrix from nested rows.
    pub fn from_rows(rows_data: &[Vec<f64>]) -> Result<Self> {
        if rows_data.is_empty() {
            return Err(MathError::InvalidArgument("from_rows: empty".into()));
        }
        let rows = rows_data.len();
        let cols = rows_data[0].len();
        if rows_data.iter().any(|r| r.len() != cols) {
            return Err(MathError::InvalidArgument("from_rows: ragged rows".into()));
        }
        let mut data = Vec::with_capacity(rows * cols);
        for r in rows_data {
            data.extend_from_slice(r);
        }
        Ok(Self { rows, cols, data })
    }

    /// Create an `rows × cols` matrix filled with `fill`.
    pub fn filled(rows: usize, cols: usize, fill: f64) -> Self {
        Self { rows, cols, data: vec![fill; rows * cols] }
    }

    /// Create an `n × n` identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::filled(n, n, 0.0);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        m
    }

    /// Create an `n × n` zero matrix.
    pub fn zeros(n: usize) -> Self {
        Self::filled(n, n, 0.0)
    }

    /// Element access (row, col).
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.cols + col]
    }

    /// Mutable element access (row, col).
    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        self.data[row * self.cols + col] = val;
    }

    /// Return a row as a slice.
    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    /// Return a column as a vector.
    pub fn col(&self, j: usize) -> Vec<f64> {
        (0..self.rows).map(|i| self[(i, j)]).collect()
    }

    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut data = vec![0.0; self.rows * self.cols];
        for i in 0..self.rows {
            for j in 0..self.cols {
                data[j * self.rows + i] = self[(i, j)];
            }
        }
        Self { rows: self.cols, cols: self.rows, data }
    }

    /// Check if the matrix is square.
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    /// Determinant via Gaussian elimination with partial pivoting.
    pub fn determinant(&self) -> Result<f64> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument(
                "determinant: matrix must be square".into(),
            ));
        }
        let n = self.rows;
        let mut a = self.data.clone();
        let mut det = 1.0;
        for col in 0..n {
            // partial pivot
            let mut pivot = col;
            let mut best = a[col * n + col].abs();
            for r in (col + 1)..n {
                let val = a[r * n + col].abs();
                if val > best {
                    best = val;
                    pivot = r;
                }
            }
            if best < 1e-14 {
                return Ok(0.0);
            }
            if pivot != col {
                a.swap_within(col * n..(col + 1) * n, pivot * n..(pivot + 1) * n);
                det = -det;
            }
            let piv = a[col * n + col];
            det *= piv;
            for r in (col + 1)..n {
                let factor = a[r * n + col] / piv;
                for c in col..n {
                    a[r * n + c] -= factor * a[col * n + c];
                }
            }
        }
        Ok(det)
    }

    /// Solve the linear system `self * x = b` via Gaussian elimination with
    /// partial pivoting. Returns the solution vector `x`.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument("solve: matrix must be square".into()));
        }
        let n = self.rows;
        if b.len() != n {
            return Err(MathError::InvalidArgument(format!(
                "solve: b has length {}, expected {}",
                b.len(),
                n
            )));
        }
        // augmented matrix
        let mut a = vec![0.0; n * (n + 1)];
        for i in 0..n {
            for j in 0..n {
                a[i * (n + 1) + j] = self[(i, j)];
            }
            a[i * (n + 1) + n] = b[i];
        }
        // forward elimination with partial pivoting
        for col in 0..n {
            let mut pivot = col;
            let mut best = a[col * (n + 1) + col].abs();
            for r in (col + 1)..n {
                let val = a[r * (n + 1) + col].abs();
                if val > best {
                    best = val;
                    pivot = r;
                }
            }
            if best < 1e-14 {
                return Err(MathError::InvalidArgument("solve: singular matrix".into()));
            }
            if pivot != col {
                for c in 0..=n {
                    let tmp = a[col * (n + 1) + c];
                    a[col * (n + 1) + c] = a[pivot * (n + 1) + c];
                    a[pivot * (n + 1) + c] = tmp;
                }
            }
            let piv = a[col * (n + 1) + col];
            for r in (col + 1)..n {
                let factor = a[r * (n + 1) + col] / piv;
                for c in col..=n {
                    a[r * (n + 1) + c] -= factor * a[col * (n + 1) + c];
                }
            }
        }
        // back substitution
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut sum = a[i * (n + 1) + n];
            for j in (i + 1)..n {
                sum -= a[i * (n + 1) + j] * x[j];
            }
            x[i] = sum / a[i * (n + 1) + i];
        }
        Ok(x)
    }

    /// Compute the inverse via Gaussian elimination on the augmented [A | I] matrix.
    /// Performs a single elimination pass instead of n separate solves.
    pub fn inverse(&self) -> Result<Matrix> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument("inverse: matrix must be square".into()));
        }
        let n = self.rows;
        // Build augmented matrix [A | I] with 2n columns
        let mut a = vec![0.0; n * 2 * n];
        for i in 0..n {
            for j in 0..n {
                a[i * 2 * n + j] = self[(i, j)];
            }
            a[i * 2 * n + n + i] = 1.0;
        }
        // Forward elimination with partial pivoting
        for col in 0..n {
            let mut pivot = col;
            let mut best = a[col * 2 * n + col].abs();
            for r in (col + 1)..n {
                let val = a[r * 2 * n + col].abs();
                if val > best {
                    best = val;
                    pivot = r;
                }
            }
            if best < 1e-14 {
                return Err(MathError::InvalidArgument("inverse: singular matrix".into()));
            }
            if pivot != col {
                for c in 0..(2 * n) {
                    a.swap(col * 2 * n + c, pivot * 2 * n + c);
                }
            }
            let piv = a[col * 2 * n + col];
            for r in (col + 1)..n {
                let factor = a[r * 2 * n + col] / piv;
                if factor != 0.0 {
                    for c in col..(2 * n) {
                        a[r * 2 * n + c] -= factor * a[col * 2 * n + c];
                    }
                }
            }
        }
        // Back substitution for all n RHS columns at once
        for col in 0..n {
            let piv = a[col * 2 * n + col];
            for c in col..(2 * n) {
                a[col * 2 * n + c] /= piv;
            }
            for i in (0..col).rev() {
                let factor = a[i * 2 * n + col];
                if factor != 0.0 {
                    for c in col..(2 * n) {
                        a[i * 2 * n + c] -= factor * a[col * 2 * n + c];
                    }
                }
            }
        }
        // Extract the right half (inverse)
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                data[i * n + j] = a[i * 2 * n + n + j];
            }
        }
        Ok(Self { rows: n, cols: n, data })
    }

    /// Matrix-vector product `self * v`.
    pub fn mul_vec(&self, v: &[f64]) -> Result<Vec<f64>> {
        if v.len() != self.cols {
            return Err(MathError::InvalidArgument(format!(
                "mul_vec: vector length {} != cols {}",
                v.len(),
                self.cols
            )));
        }
        let mut out = vec![0.0; self.rows];
        for i in 0..self.rows {
            for j in 0..self.cols {
                out[i] += self[(i, j)] * v[j];
            }
        }
        Ok(out)
    }

    /// Trace (sum of diagonal).
    pub fn trace(&self) -> Result<f64> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument("trace: matrix must be square".into()));
        }
        Ok((0..self.rows).map(|i| self[(i, i)]).sum())
    }
}

// --- Indexing ---------------------------------------------------------------

impl std::ops::Index<(usize, usize)> for Matrix {
    type Output = f64;
    fn index(&self, (r, c): (usize, usize)) -> &f64 {
        &self.data[r * self.cols + c]
    }
}

impl std::ops::IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut f64 {
        &mut self.data[r * self.cols + c]
    }
}

// --- Arithmetic -------------------------------------------------------------

impl Add for &Matrix {
    type Output = Result<Matrix>;
    fn add(self, rhs: &Matrix) -> Result<Matrix> {
        if self.rows != rhs.rows || self.cols != rhs.cols {
            return Err(MathError::InvalidArgument(format!(
                "add: shape mismatch {}x{} vs {}x{}",
                self.rows, self.cols, rhs.rows, rhs.cols
            )));
        }
        let data: Vec<f64> = self.data.iter().zip(rhs.data.iter()).map(|(a, b)| a + b).collect();
        Ok(Matrix { rows: self.rows, cols: self.cols, data })
    }
}

impl Sub for &Matrix {
    type Output = Result<Matrix>;
    fn sub(self, rhs: &Matrix) -> Result<Matrix> {
        if self.rows != rhs.rows || self.cols != rhs.cols {
            return Err(MathError::InvalidArgument(format!(
                "sub: shape mismatch {}x{} vs {}x{}",
                self.rows, self.cols, rhs.rows, rhs.cols
            )));
        }
        let data: Vec<f64> = self.data.iter().zip(rhs.data.iter()).map(|(a, b)| a - b).collect();
        Ok(Matrix { rows: self.rows, cols: self.cols, data })
    }
}

impl Mul for &Matrix {
    type Output = Result<Matrix>;
    fn mul(self, rhs: &Matrix) -> Result<Matrix> {
        if self.cols != rhs.rows {
            return Err(MathError::InvalidArgument(format!(
                "mul: shape mismatch {}x{} * {}x{}",
                self.rows, self.cols, rhs.rows, rhs.cols
            )));
        }
        let mut data = vec![0.0; self.rows * rhs.cols];
        for i in 0..self.rows {
            for j in 0..rhs.cols {
                let mut sum = 0.0;
                for k in 0..self.cols {
                    sum += self[(i, k)] * rhs[(k, j)];
                }
                data[i * rhs.cols + j] = sum;
            }
        }
        Ok(Matrix { rows: self.rows, cols: rhs.cols, data })
    }
}

/// Scalar multiplication.
impl Matrix {
    pub fn scale(&self, s: f64) -> Matrix {
        let data: Vec<f64> = self.data.iter().map(|x| x * s).collect();
        Matrix { rows: self.rows, cols: self.cols, data }
    }
}

// --- Display ----------------------------------------------------------------

impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows {
            f.write_str("[")?;
            for j in 0..self.cols {
                if j > 0 {
                    f.write_str("  ")?;
                }
                write!(f, "{:>10.4}", self[(i, j)])?;
            }
            f.write_str("]")?;
            if i + 1 < self.rows {
                f.write_str("\n")?;
            }
        }
        Ok(())
    }
}

// --- Helper trait for swap_within (not in std for slices) --------------------

trait SwapWithin {
    fn swap_within(&mut self, a: std::ops::Range<usize>, b: std::ops::Range<usize>);
}

impl SwapWithin for Vec<f64> {
    fn swap_within(&mut self, a: std::ops::Range<usize>, b: std::ops::Range<usize>) {
        debug_assert_eq!(a.len(), b.len());
        for k in 0..a.len() {
            self.swap(a.start + k, b.start + k);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    fn vec_close(a: &[f64], b: &[f64]) -> bool {
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(x, y)| close(*x, *y))
    }

    #[test]
    fn identity_determinant() {
        let m = Matrix::identity(3);
        assert!(close(m.determinant().unwrap(), 1.0));
    }

    #[test]
    fn det_2x2() {
        let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        assert!(close(m.determinant().unwrap(), -2.0));
    }

    #[test]
    fn det_3x3() {
        // | 2 0 0 |
        // | 0 3 0 |  = 6
        // | 0 0 1 |
        let m = Matrix::from_rows(&[vec![2.0, 0.0, 0.0], vec![0.0, 3.0, 0.0], vec![0.0, 0.0, 1.0]]).unwrap();
        assert!(close(m.determinant().unwrap(), 6.0));
    }

    #[test]
    fn solve_3x3() {
        // x + y = 3, 2x - y = 0, x + 2y + z = 7
        // => x=1, y=2, z=2
        let a = Matrix::from_rows(&[
            vec![1.0, 1.0, 0.0],
            vec![2.0, -1.0, 0.0],
            vec![1.0, 2.0, 1.0],
        ]).unwrap();
        let b = vec![3.0, 0.0, 7.0];
        let x = a.solve(&b).unwrap();
        assert!(vec_close(&x, &[1.0, 2.0, 2.0]));
    }

    #[test]
    fn inverse_roundtrip() {
        let a = Matrix::from_rows(&[
            vec![4.0, 7.0],
            vec![2.0, 6.0],
        ]).unwrap();
        let inv = a.inverse().unwrap();
        let prod = (&a * &inv).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(close(prod[(i, j)], expected), "({},{}): got {} want {}", i, j, prod[(i, j)], expected);
            }
        }
    }

    #[test]
    fn matrix_mul() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        let b = Matrix::from_rows(&[vec![5.0, 6.0], vec![7.0, 8.0]]).unwrap();
        let c = (&a * &b).unwrap();
        assert_eq!(c[(0, 0)], 19.0);
        assert_eq!(c[(0, 1)], 22.0);
        assert_eq!(c[(1, 0)], 43.0);
        assert_eq!(c[(1, 1)], 50.0);
    }

    #[test]
    fn transpose() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]]).unwrap();
        let t = a.transpose();
        assert_eq!(t.rows, 3);
        assert_eq!(t.cols, 2);
        assert_eq!(t[(0, 0)], 1.0);
        assert_eq!(t[(0, 1)], 4.0);
        assert_eq!(t[(2, 1)], 6.0);
    }

    #[test]
    fn trace_works() {
        let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        assert!(close(m.trace().unwrap(), 5.0));
    }

    #[test]
    fn mul_vec_works() {
        let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        let v = vec![5.0, 6.0];
        let r = m.mul_vec(&v).unwrap();
        assert!(vec_close(&r, &[17.0, 39.0]));
    }

    #[test]
    fn add_sub() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        let b = Matrix::from_rows(&[vec![5.0, 6.0], vec![7.0, 8.0]]).unwrap();
        let sum = (&a + &b).unwrap();
        assert_eq!(sum[(0, 0)], 6.0);
        assert_eq!(sum[(1, 1)], 12.0);
        let diff = (&b - &a).unwrap();
        assert_eq!(diff[(0, 0)], 4.0);
        assert_eq!(diff[(1, 1)], 4.0);
    }
}
