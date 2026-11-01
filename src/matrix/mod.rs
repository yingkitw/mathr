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
        // The pivot threshold must be *relative to the scale of the matrix*.
        // An absolute threshold (the classic `1e-14`) wrongly reports a
        // tiny-scaled but nonsingular matrix as singular: `[[1e-15, 0], [0, 1]]`
        // has determinant `1e-15`, not `0`.
        //
        // The scale is the largest absolute entry of the matrix; a pivot is
        // treated as zero only when it is negligible *relative to that*.
        let scale = self
            .data
            .iter()
            .map(|v| v.abs())
            .fold(0.0f64, f64::max);
        let pivot_tol = if scale > 0.0 {
            scale * f64::EPSILON * n as f64
        } else {
            0.0
        };
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
            if best <= pivot_tol {
                // The column is (numerically) zero — the determinant is zero.
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
                    a.swap(col * (n + 1) + c, pivot * (n + 1) + c);
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

    /// Rank of the matrix, approximated by counting singular values above a
    /// tolerance using Gaussian elimination (a basic, O(n³) estimate).
    pub fn rank(&self, tol: f64) -> usize {
        let mut a = self.data.clone();
        let n = self.rows;
        let m = self.cols;
        let tol = if tol <= 0.0 { 1e-10 } else { tol };
        let mut rank = 0;
        let mut row = 0;
        for col in 0..m {
            if row >= n {
                break;
            }
            // partial pivot
            let mut pivot = row;
            let mut best = a[row * m + col].abs();
            for r in (row + 1)..n {
                let val = a[r * m + col].abs();
                if val > best {
                    best = val;
                    pivot = r;
                }
            }
            if best < tol {
                continue;
            }
            if pivot != row {
                for c in col..m {
                    a.swap(row * m + c, pivot * m + c);
                }
            }
            let piv = a[row * m + col];
            for r in (row + 1)..n {
                let factor = a[r * m + col] / piv;
                for c in col..m {
                    a[r * m + c] -= factor * a[row * m + c];
                }
            }
            row += 1;
            rank += 1;
        }
        rank
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

// --- Hilbert matrix and regularised solvers ---------------------------------


mod decomp;
mod eigen;

pub use decomp::{Cholesky, Lu, Qr, Svd};
pub use eigen::{EigenPair, PowerIterOptions};

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    /// Relative comparison — for values whose magnitude may be tiny.
    fn close_rel(a: f64, b: f64) -> bool {
        if a == b {
            return true;
        }
        let scale = a.abs().max(b.abs());
        (a - b).abs() <= 1e-9 * scale
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
    fn det_tiny_scaled_nonsingular_is_not_zero() {
        // Regression: the pivot threshold used to be absolute (1e-14), so a
        // tiny-scaled but nonsingular matrix was reported as singular.
        let m = Matrix::from_rows(&[vec![1e-15, 0.0], vec![0.0, 1.0]]).unwrap();
        let d = m.determinant().unwrap();
        assert!(d != 0.0, "expected a nonzero determinant, got {d}");
        assert!(close_rel(d, 1e-15), "got {d}");
        // A uniformly tiny matrix too.
        let m = Matrix::from_rows(&[vec![1e-20, 0.0], vec![0.0, 1e-20]]).unwrap();
        assert!(close_rel(m.determinant().unwrap(), 1e-40));
    }

    #[test]
    fn det_singular_matrices_are_still_zero() {
        // The scale-relative threshold must not make singular matrices
        // nonsingular.
        let m = Matrix::from_rows(&[vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
        assert_eq!(m.determinant().unwrap(), 0.0);
        let m = Matrix::from_rows(&[vec![0.0, 0.0], vec![0.0, 0.0]]).unwrap();
        assert_eq!(m.determinant().unwrap(), 0.0);
        // A singular *tiny* matrix is still singular.
        let m = Matrix::from_rows(&[vec![1e-15, 1e-15], vec![1e-15, 1e-15]]).unwrap();
        assert_eq!(m.determinant().unwrap(), 0.0);
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

    #[test]
    fn lu_round_trip() {
        let a = Matrix::from_rows(&[
            vec![2.0, 1.0, 1.0],
            vec![4.0, -6.0, 0.0],
            vec![-2.0, 7.0, 2.0],
        ]).unwrap();
        let fact = a.lu().unwrap();
        let reconstructed = fact.reconstruct();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(a[(i, j)], reconstructed[(i, j)]), "(i={}, j={}): {} vs {}", i, j, a[(i, j)], reconstructed[(i, j)]);
            }
        }
    }

    #[test]
    fn lu_determinant() {
        let a = Matrix::from_rows(&[
            vec![2.0, 1.0, 1.0],
            vec![4.0, -6.0, 0.0],
            vec![-2.0, 7.0, 2.0],
        ]).unwrap();
        let fact = a.lu().unwrap();
        // det = 2 * ((-6)*2 - 0*7) - 1*(4*2 - 0*-2) + 1*(4*7 - (-6)*-2)
        //    = 2*(-12) - 1*(8) + 1*(28 - 12)
        //    = -24 - 8 + 16 = -16
        assert!(close(fact.determinant(), -16.0));
    }

    #[test]
    fn lu_solve_matches_gauss() {
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![2.0, 5.0, 2.0],
            vec![3.0, 1.0, 4.0],
        ]).unwrap();
        let b = vec![14.0, 18.0, 20.0];
        let x_lu = a.lu().unwrap().solve(&b).unwrap();
        let x_gauss = a.solve(&b).unwrap();
        assert!(vec_close(&x_lu, &x_gauss));
    }

    #[test]
    fn lu_inverse_matches() {
        let a = Matrix::from_rows(&[
            vec![4.0, 7.0, 2.0],
            vec![2.0, 6.0, 1.0],
            vec![1.0, 2.0, 3.0],
        ]).unwrap();
        let inv_lu = a.lu().unwrap().inverse().unwrap();
        let inv_gauss = a.inverse().unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(inv_lu[(i, j)], inv_gauss[(i, j)]), "{} vs {}", inv_lu[(i, j)], inv_gauss[(i, j)]);
            }
        }
    }

    #[test]
    fn lu_with_pivoting() {
        // Zero in the pivot position forces a row swap.
        let a = Matrix::from_rows(&[
            vec![0.0, 2.0],
            vec![1.0, 3.0],
        ]).unwrap();
        let fact = a.lu().unwrap();
        assert!(close(fact.determinant(), -2.0));
        let x = fact.solve(&[4.0, 5.0]).unwrap();
        // 2*x2 = 4 => x2 = 2; x1 + 3*x2 = 5 => x1 = -1
        assert!(vec_close(&x, &[-1.0, 2.0]));
    }

    #[test]
    fn matrix_rank() {
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![2.0, 4.0, 6.0],
            vec![0.0, 1.0, 0.0],
        ]).unwrap();
        assert_eq!(a.rank(1e-10), 2);
    }

    #[test]
    fn cholesky_round_trip() {
        let a = Matrix::from_rows(&[
            vec![4.0, 12.0, -16.0],
            vec![12.0, 37.0, -43.0],
            vec![-16.0, -43.0, 98.0],
        ]).unwrap();
        let c = a.cholesky().unwrap();
        let reconstructed = c.reconstruct();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(a[(i, j)], reconstructed[(i, j)]),
                        "({},{}): got {} want {}", i, j, reconstructed[(i, j)], a[(i, j)]);
            }
        }
    }

    #[test]
    fn cholesky_solve_known() {
        // Solve A x = b with the Hilbert-like SPD matrix above.
        let a = Matrix::from_rows(&[
            vec![4.0, 12.0, -16.0],
            vec![12.0, 37.0, -43.0],
            vec![-16.0, -43.0, 98.0],
        ]).unwrap();
        let b = vec![1.0, 2.0, 3.0];
        let x_cho = a.cholesky().unwrap().solve(&b).unwrap();
        let x_ref = a.solve(&b).unwrap();
        assert!(vec_close(&x_cho, &x_ref));
    }

    #[test]
    fn cholesky_rejects_non_spd() {
        // Not positive-definite: diagonal element becomes negative.
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0],
            vec![2.0, 1.0],
        ]).unwrap();
        assert!(a.cholesky().is_err());
    }

    #[test]
    fn power_iteration_dominant() {
        // A simple 2x2 matrix with known eigenvalues.  A = [[2, 1], [1, 2]]
        // has eigenvalues 3 and 1; the dominant one is 3.
        let a = Matrix::from_rows(&[vec![2.0, 1.0], vec![1.0, 2.0]]).unwrap();
        let result = a.power_iteration(PowerIterOptions::default()).unwrap();
        assert!(close(result.value, 3.0), "got {}", result.value);
        // Eigenvector should be proportional to (1, 1) (or its negation).
        let prod = result.vector[0] * result.vector[1];
        assert!(prod > 0.0, "expected same-sign components, got {:?}", result.vector);
        // Normalised; magnitude ~ 1/√2.
        let mag = (result.vector[0] * result.vector[0] + result.vector[1] * result.vector[1]).sqrt();
        assert!(close(mag, 1.0));
    }

    #[test]
    fn power_iteration_3x3() {
        // A diagonal matrix has eigenvalues equal to its diagonals.
        let a = Matrix::from_rows(&[
            vec![2.0, 0.0, 0.0],
            vec![0.0, 5.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ]).unwrap();
        let result = a.power_iteration(PowerIterOptions::default()).unwrap();
        assert!(close(result.value, 5.0), "got {}", result.value);
        // Eigenvector for λ=5 should be (0, ±1, 0).
        assert!(result.vector[1].abs() > 0.99, "got {:?}", result.vector);
    }

    #[test]
    fn svd_reconstructs() {
        // A = U Σ Vᵀ where Σ has the singular values on its diagonal.
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0],
            vec![3.0, 4.0],
            vec![5.0, 6.0],
        ]).unwrap();
        let svd = a.svd().unwrap();
        // Σ is m × n = 3 × 2.
        let mut sigma = vec![0.0; a.rows * a.cols];
        for i in 0..a.cols.min(a.rows) {
            sigma[i * a.cols + i] = svd.singular_values[i];
        }
        let s_mat = Matrix::from_row_major(a.rows, a.cols, sigma).unwrap();
        let ut = svd.u.transpose();
        let product = (&(&ut * &a).unwrap() * &svd.v).unwrap(); // Uᵀ A V should equal Σ
        // Compare product to Σ.
        for r in 0..a.rows {
            for c in 0..a.cols {
                let expected = if r == c { svd.singular_values[r] } else { 0.0 };
                assert!((product[(r, c)] - expected).abs() < 1e-8,
                        "({}, {}): got {} want {}", r, c, product[(r, c)], expected);
            }
        }
        let _ = s_mat;
    }

    #[test]
    fn svd_diagonal_singular_values() {
        // For a diagonal matrix, singular values equal |diagonal entries|.
        let a = Matrix::from_rows(&[
            vec![3.0, 0.0, 0.0],
            vec![0.0, -2.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ]).unwrap();
        let svd = a.svd().unwrap();
        assert!((svd.singular_values[0] - 3.0).abs() < 1e-8);
        assert!((svd.singular_values[1] - 2.0).abs() < 1e-8);
        assert!((svd.singular_values[2] - 1.0).abs() < 1e-8);
    }

    #[test]
    fn condition_number_known_values() {
        // Identity: κ = 1
        let i = Matrix::identity(3);
        assert!((i.condition_number().unwrap() - 1.0).abs() < 1e-10);
        // diag(3, 1): κ = 3
        let d = Matrix::from_rows(&[vec![3.0, 0.0], vec![0.0, 1.0]]).unwrap();
        assert!((d.condition_number().unwrap() - 3.0).abs() < 1e-8);
        // Singular matrix: κ = inf
        let s = Matrix::from_rows(&[vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
        assert_eq!(s.condition_number().unwrap(), f64::INFINITY);
    }

    #[test]
    fn condition_number_hilbert_is_large() {
        // 6×6 Hilbert matrix: κ ≈ 1.5e7 — well-conditioned check boundary
        let h = Matrix::hilbert(6);
        let k = h.condition_number().unwrap();
        assert!(k > 1e6, "κ(H₆) = {}", k);
        assert!(k < 1e9, "κ(H₆) = {}", k);
    }

    #[test]
    fn nullspace_rank_deficient() {
        // [[1, 1], [1, 1]] → nullspace spanned by (−1, 1)/√2
        let s = Matrix::from_rows(&[vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
        let basis = s.nullspace(0.0).unwrap();
        assert_eq!(basis.len(), 1);
        // A·v = 0
        let av = s.mul_vec(&basis[0]).unwrap();
        assert!(av.iter().all(|&x| x.abs() < 1e-9));
    }

    #[test]
    fn nullspace_wide_matrix() {
        // 1×3 [1 2 3]: nullspace is a 2D plane
        let w = Matrix::from_rows(&[vec![1.0, 2.0, 3.0]]).unwrap();
        let basis = w.nullspace(0.0).unwrap();
        assert_eq!(basis.len(), 2);
        for v in &basis {
            assert!((v[0] + 2.0 * v[1] + 3.0 * v[2]).abs() < 1e-9);
            // orthonormal
            let n: f64 = v.iter().map(|x| x * x).sum();
            assert!((n - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn nullspace_full_rank_is_empty() {
        let a = Matrix::from_rows(&[vec![1.0, 0.0], vec![0.0, 1.0]]).unwrap();
        assert!(a.nullspace(0.0).unwrap().is_empty());
    }

    #[test]
    fn nullspace_vector_satisfies_homogeneous_system() {
        // 3×4 rank-2 system
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 0.0, 3.0],
            vec![2.0, 4.0, 0.0, 6.0],
            vec![1.0, 1.0, 1.0, 0.0],
        ]).unwrap();
        let basis = a.nullspace(0.0).unwrap();
        assert!(basis.len() >= 2);
        for v in &basis {
            let av = a.mul_vec(v).unwrap();
            assert!(av.iter().all(|&x| x.abs() < 1e-8));
        }
    }

    #[test]
    fn symmetric_eig_2x2() {
        // [[2, 1], [1, 2]] → eigenvalues 1, 3
        let a = Matrix::from_rows(&[vec![2.0, 1.0], vec![1.0, 2.0]]).unwrap();
        let (vals, vecs) = a.symmetric_eig().unwrap();
        assert!(close(vals[0], 1.0), "λ₀ = {}", vals[0]);
        assert!(close(vals[1], 3.0), "λ₁ = {}", vals[1]);
        // Verify A v = λ v for each eigenvector.
        for j in 0..2 {
            let v: Vec<f64> = (0..2).map(|i| vecs[(i, j)]).collect();
            let av = a.mul_vec(&v).unwrap();
            for i in 0..2 {
                assert!((av[i] - vals[j] * v[i]).abs() < 1e-9,
                        "A v_{} != λ_{} v_{} at i={}", j, j, j, i);
            }
        }
    }

    #[test]
    fn symmetric_eig_3x3_diagonal() {
        let a = Matrix::from_rows(&[
            vec![5.0, 0.0, 0.0],
            vec![0.0, 2.0, 0.0],
            vec![0.0, 0.0, 8.0],
        ]).unwrap();
        let (vals, _) = a.symmetric_eig().unwrap();
        assert!(close(vals[0], 2.0));
        assert!(close(vals[1], 5.0));
        assert!(close(vals[2], 8.0));
    }

    #[test]
    fn symmetric_eig_3x3_general() {
        // Symmetric 3×3 with known eigenvalues.
        let a = Matrix::from_rows(&[
            vec![4.0, 1.0, 2.0],
            vec![1.0, 3.0, 0.0],
            vec![2.0, 0.0, 5.0],
        ]).unwrap();
        let (vals, vecs) = a.symmetric_eig().unwrap();
        // Verify A v = λ v for each eigenvector.
        for j in 0..3 {
            let v: Vec<f64> = (0..3).map(|i| vecs[(i, j)]).collect();
            let av = a.mul_vec(&v).unwrap();
            for i in 0..3 {
                assert!((av[i] - vals[j] * v[i]).abs() < 1e-8,
                        "A v_{} != λ_{} v_{} at i={}: {} vs {}", j, j, j, i, av[i], vals[j] * v[i]);
            }
        }
        // Verify eigenvalues are sorted ascending.
        assert!(vals[0] <= vals[1] && vals[1] <= vals[2]);
        // Verify eigenvectors are orthonormal.
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| vecs[(k, i)] * vecs[(k, j)]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-8, "v_{}·v_{} = {}, expected {}", i, j, dot, expected);
            }
        }
    }

    #[test]
    fn symmetric_eig_4x4() {
        // Hilbert-like symmetric matrix.
        let a = Matrix::from_rows(&[
            vec![1.0, 0.5, 0.3333333333333333, 0.25],
            vec![0.5, 1.0, 0.5, 0.3333333333333333],
            vec![0.3333333333333333, 0.5, 1.0, 0.5],
            vec![0.25, 0.3333333333333333, 0.5, 1.0],
        ]).unwrap();
        let (vals, vecs) = a.symmetric_eig().unwrap();
        // All eigenvalues of this matrix should be positive.
        for &v in &vals {
            assert!(v > 0.0, "expected positive eigenvalue, got {}", v);
        }
        // Verify A v = λ v.
        for j in 0..4 {
            let v: Vec<f64> = (0..4).map(|i| vecs[(i, j)]).collect();
            let av = a.mul_vec(&v).unwrap();
            for i in 0..4 {
                assert!((av[i] - vals[j] * v[i]).abs() < 1e-7,
                        "A v_{} != λ_{} v_{} at i={}", j, j, j, i);
            }
        }
    }

    #[test]
    fn symmetric_eig_identity() {
        let a = Matrix::identity(5);
        let (vals, vecs) = a.symmetric_eig().unwrap();
        for &v in &vals {
            assert!(close(v, 1.0));
        }
        // Eigenvectors should form an orthonormal basis.
        for i in 0..5 {
            for j in 0..5 {
                let dot: f64 = (0..5).map(|k| vecs[(k, i)] * vecs[(k, j)]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-8);
            }
        }
    }

    #[test]
    fn hessenberg_3x3() {
        let a = Matrix::from_rows(&[
            vec![4.0, 1.0, 2.0],
            vec![1.0, 3.0, 0.0],
            vec![2.0, 0.0, 5.0],
        ]).unwrap();
        let (h, q) = a.hessenberg().unwrap();
        // H should be upper Hessenberg: zero below sub-diagonal.
        for i in 2..3 {
            for j in 0..i - 1 {
                assert!(h[(i, j)].abs() < 1e-10, "H[{}][{}] = {} should be 0", i, j, h[(i, j)]);
            }
        }
        // Verify A = Q H Q^T.
        let qt = q.transpose();
        let qhq = (&q * &h).unwrap();
        let qhqt = (&qhq * &qt).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((qhqt[(i, j)] - a[(i, j)]).abs() < 1e-10,
                    "QHQ^T[{}][{}] = {} != A[{}][{}] = {}", i, j, qhqt[(i, j)], i, j, a[(i, j)]);
            }
        }
        // Q should be orthogonal.
        let qtq = (&qt * &q).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((qtq[(i, j)] - expected).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn hessenberg_4x4_general() {
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0, 4.0],
            vec![5.0, 6.0, 7.0, 8.0],
            vec![9.0, 10.0, 11.0, 12.0],
            vec![13.0, 14.0, 15.0, 16.0],
        ]).unwrap();
        let (h, q) = a.hessenberg().unwrap();
        // Check upper Hessenberg structure.
        for i in 2..4 {
            for j in 0..i - 1 {
                assert!(h[(i, j)].abs() < 1e-10, "H[{}][{}] = {} should be 0", i, j, h[(i, j)]);
            }
        }
        // Verify A = Q H Q^T.
        let qt = q.transpose();
        let qhq = (&q * &h).unwrap();
        let qhqt = (&qhq * &qt).unwrap();
        for i in 0..4 {
            for j in 0..4 {
                assert!((qhqt[(i, j)] - a[(i, j)]).abs() < 1e-10,
                    "QHQ^T[{}][{}] = {} != A[{}][{}] = {}", i, j, qhqt[(i, j)], i, j, a[(i, j)]);
            }
        }
    }

    #[test]
    fn schur_symmetric_3x3() {
        // For a symmetric matrix, Schur form should be diagonal.
        let a = Matrix::from_rows(&[
            vec![4.0, 1.0, 2.0],
            vec![1.0, 3.0, 0.0],
            vec![2.0, 0.0, 5.0],
        ]).unwrap();
        let (t, q) = a.schur().unwrap();
        // Verify A = Q T Q^T.
        let qt = q.transpose();
        let qtq = (&qt * &q).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((qtq[(i, j)] - expected).abs() < 1e-8, "Q not orthogonal at ({},{})", i, j);
            }
        }
        let q_t_qt = (&(&q * &t).unwrap() * &qt).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((q_t_qt[(i, j)] - a[(i, j)]).abs() < 1e-8,
                    "QTQ^T[{}][{}] = {} != A[{}][{}] = {}", i, j, q_t_qt[(i, j)], i, j, a[(i, j)]);
            }
        }
        // For symmetric A, T should be (nearly) diagonal.
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    assert!(t[(i, j)].abs() < 1e-6, "T[{}][{}] = {} should be ~0", i, j, t[(i, j)]);
                }
            }
        }
        // Eigenvalues on diagonal should match symmetric_eig.
        let (eigvals, _) = a.symmetric_eig().unwrap();
        let mut schur_vals: Vec<f64> = (0..3).map(|i| t[(i, i)]).collect();
        schur_vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for i in 0..3 {
            assert!((schur_vals[i] - eigvals[i]).abs() < 1e-6,
                "Schur eigenvalue {} vs symmetric_eig {}", schur_vals[i], eigvals[i]);
        }
    }

    #[test]
    fn schur_general_3x3() {
        // Non-symmetric matrix with real eigenvalues.
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 10.0],
        ]).unwrap();
        let (t, q) = a.schur().unwrap();
        // Verify A = Q T Q^T.
        let qt = q.transpose();
        let q_t_qt = (&(&q * &t).unwrap() * &qt).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((q_t_qt[(i, j)] - a[(i, j)]).abs() < 1e-8,
                    "QTQ^T[{}][{}] = {} != A[{}][{}] = {}", i, j, q_t_qt[(i, j)], i, j, a[(i, j)]);
            }
        }
        // Q should be orthogonal.
        let qtq = (&qt * &q).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((qtq[(i, j)] - expected).abs() < 1e-8);
            }
        }
        // T should be upper triangular (quasi-triangular with 1x1 blocks for real eigenvalues).
        for i in 1..3 {
            assert!(t[(i, i - 1)].abs() < 1e-6, "T[{}][{}] = {} should be ~0", i, i - 1, t[(i, i - 1)]);
        }
    }

    #[test]
    fn schur_2x2_complex() {
        // Rotation matrix has complex eigenvalues e^{±iθ}.
        let theta: f64 = 0.5;
        let a = Matrix::from_rows(&[
            vec![theta.cos(), -theta.sin()],
            vec![theta.sin(), theta.cos()],
        ]).unwrap();
        let (t, q) = a.schur().unwrap();
        // Verify A = Q T Q^T.
        let qt = q.transpose();
        let q_t_qt = (&(&q * &t).unwrap() * &qt).unwrap();
        for i in 0..2 {
            for j in 0..2 {
                assert!((q_t_qt[(i, j)] - a[(i, j)]).abs() < 1e-8,
                    "QTQ^T[{}][{}] = {} != A[{}][{}] = {}", i, j, q_t_qt[(i, j)], i, j, a[(i, j)]);
            }
        }
        // T should be a 2×2 block with complex eigenvalues.
        let trace = t[(0, 0)] + t[(1, 1)];
        assert!((trace - 2.0 * theta.cos()).abs() < 1e-8, "trace = {} expected {}", trace, 2.0 * theta.cos());
        let det = t[(0, 0)] * t[(1, 1)] - t[(0, 1)] * t[(1, 0)];
        assert!((det - 1.0).abs() < 1e-8, "det = {} expected 1", det);
    }

    // --- Hilbert matrix tests ---

    #[test]
    fn hilbert_construction() {
        let h = Matrix::hilbert(3);
        assert!(close(h[(0, 0)], 1.0));
        assert!(close(h[(0, 1)], 0.5));
        assert!(close(h[(0, 2)], 1.0 / 3.0));
        assert!(close(h[(1, 0)], 0.5));
        assert!(close(h[(1, 1)], 1.0 / 3.0));
        assert!(close(h[(2, 2)], 1.0 / 5.0));
    }

    #[test]
    fn hilbert_symmetric() {
        let h = Matrix::hilbert(5);
        for i in 0..5 {
            for j in 0..5 {
                assert!(close(h[(i, j)], h[(j, i)]));
            }
        }
    }

    #[test]
    fn tikhonov_well_conditioned() {
        // For a well-conditioned system, Tikhonov with λ=0 should match plain solve
        let a = Matrix::from_rows(&[
            vec![4.0, 3.0, 2.0],
            vec![1.0, 5.0, 3.0],
            vec![2.0, 1.0, 6.0],
        ]).unwrap();
        let b = vec![20.0, 14.0, 15.0];
        let x_plain = a.solve(&b).unwrap();
        let x_tikh = a.solve_tikhonov(&b, 0.0).unwrap();
        assert!(vec_close(&x_plain, &x_tikh));
    }

    #[test]
    fn tikhonov_hilbert_stable() {
        // Hilbert matrix is notoriously ill-conditioned. Tikhonov regularisation
        // should produce a solution with smaller norm than the unregularised solve.
        let h = Matrix::hilbert(5);
        let x_true = vec![1.0, 1.0, 1.0, 1.0, 1.0];
        let b = h.mul_vec(&x_true).unwrap();

        // Plain solve on 5×5 Hilbert amplifies rounding — check regularisation
        // produces a bounded solution
        let x_reg = h.solve_tikhonov(&b, 1e-2).unwrap();
        let norm_reg: f64 = x_reg.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(norm_reg < 10.0, "regularised solution should be bounded: ||x|| = {}", norm_reg);
        // Each component should be within a reasonable range
        for &xi in &x_reg {
            assert!(xi.abs() < 5.0, "component {} too large", xi);
        }
    }

    #[test]
    fn tikhonov_shrinks_solution() {
        // Larger λ should produce smaller ||x||
        let h = Matrix::hilbert(4);
        let b = vec![1.0, 0.5, 0.333, 0.25];
        let x_small = h.solve_tikhonov(&b, 1e-10).unwrap();
        let x_large = h.solve_tikhonov(&b, 1.0).unwrap();
        let norm_small: f64 = x_small.iter().map(|x| x * x).sum::<f64>().sqrt();
        let norm_large: f64 = x_large.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(norm_large < norm_small,
            "larger λ should shrink ||x||: {} vs {}", norm_large, norm_small);
    }

    #[test]
    fn tikhonov_invalid_args() {
        let a = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        assert!(a.solve_tikhonov(&[1.0], 0.0).is_err()); // b length mismatch
        assert!(a.solve_tikhonov(&[1.0, 2.0], -1.0).is_err()); // negative lambda
    }

    #[test]
    fn tikhonov_rectangular() {
        // Overdetermined least-squares: A is 3×2, solve via Tikhonov with λ=0
        let a = Matrix::from_rows(&[
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
        ]).unwrap();
        // Data: y = 1 + 2x → b = [3, 5, 7]
        let b = vec![3.0, 5.0, 7.0];
        let x = a.solve_tikhonov(&b, 0.0).unwrap();
        assert!(close(x[0], 1.0), "intercept should be 1, got {}", x[0]);
        assert!(close(x[1], 2.0), "slope should be 2, got {}", x[1]);
    }

    // ===== QR Decomposition Tests =====

    #[test]
    fn qr_square_reconstruct() {
        let a = Matrix::from_rows(&[
            vec![12.0, -51.0, 4.0],
            vec![6.0, 167.0, -68.0],
            vec![-4.0, 24.0, -41.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let recon = qr.reconstruct();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(a[(i, j)], recon[(i, j)]),
                    "QR reconstruct mismatch at ({},{}): {} vs {}", i, j, a[(i,j)], recon[(i,j)]);
            }
        }
    }

    #[test]
    fn qr_q_orthogonal() {
        let a = Matrix::from_rows(&[
            vec![12.0, -51.0, 4.0],
            vec![6.0, 167.0, -68.0],
            vec![-4.0, 24.0, -41.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let q = qr.q();
        // QᵀQ should be identity
        let qt = q.transpose();
        let qtq = (&qt * &q).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(close(qtq[(i, j)], expected),
                    "QᵀQ not identity at ({},{}): {}", i, j, qtq[(i,j)]);
            }
        }
    }

    #[test]
    fn qr_r_upper_triangular() {
        let a = Matrix::from_rows(&[
            vec![12.0, -51.0, 4.0],
            vec![6.0, 167.0, -68.0],
            vec![-4.0, 24.0, -41.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let r = qr.r();
        for i in 0..3 {
            for j in 0..i {
                assert!(r[(i, j)].abs() < 1e-10,
                    "R not upper-triangular at ({},{}): {}", i, j, r[(i,j)]);
            }
        }
    }

    #[test]
    fn qr_rectangular_tall() {
        // 4×2 matrix
        let a = Matrix::from_rows(&[
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
            vec![1.0, 4.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let q = qr.q();
        let r = qr.r();
        assert_eq!(q.rows, 4);
        assert_eq!(q.cols, 4);
        assert_eq!(r.rows, 4);
        assert_eq!(r.cols, 2);
        // Reconstruct
        let recon = qr.reconstruct();
        for i in 0..4 {
            for j in 0..2 {
                assert!(close(a[(i, j)], recon[(i, j)]));
            }
        }
    }

    #[test]
    fn qr_rectangular_wide() {
        // 2×4 matrix (m < n)
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0, 4.0],
            vec![5.0, 6.0, 7.0, 8.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let q = qr.q();
        let r = qr.r();
        assert_eq!(q.rows, 2);
        assert_eq!(q.cols, 2);
        assert_eq!(r.rows, 2);
        assert_eq!(r.cols, 4);
        // Reconstruct
        let recon = qr.reconstruct();
        for i in 0..2 {
            for j in 0..4 {
                assert!(close(a[(i, j)], recon[(i, j)]));
            }
        }
    }

    #[test]
    fn qr_least_squares_linear_fit() {
        // Fit y = 1 + 2x to data (1,3), (2,5), (3,7)
        let a = Matrix::from_rows(&[
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
        ]).unwrap();
        let b = vec![3.0, 5.0, 7.0];
        let x = a.qr().unwrap().solve(&b).unwrap();
        assert!(close(x[0], 1.0), "intercept: {}", x[0]);
        assert!(close(x[1], 2.0), "slope: {}", x[1]);
    }

    #[test]
    fn qr_least_squares_noisy() {
        // Fit y = 2 + 0.5x with noise: (0, 2.1), (1, 2.4), (2, 3.1), (3, 3.4)
        // True line: y = 2 + 0.5x → b ≈ [2, 2.5, 3, 3.5]
        let a = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
        ]).unwrap();
        let b = vec![2.1, 2.4, 3.1, 3.4];
        let x = a.qr().unwrap().solve(&b).unwrap();
        // Should be close to [2, 0.5]
        assert!((x[0] - 2.0).abs() < 0.1, "intercept: {}", x[0]);
        assert!((x[1] - 0.5).abs() < 0.1, "slope: {}", x[1]);
    }

    #[test]
    fn qr_solve_matches_lu_for_square() {
        // For a square invertible system, QR and LU should give the same answer.
        let a = Matrix::from_rows(&[
            vec![4.0, 3.0],
            vec![6.0, 3.0],
        ]).unwrap();
        let b = vec![10.0, 12.0];
        let x_qr = a.qr().unwrap().solve(&b).unwrap();
        let x_lu = a.lu().unwrap().solve(&b).unwrap();
        for i in 0..2 {
            assert!(close(x_qr[i], x_lu[i]), "QR vs LU mismatch at {}: {} vs {}", i, x_qr[i], x_lu[i]);
        }
    }

    #[test]
    fn qr_identity_matrix() {
        let a = Matrix::from_rows(&[
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let recon = qr.reconstruct();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(a[(i, j)], recon[(i, j)]));
            }
        }
        // R should be diagonal with |1| on the diagonal (sign may vary)
        let r = qr.r();
        for i in 0..3 {
            assert!(close(r[(i, i)].abs(), 1.0), "R diagonal abs should be 1, got {}", r[(i,i)]);
        }
    }

    #[test]
    fn qr_already_upper_triangular() {
        let a = Matrix::from_rows(&[
            vec![2.0, 1.0, 3.0],
            vec![0.0, 5.0, 1.0],
            vec![0.0, 0.0, 4.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let recon = qr.reconstruct();
        for i in 0..3 {
            for j in 0..3 {
                assert!(close(a[(i, j)], recon[(i, j)]));
            }
        }
    }

    #[test]
    fn qr_rank_deficient_errors() {
        // Column 2 is 2× column 1 → rank deficient
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0],
            vec![1.0, 2.0],
            vec![1.0, 2.0],
        ]).unwrap();
        let b = vec![1.0, 2.0, 3.0];
        let result = a.qr().unwrap().solve(&b);
        assert!(result.is_err(), "rank-deficient should error");
    }

    #[test]
    fn qr_wrong_b_length_errors() {
        let a = Matrix::from_rows(&[
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
        ]).unwrap();
        let b = vec![1.0, 2.0]; // should be 3
        assert!(a.qr().unwrap().solve(&b).is_err());
    }

    #[test]
    fn qr_underdetermined_errors() {
        // 2×3 matrix (m < n) — underdetermined
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
        ]).unwrap();
        let b = vec![1.0, 2.0];
        assert!(a.qr().unwrap().solve(&b).is_err());
    }

    #[test]
    fn qr_q_orthogonal_rectangular() {
        let a = Matrix::from_rows(&[
            vec![1.0, 1.0],
            vec![1.0, 2.0],
            vec![1.0, 3.0],
            vec![1.0, 4.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let q = qr.q();
        let qt = q.transpose();
        let qtq = (&qt * &q).unwrap();
        // QᵀQ = I (4×4)
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(close(qtq[(i, j)], expected),
                    "QᵀQ not identity at ({},{}): {}", i, j, qtq[(i,j)]);
            }
        }
    }

    #[test]
    fn qr_classical_example() {
        // Classic textbook example (Trefethen & Bau)
        // A = [[1, 2], [2, 4], [3, 6]] is rank-deficient, so use a full-rank one:
        // A = [[1, 2], [2, 3], [3, 4]]
        let a = Matrix::from_rows(&[
            vec![1.0, 2.0],
            vec![2.0, 3.0],
            vec![3.0, 4.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let recon = qr.reconstruct();
        for i in 0..3 {
            for j in 0..2 {
                assert!(close(a[(i, j)], recon[(i, j)]));
            }
        }
        // R should be 3×2 upper-triangular (zeros below diagonal in 2×2 block)
        let r = qr.r();
        assert!(r[(1, 0)].abs() < 1e-10);
        assert!(r[(2, 0)].abs() < 1e-10);
        assert!(r[(2, 1)].abs() < 1e-10);
    }

    #[test]
    fn qr_2x2_simple() {
        // Simple 2×2: A = [[4, 1], [0, 3]]
        // QR: Q = [[1, 0], [0, 1]], R = A (already upper-triangular)
        let a = Matrix::from_rows(&[
            vec![4.0, 1.0],
            vec![0.0, 3.0],
        ]).unwrap();
        let qr = a.qr().unwrap();
        let recon = qr.reconstruct();
        assert!(close(recon[(0, 0)], 4.0));
        assert!(close(recon[(0, 1)], 1.0));
        assert!(close(recon[(1, 0)], 0.0));
        assert!(close(recon[(1, 1)], 3.0));
    }
}
