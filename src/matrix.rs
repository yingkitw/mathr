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
                    let tmp = a[row * m + c];
                    a[row * m + c] = a[pivot * m + c];
                    a[pivot * m + c] = tmp;
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

    /// LU decomposition with partial pivoting: returns a [`Lu`] object
    /// holding the combined L/U factors and a permutation vector such that
    /// `P · A = L · U`.
    ///
    /// L is unit-lower-triangular (1s on the diagonal), U is upper-triangular,
    /// and `P` is the row-permutation matrix induced by `piv`.
    pub fn lu(&self) -> Result<Lu> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument("lu: matrix must be square".into()));
        }
        let n = self.rows;
        let mut lu = self.data.clone();
        let mut piv: Vec<usize> = (0..n).collect();
        let mut sign = 1.0_f64;

        for k in 0..n {
            // Find pivot row
            let mut best = lu[k * n + k].abs();
            let mut pivot = k;
            for r in (k + 1)..n {
                let val = lu[r * n + k].abs();
                if val > best {
                    best = val;
                    pivot = r;
                }
            }
            if best < 1e-14 {
                return Err(MathError::InvalidArgument(
                    "lu: singular or near-singular matrix".into(),
                ));
            }
            if pivot != k {
                for c in 0..n {
                    lu.swap(k * n + c, pivot * n + c);
                }
                piv.swap(k, pivot);
                sign = -sign;
            }
            let piv_val = lu[k * n + k];
            for r in (k + 1)..n {
                lu[r * n + k] /= piv_val;
                let factor = lu[r * n + k];
                for c in (k + 1)..n {
                    lu[r * n + c] -= factor * lu[k * n + c];
                }
            }
        }
        Ok(Lu {
            n,
            lu,
            piv,
            sign,
        })
    }

    /// Cholesky decomposition: `A = L · Lᵀ` for symmetric positive-definite
    /// `A`. Returns the lower-triangular factor `L` wrapped in [`Cholesky`].
    /// Returns [`MathError::InvalidArgument`] if `A` is not square or if it
    /// is not positive-definite (a non-positive pivot or NaN is detected).
    pub fn cholesky(&self) -> Result<Cholesky> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument("cholesky: matrix must be square".into()));
        }
        let n = self.rows;
        let mut l = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let mut sum = self[(i, j)];
                for k in 0..j {
                    sum -= l[i * n + k] * l[j * n + k];
                }
                let value = if i == j {
                    if sum <= 0.0 || !sum.is_finite() {
                        return Err(MathError::InvalidArgument(
                            "cholesky: matrix is not positive-definite".into(),
                        ));
                    }
                    sum.sqrt()
                } else {
                    sum / l[j * n + j]
                };
                l[i * n + j] = value;
            }
        }
        Ok(Cholesky { n, l })
    }

    /// Hessenberg decomposition: `A = Q · H · Qᵀ` where `H` is upper
    /// Hessenberg (zeros below the first sub-diagonal) and `Q` is orthogonal.
    ///
    /// Uses Householder reflections to introduce zeros column by column.
    /// This is the standard first step before applying the QR algorithm for
    /// general (non-symmetric) eigenvalue problems.
    pub fn hessenberg(&self) -> Result<(Matrix, Matrix)> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument(
                "hessenberg: matrix must be square".into(),
            ));
        }
        let n = self.rows;
        if n <= 1 {
            return Ok((self.clone(), Matrix::identity(n)));
        }

        let mut h = self.data.clone();
        let mut q = vec![0.0; n * n];
        for i in 0..n { q[i * n + i] = 1.0; }

        for k in 0..n - 2 {
            // Norm of the sub-column below the diagonal (rows k+1..n).
            let mut sigma = 0.0;
            for i in (k + 1)..n {
                sigma += h[i * n + k] * h[i * n + k];
            }
            if sigma < 1e-30 {
                continue;
            }
            let alpha = if h[(k + 1) * n + k] >= 0.0 {
                -sigma.sqrt()
            } else {
                sigma.sqrt()
            };
            // Householder vector (indices k+1..n).
            let mut v = vec![0.0; n];
            v[k + 1] = h[(k + 1) * n + k] - alpha;
            for i in (k + 2)..n {
                v[i] = h[i * n + k];
            }
            let v_norm_sq: f64 = (k + 1..n).map(|i| v[i] * v[i]).sum();
            if v_norm_sq < 1e-30 {
                continue;
            }
            let beta = 2.0 / v_norm_sq;

            // Apply H from the left: H_left = I - beta v v^T
            // H ← H_left · H  (affects rows k+1..n, all columns)
            for j in 0..n {
                let mut s = 0.0;
                for i in (k + 1)..n {
                    s += v[i] * h[i * n + j];
                }
                let factor = beta * s;
                for i in (k + 1)..n {
                    h[i * n + j] -= factor * v[i];
                }
            }
            // Apply H from the right: H ← H · H_left  (affects all rows, cols k+1..n)
            for i in 0..n {
                let mut s = 0.0;
                for j in (k + 1)..n {
                    s += h[i * n + j] * v[j];
                }
                let factor = beta * s;
                for j in (k + 1)..n {
                    h[i * n + j] -= factor * v[j];
                }
            }
            // Accumulate Q ← Q · H_left (right-multiply, cols k+1..n)
            for i in 0..n {
                let mut s = 0.0;
                for j in (k + 1)..n {
                    s += q[i * n + j] * v[j];
                }
                let factor = beta * s;
                for j in (k + 1)..n {
                    q[i * n + j] -= factor * v[j];
                }
            }
        }

        // Clean up negligible elements below the sub-diagonal.
        for i in 2..n {
            for j in 0..i - 1 {
                h[i * n + j] = 0.0;
            }
        }

        Ok((
            Matrix::from_row_major(n, n, h)?,
            Matrix::from_row_major(n, n, q)?,
        ))
    }

    /// Real Schur decomposition: `A = Q · T · Qᵀ` where `T` is quasi-upper
    /// triangular (1×1 blocks for real eigenvalues, 2×2 blocks for complex
    /// conjugate eigenvalue pairs) and `Q` is orthogonal.
    ///
    /// Reduces the matrix to Hessenberg form, then applies the shifted QR
    /// algorithm with Wilkinson shifts and deflation.  Returns the Schur
    /// form `T` and the orthogonal matrix `Q`.
    pub fn schur(&self) -> Result<(Matrix, Matrix)> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument(
                "schur: matrix must be square".into(),
            ));
        }
        let n = self.rows;
        if n == 0 {
            return Ok((Matrix::identity(0), Matrix::identity(0)));
        }
        if n == 1 {
            return Ok((self.clone(), Matrix::identity(1)));
        }

        // Step 1: Reduce to Hessenberg form.
        let (mut t, mut q) = self.hessenberg()?;
        let mut t_data = std::mem::take(&mut t.data);
        let mut q_data = std::mem::take(&mut q.data);

        let eps = 1e-14;
        let max_iter = 300 * n;
        let mut iter = 0;
        let mut m = n; // active sub-matrix [0, m)

        while m > 1 {
            // Find the smallest l such that the sub-diagonal t[l][l-1] is negligible.
            let mut l = m - 1;
            loop {
                if l == 0 { break; }
                let off = t_data[l * n + l - 1].abs();
                let diag = t_data[(l - 1) * n + (l - 1)].abs() + t_data[l * n + l].abs();
                if off <= eps * diag {
                    t_data[l * n + l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }

            if l == m - 1 {
                // 1×1 block converged.
                m -= 1;
                continue;
            }

            // Check for 2×2 block convergence at the bottom.
            if l == m - 2 {
                // Check if the 2×2 block has complex eigenvalues.
                let a = t_data[(m - 2) * n + (m - 2)];
                let b = t_data[(m - 2) * n + (m - 1)];
                let c = t_data[(m - 1) * n + (m - 2)];
                let d = t_data[(m - 1) * n + (m - 1)];
                let disc = (a - d) * (a - d) + 4.0 * b * c;
                if disc < 0.0 {
                    // Complex conjugate pair — 2×2 block is already in Schur form.
                    m -= 2;
                    continue;
                }
                // Real eigenvalues — compute the Wilkinson shift and do one more sweep.
                // (Fall through to the QR step below.)
            }

            iter += 1;
            if iter > max_iter {
                return Err(MathError::NotConvergent(
                    "schur: QR iteration did not converge".into(),
                ));
            }

            // Wilkinson shift from the trailing 2×2 block.
            let a = t_data[(m - 2) * n + (m - 2)];
            let b = t_data[(m - 2) * n + (m - 1)];
            let c = t_data[(m - 1) * n + (m - 2)];
            let d = t_data[(m - 1) * n + (m - 1)];
            let disc = (a - d) * (a - d) + 4.0 * b * c;
            let mu = if disc < 0.0 {
                // Complex eigenvalues — use the real part as shift.
                (a + d) / 2.0
            } else {
                let sqrt_disc = disc.sqrt();
                let lambda1 = (a + d + sqrt_disc) / 2.0;
                let lambda2 = (a + d - sqrt_disc) / 2.0;
                // Choose the eigenvalue closer to d (the bottom-right entry).
                if (lambda1 - d).abs() < (lambda2 - d).abs() {
                    lambda1
                } else {
                    lambda2
                }
            };

            // Shift: T ← T - μI
            for i in l..m {
                t_data[i * n + i] -= mu;
            }

            // QR factorization of the active Hessenberg block [l, m).
            // Since the matrix is upper Hessenberg, only one Givens rotation per
            // column is needed (to zero the single sub-diagonal element).
            let mut rots: Vec<(usize, f64, f64)> = Vec::with_capacity(m - l - 1);
            for k in l..m - 1 {
                let x = t_data[k * n + k];
                let z = t_data[(k + 1) * n + k];
                let (c, s) = givens(x, z);
                rots.push((k, c, s));
                // Apply G to rows k, k+1 (all columns, since Hessenberg has fill-in above).
                for j in k..n {
                    let ak = t_data[k * n + j];
                    let ak1 = t_data[(k + 1) * n + j];
                    t_data[k * n + j] = c * ak + s * ak1;
                    t_data[(k + 1) * n + j] = -s * ak + c * ak1;
                }
            }

            // RQ: apply G_k^T from the right (columns k, k+1), forward order.
            for &(k, c, s) in &rots {
                for i in 0..(k + 2).min(m) {
                    let aik = t_data[i * n + k];
                    let aik1 = t_data[i * n + k + 1];
                    t_data[i * n + k] = c * aik + s * aik1;
                    t_data[i * n + k + 1] = -s * aik + c * aik1;
                }
            }

            // Undo shift: T ← T + μI
            for i in l..m {
                t_data[i * n + i] += mu;
            }

            // Accumulate eigenvectors: Q ← Q · G_k^T (forward order, all rows).
            for &(k, c, s) in &rots {
                for i in 0..n {
                    let qik = q_data[i * n + k];
                    let qik1 = q_data[i * n + k + 1];
                    q_data[i * n + k] = c * qik + s * qik1;
                    q_data[i * n + k + 1] = -s * qik + c * qik1;
                }
            }
        }

        Ok((
            Matrix::from_row_major(n, n, t_data)?,
            Matrix::from_row_major(n, n, q_data)?,
        ))
    }

    /// Symmetric eigenvalue decomposition via the QR algorithm.
    ///
    /// Reduces the symmetric matrix to tridiagonal form using Householder
    /// reflections, then applies the shifted QR iteration (Wilkinson shift)
    /// with Givens rotations until all off-diagonal elements are negligible.
    /// Returns the eigenvalues in ascending order and the corresponding
    /// eigenvectors as columns of the returned matrix.
    ///
    /// The input is **not** checked for symmetry; if the matrix is not
    /// (approximately) symmetric the result is undefined.
    pub fn symmetric_eig(&self) -> Result<(Vec<f64>, Matrix)> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument(
                "symmetric_eig: matrix must be square".into(),
            ));
        }
        let n = self.rows;
        if n == 0 {
            return Ok((Vec::new(), Matrix::identity(0)));
        }
        if n == 1 {
            return Ok((vec![self[(0, 0)]], Matrix::identity(1)));
        }

        // --- Step 1: Householder tridiagonalisation ---
        let mut a = self.data.clone();
        let mut q = vec![0.0; n * n];
        for i in 0..n { q[i * n + i] = 1.0; }

        for k in 0..n - 2 {
            // Compute the norm of the sub-column below the diagonal.
            let mut sigma = 0.0;
            for i in (k + 1)..n {
                sigma += a[i * n + k] * a[i * n + k];
            }
            if sigma < 1e-30 {
                continue;
            }
            let alpha = if a[(k + 1) * n + k] >= 0.0 {
                -sigma.sqrt()
            } else {
                sigma.sqrt()
            };
            // Householder vector v (only indices k+1..n are nonzero).
            let mut v = vec![0.0; n];
            v[k + 1] = a[(k + 1) * n + k] - alpha;
            for i in (k + 2)..n {
                v[i] = a[i * n + k];
            }
            let v_norm_sq: f64 = (k + 1..n).map(|i| v[i] * v[i]).sum();
            if v_norm_sq < 1e-30 {
                continue;
            }
            let beta = 2.0 / v_norm_sq;

            // Apply H = I - beta * v * v^T  from both sides: A ← H A H.
            // Since H is symmetric, this is A ← H A H.
            // Compute p = beta * A * v  (only rows 0..n, cols k+1..n matter).
            let mut p = vec![0.0; n];
            for i in 0..n {
                let mut s = 0.0;
                for j in (k + 1)..n {
                    s += a[i * n + j] * v[j];
                }
                p[i] = beta * s;
            }
            // K = beta * v^T * p / 2  (so that q = p - K*v makes the update symmetric)
            let k_dot: f64 = beta * (k + 1..n).map(|i| v[i] * p[i]).sum::<f64>() / 2.0;
            // q_vec = p - K * v
            let mut q_vec = vec![0.0; n];
            for i in 0..n {
                q_vec[i] = p[i] - k_dot * v[i];
            }
            // A ← A - v * q^T - q * v^T
            for i in 0..n {
                for j in 0..n {
                    a[i * n + j] -= v[i] * q_vec[j] + q_vec[i] * v[j];
                }
            }
            // Accumulate Q ← Q * H (right-multiply).
            for i in 0..n {
                let mut s = 0.0;
                for j in (k + 1)..n {
                    s += q[i * n + j] * v[j];
                }
                let factor = beta * s;
                for j in (k + 1)..n {
                    q[i * n + j] -= factor * v[j];
                }
            }
        }

        // --- Step 2: QR iteration on the tridiagonal matrix ---
        // Extract diagonal d[0..n] and off-diagonal e[0..n-1].
        let mut d: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
        let mut e: Vec<f64> = (0..n - 1).map(|i| a[(i + 1) * n + i]).collect();
        let mut vecs = q;

        let eps = 1e-14;
        let max_iter = 200 * n;
        let mut m = n;
        let mut iter = 0;

        while m > 1 {
            // Find the smallest l such that e[l-1] is negligible.
            let mut l = m - 1;
            loop {
                if l == 0 { break; }
                if e[l - 1].abs() <= eps * (d[l - 1].abs() + d[l].abs()) {
                    e[l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }

            if l == m - 1 {
                m -= 1;
                continue;
            }

            iter += 1;
            if iter > max_iter {
                return Err(MathError::NotConvergent(
                    "symmetric_eig: QR iteration did not converge".into(),
                ));
            }

            // Wilkinson shift from the trailing 2×2 block.
            let dd = d[m - 2] - d[m - 1];
            let mu = if dd.abs() < 1e-300 {
                d[m - 1] - e[m - 2].abs()
            } else {
                let t = e[m - 2] / dd;
                let sgn = if dd >= 0.0 { 1.0 } else { -1.0 };
                d[m - 1] - sgn * e[m - 2].abs() * (t / (1.0 + (t * t).sqrt()))
            };

            // Shift.
            for i in l..m {
                d[i] -= mu;
            }

            // Build the active block (already shifted since d was shifted).
            let bs = m - l;
            let mut block = vec![0.0; bs * bs];
            for i in 0..bs {
                block[i * bs + i] = d[l + i];
                if i + 1 < bs {
                    block[i * bs + i + 1] = e[l + i];
                    block[(i + 1) * bs + i] = e[l + i];
                }
            }

            // QR factorization via Givens (one per sub-diagonal element, since tridiagonal).
            let mut rots2: Vec<(usize, f64, f64)> = Vec::with_capacity(bs - 1);
            for k in 0..bs - 1 {
                let (c, s) = givens(block[k * bs + k], block[(k + 1) * bs + k]);
                rots2.push((k, c, s));
                for j in k..bs {
                    let ak = block[k * bs + j];
                    let ak1 = block[(k + 1) * bs + j];
                    block[k * bs + j] = c * ak + s * ak1;
                    block[(k + 1) * bs + j] = -s * ak + c * ak1;
                }
            }

            // RQ: apply G_k^T from the right (columns k, k+1), forward order.
            for &(k, c, s) in &rots2 {
                for i in 0..bs {
                    let aik = block[i * bs + k];
                    let aik1 = block[i * bs + k + 1];
                    block[i * bs + k] = c * aik + s * aik1;
                    block[i * bs + k + 1] = -s * aik + c * aik1;
                }
            }

            // Undo shift on block.
            for i in 0..bs {
                block[i * bs + i] += mu;
            }

            // Extract tridiagonal from block.
            for i in 0..bs {
                d[l + i] = block[i * bs + i];
                if i + 1 < bs {
                    e[l + i] = block[i * bs + i + 1];
                }
            }

            // Accumulate eigenvectors: V ← V * Q (forward order).
            for &(k, c, s) in &rots2 {
                let col_k = l + k;
                let col_k1 = l + k + 1;
                for i in 0..n {
                    let vik = vecs[i * n + col_k];
                    let vik1 = vecs[i * n + col_k1];
                    vecs[i * n + col_k] = c * vik + s * vik1;
                    vecs[i * n + col_k1] = -s * vik + c * vik1;
                }
            }
        }

        // --- Step 3: Sort eigenvalues ascending, reorder eigenvectors ---
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap_or(std::cmp::Ordering::Equal));
        let eigenvalues: Vec<f64> = idx.iter().map(|&i| d[i]).collect();
        let mut eigvecs = vec![0.0; n * n];
        for (col, &i) in idx.iter().enumerate() {
            for r in 0..n {
                eigvecs[r * n + col] = vecs[r * n + i];
            }
        }

        Ok((eigenvalues, Matrix::from_row_major(n, n, eigvecs)?))
    }

    /// Power iteration: returns the dominant eigenvalue and a corresponding
    /// eigenvector using the iterative scheme `v ← A · v`, `v ← v / ‖v‖`.
    /// The starting vector is `(1, 1, …, 1)` (in the standard basis).
    ///
    /// Converges when there is a single strictly dominant eigenvalue; for
    /// matrices with tied leading eigenvalues the iterates may oscillate
    /// and the answer is best-effort within `tol` / `max_iter`.
    pub fn power_iteration(&self, options: PowerIterOptions) -> Result<EigenPair> {
        if !self.is_square() {
            return Err(MathError::InvalidArgument(
                "power_iteration: matrix must be square".into(),
        ));
        }
        let n = self.rows;
        let mut v = vec![1.0_f64; n];
        let mut lambda = 0.0_f64;
        for _ in 0..options.max_iter {
            // Av = self * v
            let av = self.mul_vec(&v)?;
            let norm = av.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                return Err(MathError::NotConvergent(
                    "power_iteration: zero norm (singular matrix?)".into(),
                ));
            }
            v = av.iter().map(|x| x / norm).collect();
            // Rayleigh quotient: λ ≈ vᵀ A v / vᵀ v = vᵀ (A v) since ‖v‖ = 1.
            let av_again = self.mul_vec(&v)?;
            let new_lambda: f64 = v.iter().zip(av_again.iter()).map(|(a, b)| a * b).sum();
            if (new_lambda - lambda).abs() < options.tol {
                lambda = new_lambda;
                break;
            }
            lambda = new_lambda;
        }
        Ok(EigenPair { value: lambda, vector: v })
    }
}

/// Cholesky decomposition `A = L · Lᵀ` for a symmetric positive-definite matrix.
/// Holds the lower-triangular `L` factor.
#[derive(Debug, Clone)]
pub struct Cholesky {
    n: usize,
    /// Lower-triangular Cholesky factor in column-major-friendly flat storage.
    l: Vec<f64>,
}

impl Cholesky {
    /// Solve `A · x = b` using the precomputed Cholesky factor `L`
    /// (`A · x = b  ⇔  L · y = b,  Lᵀ · x = y`).
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        let n = self.n;
        if b.len() != n {
            return Err(MathError::InvalidArgument(format!(
                "cholesky.solve: b has length {}, expected {}",
                b.len(), n
            )));
        }
        // Forward sub: L · y = b. L is lower-triangular with L[i][i] = l[i*n+i].
        let mut y = vec![0.0; n];
        for i in 0..n {
            let mut sum = b[i];
            for k in 0..i {
                sum -= self.l[i * n + k] * y[k];
            }
            y[i] = sum / self.l[i * n + i];
        }
        // Back sub: Lᵀ · x = y.
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut sum = y[i];
            for k in (i + 1)..n {
                sum -= self.l[k * n + i] * x[k];
            }
            x[i] = sum / self.l[i * n + i];
        }
        Ok(x)
    }

    /// Recover the original matrix `L · Lᵀ`.
    pub fn reconstruct(&self) -> Matrix {
        let n = self.n;
        let mut a = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..=r {
                let mut sum = 0.0;
                for k in 0..=r.min(c) {
                    sum += self.l[r * n + k] * self.l[c * n + k];
                }
                a[r * n + c] = sum;
                a[c * n + r] = sum; // symmetric
            }
        }
        Matrix::from_row_major(n, n, a).unwrap()
    }

    /// Return the lower-triangular factor as its own matrix.
    pub fn l_factor(&self) -> Matrix {
        let n = self.n;
        let mut l_data = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..=r {
                l_data[r * n + c] = self.l[r * n + c];
            }
        }
        Matrix::from_row_major(n, n, l_data).unwrap()
    }
}

/// Options for [`Matrix::power_iteration`].
#[derive(Debug, Clone, Copy)]
pub struct PowerIterOptions {
    /// Max iterations.
    pub max_iter: usize,
    /// Convergence tolerance on `|λᵏ⁺¹ − λᵏ|`.
    pub tol: f64,
}

impl Default for PowerIterOptions {
    fn default() -> Self {
        Self { max_iter: 1000, tol: 1e-10 }
    }
}

/// A converged power-iteration result: `(λ, v)` where `v` is the eigenvector
/// corresponding to the dominant eigenvalue `λ`.
#[derive(Debug, Clone)]
pub struct EigenPair {
    pub value: f64,
    pub vector: Vec<f64>,
}

impl std::ops::Index<usize> for EigenPair {
    type Output = f64;
    fn index(&self, i: usize) -> &f64 { &self.vector[i] }
}

/// Singular Value Decomposition: `A = U · Σ · Vᵀ` for any `m × n` matrix
/// (possibly rectangular). The columns of `u` are the left singular vectors,
/// `singular_values` are the diagonal of `Σ` (in descending order), and the
/// columns of `v` are the right singular vectors.
#[derive(Debug, Clone)]
pub struct Svd {
    pub u: Matrix,
    pub singular_values: Vec<f64>,
    pub v: Matrix,
}

impl Matrix {
    /// Compute the singular value decomposition `A = U · Σ · Vᵀ` using the
    /// one-sided Jacobi rotation algorithm. Handles rectangular `m × n`
    /// matrices where `m ≥ n` (the more general case reduces to this by
    /// transposition, which doesn't change singular values).
    ///
    /// The largest off-diagonal element of `Aᵀ · A` is zeroed out by a
    /// Givens rotation, and the accumulated rotations form `V`.
    /// `U = A · V · Σ⁻¹` (with σ⁻¹ = 0 for any zero singular value).
    pub fn svd(&self) -> Result<Svd> {
        let m = self.rows;
        let n = self.cols;
        if m == 0 || n == 0 {
            return Err(MathError::InvalidArgument("svd: empty matrix".into()));
        }
        // Use Aᵀ · A (an n×n symmetric matrix) to drive Jacobi rotations.
        // For efficiency when m < n, work on A · Aᵀ (m×m) and swap roles.
        let (work_mat, work_m, work_n, left_is_u) = if m >= n {
            (self.clone(), m, n, true)
        } else {
            // Use the transpose: singular values are the same.
            let at = self.transpose();
            (at, n, m, false)
        };
        // B = work_matᵀ · work_mat  (work_n × work_n symmetric)
        let mut b = vec![0.0; work_n * work_n];
        for r in 0..work_n {
            for c in r..work_n {
                let mut s = 0.0;
                for k in 0..work_m {
                    s += work_mat[(k, r)] * work_mat[(k, c)];
                }
                b[r * work_n + c] = s;
                b[c * work_n + r] = s;
            }
        }
        // V accumulates the rotations; start as identity.
        let mut v = vec![0.0; work_n * work_n];
        for i in 0..work_n { v[i * work_n + i] = 1.0; }
        let max_sweeps = 100usize;
        let tol = 1e-14;
        let mut singular_values = vec![0.0; work_n];
        for _ in 0..max_sweeps {
            // Find largest off-diagonal element.
            let mut p = 0;
            let mut q = 1;
            let mut max_off = 0.0;
            for r in 0..work_n {
                for c in (r + 1)..work_n {
                    let v = b[r * work_n + c].abs();
                    if v > max_off {
                        max_off = v;
                        p = r;
                        q = c;
                    }
                }
            }
            if max_off < tol {
                break;
            }
            // Compute rotation: rotate in (p, q) plane to zero b[p][q].
            let app = b[p * work_n + p];
            let aqq = b[q * work_n + q];
            let apq = b[p * work_n + q];
            let theta = if (app - aqq).abs() < 1e-30 {
                std::f64::consts::FRAC_PI_4
            } else {
                0.5 * ((2.0 * apq) / (aqq - app)).atan()
            };
            let (s, c) = theta.sin_cos();
            // Apply rotation to B: B ← R · B · Rᵀ.
            for i in 0..work_n {
                if i != p && i != q {
                    let bip = b[i * work_n + p];
                    let biq = b[i * work_n + q];
                    b[i * work_n + p] = c * bip - s * biq;
                    b[p * work_n + i] = b[i * work_n + p];
                    b[i * work_n + q] = s * bip + c * biq;
                    b[q * work_n + i] = b[i * work_n + q];
                }
            }
            let bpp = c * c * app - 2.0 * s * c * apq + s * s * aqq;
            let bqq = s * s * app + 2.0 * s * c * apq + c * c * aqq;
            b[p * work_n + p] = bpp;
            b[q * work_n + q] = bqq;
            b[p * work_n + q] = 0.0;
            b[q * work_n + p] = 0.0;
            // Accumulate V ← V · Rᵀ (= V · Rᵀ since R is in cols p, q with Rᵀ).
            // Equivalent update: column p ← c · v_p − s · v_q, column q ← s · v_p + c · v_q.
            for i in 0..work_n {
                let vip = v[i * work_n + p];
                let viq = v[i * work_n + q];
                v[i * work_n + p] = c * vip - s * viq;
                v[i * work_n + q] = s * vip + c * viq;
            }
        }
        // Singular values are sqrt(eigenvalues), sorted descending.
        let mut eig_pairs: Vec<(f64, usize)> = (0..work_n)
            .map(|i| (b[i * work_n + i], i))
            .collect();
        eig_pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        for (i, (val, _)) in eig_pairs.iter().enumerate() {
            singular_values[i] = if *val < 0.0 { 0.0 } else { val.sqrt() };
        }
        // Reorder V columns by the same ordering.
        let mut v_re = vec![0.0; work_n * work_n];
        for (new_col, (_, old_col)) in eig_pairs.iter().enumerate() {
            for r in 0..work_n {
                v_re[r * work_n + new_col] = v[r * work_n + *old_col];
            }
        }
        let v_mat = Matrix::from_row_major(work_n, work_n, v_re).unwrap();
        // U = A · V · diag(1/σ_i) (clipped by σ_min threshold).
        let _ = left_is_u;
        // Recompute U from the original A (we may have transposed). Use self
        // directly so result dimensions are m × m.
        let u_rows = self.rows;
        let u_cols = work_n;
        let mut u_data = vec![0.0; u_rows * u_cols];
        // Anything smaller than this is treated as zero.
        let sigma_thresh = 1e-12;
        for j in 0..u_cols {
            let sigma = singular_values[j];
            for r in 0..u_rows {
                let mut sum = 0.0;
                for k in 0..self.cols {
                    sum += self[(r, k)] * v_mat[(k, j)];
                }
                u_data[r * u_cols + j] = if sigma > sigma_thresh { sum / sigma } else { 0.0 };
            }
        }
        // Pad u_data to be square (u_rows × u_rows). The extra columns for
        // m > n are orthonormal complements (not fully computed here; we
        // pad with zeros — sufficient for A·V·Σ⁻¹ = U[:, :k]).
        let mut u_full = vec![0.0; u_rows * u_rows];
        for r in 0..u_rows {
            for c in 0..u_cols {
                u_full[r * u_rows + c] = u_data[r * u_cols + c];
            }
            // The remaining columns are zero — these are not used by our
            // downstream consumers since Σ's tail entries are also zero.
        }
        Ok(Svd {
            u: Matrix::from_row_major(u_rows, u_rows, u_full).unwrap(),
            singular_values,
            v: v_mat,
        })
    }
}

/// A row-major LU decomposition of a square matrix.
///
/// Holds `lu` such that the strict lower triangle (with the 1-implicit
/// diagonal factored out) is L and the upper triangle (including diagonal)
/// is U. Together with `piv`, `P · A = L · U` where `P` is the implicit
/// row permutation.
#[derive(Debug, Clone)]
pub struct Lu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
    sign: f64,
}

impl Lu {
    /// Determinant of the original matrix from the LU factors.
    pub fn determinant(&self) -> f64 {
        let mut d = self.sign;
        for i in 0..self.n {
            d *= self.lu[i * self.n + i];
        }
        d
    }

    /// Recover the original matrix `A = P⁻¹ · L · U`.
    pub fn reconstruct(&self) -> Matrix {
        let n = self.n;
        // Compute L*U (which equals P*A by construction).
        let mut pa = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                let mut sum = 0.0;
                let hi = r.min(c);
                for k in 0..=hi {
                    let l = if k == r { 1.0 } else { self.lu[r * n + k] };
                    let u = self.lu[k * n + c];
                    sum += l * u;
                }
                pa[r * n + c] = sum;
            }
        }
        // Invert piv: inv_piv[i] is the row whose pivot ended up at position i.
        let mut inv_piv = vec![0usize; n];
        for r in 0..n {
            inv_piv[self.piv[r]] = r;
        }
        // A[inv_piv_target][c] = (P*A)[source][c]
        let mut data = vec![0.0; n * n];
        for r in 0..n {
            for c in 0..n {
                data[inv_piv[r] * n + c] = pa[r * n + c];
            }
        }
        Matrix::from_row_major(n, n, data).unwrap()
    }

    /// Solve `A · x = b` using the LU factors. `b` may be a vector or a
    /// matrix whose columns are multiple right-hand sides.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        if b.len() != self.n {
            return Err(MathError::InvalidArgument(format!(
                "lu.solve: b has length {}, expected {}",
                b.len(),
                self.n
            )));
        }
        // Apply permutation: P·b
        let mut pb: Vec<f64> = (0..self.n).map(|i| b[self.piv[i]]).collect();
        // Forward substitution: L·y = P·b  (L is unit-lower, skip diagonal)
        for r in 1..self.n {
            let mut sum = pb[r];
            for c in 0..r {
                sum -= self.lu[r * self.n + c] * pb[c];
            }
            pb[r] = sum;
        }
        // Back substitution: U·x = y
        let mut x = vec![0.0; self.n];
        for r in (0..self.n).rev() {
            let mut sum = pb[r];
            for c in (r + 1)..self.n {
                sum -= self.lu[r * self.n + c] * x[c];
            }
            x[r] = sum / self.lu[r * self.n + r];
        }
        Ok(x)
    }

    /// Compute the inverse of the original matrix by solving `A · x = e_i`
    /// for each column of the identity.
    pub fn inverse(&self) -> Result<Matrix> {
        let n = self.n;
        let mut inv = vec![0.0; n * n];
        let mut e = vec![0.0; n];
        for i in 0..n {
            e[i] = 1.0;
            let col = self.solve(&e)?;
            for r in 0..n {
                inv[r * n + i] = col[r];
            }
            e[i] = 0.0;
        }
        Ok(Matrix::from_row_major(n, n, inv)?)
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

/// Compute Givens rotation (c, s) such that `[[c, s], [-s, c]] · [x, z]ᵀ = [r, 0]ᵀ`.
fn givens(x: f64, z: f64) -> (f64, f64) {
    if z.abs() < 1e-300 {
        (1.0, 0.0)
    } else {
        let r = (x * x + z * z).sqrt();
        (x / r, z / r)
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

impl Matrix {
    /// Construct an `n × n` Hilbert matrix with entries `H[i][j] = 1/(i+j+1)`.
    pub fn hilbert(n: usize) -> Self {
        let mut data = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                data[i * n + j] = 1.0 / ((i + j + 1) as f64);
            }
        }
        Self { rows: n, cols: n, data }
    }

    /// Solve `Ax = b` with Tikhonov (L2) regularisation: minimise
    /// `||Ax - b||² + λ||x||²`, which gives `(AᵀA + λI)x = Aᵀb`.
    ///
    /// Useful for ill-conditioned systems (e.g. Hilbert matrices) where
    /// the plain solve amplifies noise.
    pub fn solve_tikhonov(&self, b: &[f64], lambda: f64) -> Result<Vec<f64>> {
        if self.rows != b.len() {
            return Err(MathError::InvalidArgument(format!(
                "solve_tikhonov: b length {} != rows {}",
                b.len(),
                self.rows
            )));
        }
        if lambda < 0.0 {
            return Err(MathError::InvalidArgument(
                "solve_tikhonov: lambda must be non-negative".into(),
            ));
        }
        let n = self.cols;
        let at = self.transpose();
        // AᵀA
        let ata = (&at * self)?;
        // Aᵀb
        let atb = at.mul_vec(b)?;
        // (AᵀA + λI)
        let mut reg = ata;
        for i in 0..n {
            reg[(i, i)] += lambda;
        }
        reg.solve(&atb)
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
        let x = fact.solve(&vec![4.0, 5.0]).unwrap();
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
}
