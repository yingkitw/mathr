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
}
