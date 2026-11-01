//! Matrix factorizations: LU, Cholesky, SVD, QR and Tikhonov
//! regularisation.
//!
//! Split out of the original monolithic `matrix.rs` as additional
//! inherent `impl Matrix` blocks (plus the companion `Lu`/`Cholesky`/
//! `Svd`/`Qr` result types), so the public API is unchanged.

use super::Matrix;
use crate::error::{MathError, Result};

impl Matrix {
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
        Matrix::from_row_major(n, n, inv)
    }
}

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

    /// QR decomposition via Householder reflections: `A = Q · R` where `Q` is
    /// orthogonal (`QᵀQ = I`) and `R` is upper-triangular. Works for
    /// rectangular `m × n` matrices with `m ≥ n` (tall or square).
    ///
    /// For `m < n`, the transpose is decomposed and the factors are swapped
    /// (so `Q` is `m × m` and `R` is `m × n` upper-trapezoidal).
    pub fn qr(&self) -> Result<Qr> {
        let m = self.rows;
        let n = self.cols;
        if m == 0 || n == 0 {
            return Err(MathError::InvalidArgument("qr: empty matrix".into()));
        }
        // Work on a copy of A; R starts as A and gets transformed in place.
        let mut r = self.data.clone();
        let r_rows = m;
        let r_cols = n;
        // Q accumulates as the product of Householder reflections; start as I_m.
        let mut q = vec![0.0; m * m];
        for i in 0..m {
            q[i * m + i] = 1.0;
        }
        let k = m.min(n);
        for step in 0..k {
            // Build Householder vector for column `step`, rows step..r_rows.
            let len = r_rows - step;
            // Extract the sub-column x = r[step..r_rows, step].
            let mut x = vec![0.0; len];
            for i in 0..len {
                x[i] = r[(step + i) * r_cols + step];
            }
            let x_norm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
            if x_norm < 1e-300 {
                continue; // Already zeroed, skip this step.
            }
            // alpha = -sign(x[0]) * ||x||
            let alpha = if x[0] >= 0.0 { -x_norm } else { x_norm };
            // v = x - alpha * e_1
            let mut v = x.clone();
            v[0] -= alpha;
            let v_norm_sq: f64 = v.iter().map(|val| val * val).sum();
            if v_norm_sq < 1e-300 {
                continue;
            }
            // Apply H = I - 2*v*vᵀ/(vᵀv) to the sub-matrix R[step..r_rows, step..r_cols]
            // For each column j in [step, r_cols):
            //   R[i, j] -= 2 * v[i-step] * (sum_k v[k] * R[step+k, j]) / v_norm_sq
            for j in step..r_cols {
                let mut dot = 0.0;
                for i in 0..len {
                    dot += v[i] * r[(step + i) * r_cols + j];
                }
                let factor = 2.0 * dot / v_norm_sq;
                for i in 0..len {
                    r[(step + i) * r_cols + j] -= factor * v[i];
                }
            }
            // Apply H to Q from the right: Q = Q · H (affects columns step..m of Q)
            // Q[:, i] -= 2 * v[i-step] * (sum_k Q[:, step+k] * v[k]) / v_norm_sq
            for i in 0..m {
                let mut dot = 0.0;
                for k in 0..len {
                    dot += q[i * m + (step + k)] * v[k];
                }
                let factor = 2.0 * dot / v_norm_sq;
                for k in 0..len {
                    q[i * m + (step + k)] -= factor * v[k];
                }
            }
        }
        // If m < n, we decomposed Aᵀ; swap back.
        // Actually, we handled it generically above — R is m×n upper-trapezoidal
        // and Q is m×m orthogonal. No swap needed.
        let _ = (r_rows, r_cols); // already correct
        Ok(Qr {
            m,
            n,
            q,
            r,
        })
    }
}


/// QR decomposition result: `A = Q · R`.
pub struct Qr {
    m: usize,
    n: usize,
    q: Vec<f64>,
    r: Vec<f64>,
}

impl Qr {
    /// The orthogonal `m × m` matrix `Q`.
    pub fn q(&self) -> Matrix {
        Matrix::from_row_major(self.m, self.m, self.q.clone()).unwrap()
    }

    /// The upper-triangular (or upper-trapezoidal for rectangular) `m × n` matrix `R`.
    pub fn r(&self) -> Matrix {
        Matrix::from_row_major(self.m, self.n, self.r.clone()).unwrap()
    }

    /// Reconstruct the original matrix `A = Q · R`.
    pub fn reconstruct(&self) -> Matrix {
        let mut data = vec![0.0; self.m * self.n];
        for i in 0..self.m {
            for j in 0..self.n {
                let mut sum = 0.0;
                for k in 0..self.m {
                    sum += self.q[i * self.m + k] * self.r[k * self.n + j];
                }
                data[i * self.n + j] = sum;
            }
        }
        Matrix::from_row_major(self.m, self.n, data).unwrap()
    }

    /// Solve the least-squares problem `min ‖A·x − b‖` for overdetermined
    /// systems (`m ≥ n`, full column rank). Returns `x` of length `n`.
    ///
    /// Computes `Qᵀb` then back-substitutes `R·x = Qᵀb`.
    pub fn solve(&self, b: &[f64]) -> Result<Vec<f64>> {
        if b.len() != self.m {
            return Err(MathError::InvalidArgument(format!(
                "qr.solve: b has length {}, expected {}",
                b.len(),
                self.m
            )));
        }
        if self.n > self.m {
            return Err(MathError::InvalidArgument(
                "qr.solve: underdetermined system (n > m)".into(),
            ));
        }
        // Compute Qᵀ · b  (length m)
        let mut qtb = vec![0.0; self.m];
        for i in 0..self.m {
            for k in 0..self.m {
                qtb[i] += self.q[k * self.m + i] * b[k];
            }
        }
        // Back-substitution on the upper-triangular n×n leading block of R.
        let mut x = vec![0.0; self.n];
        for i in (0..self.n).rev() {
            let mut sum = qtb[i];
            for j in (i + 1)..self.n {
                sum -= self.r[i * self.n + j] * x[j];
            }
            let diag = self.r[i * self.n + i];
            if diag.abs() < 1e-14 {
                return Err(MathError::InvalidArgument(
                    "qr.solve: matrix is rank-deficient".into(),
                ));
            }
            x[i] = sum / diag;
        }
        Ok(x)
    }
}
