//! Eigenvalue solvers: Hessenberg reduction, real Schur decomposition,
//! symmetric eigen decomposition (QR algorithm), power iteration,
//! condition number and nullspace.
//!
//! Split out of the original monolithic `matrix.rs` as additional
//! inherent `impl Matrix` blocks, so the public API is unchanged.

use super::Matrix;
use crate::error::{MathError, Result};

impl Matrix {
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

impl Matrix {
    /// Condition number in the 2-norm: `σ_max / σ_min` from the singular
    /// value decomposition. Returns `inf` for singular (or zero) matrices.
    ///
    /// ```text
    /// κ(A) = ‖A‖₂ · ‖A⁻¹‖₂
    /// ```
    ///
    /// Values above ~1e12 indicate the linear system `A·x = b` is
    /// effectively ill-conditioned in `f64` arithmetic.
    pub fn condition_number(&self) -> Result<f64> {
        let svd = self.svd()?;
        let max = svd
            .singular_values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let min = svd
            .singular_values
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        if min <= 0.0 || !min.is_finite() {
            return Ok(f64::INFINITY);
        }
        Ok(max / min)
    }

    /// Orthonormal basis of the (right) nullspace `{x : A·x = 0}`, computed
    /// from the eigendecomposition of `Aᵀ·A`: eigenvectors whose eigenvalue
    /// magnitude is at most `tol` (relative to the largest eigenvalue) span
    /// the nullspace. `tol <= 0` selects the default `1e-10`.
    ///
    /// Returns one basis vector (length `cols`) per nullspace dimension,
    /// or an empty vector when `A` has full column rank.
    pub fn nullspace(&self, tol: f64) -> Result<Vec<Vec<f64>>> {
        let tol = if tol <= 0.0 { 1e-10 } else { tol };
        let at = self.transpose();
        let ata = (&at * self)?;
        let (values, vectors) = ata.symmetric_eig()?;
        let scale = values.last().copied().unwrap_or(0.0).abs().max(1.0);
        // Eigenvalues are ascending: everything below the threshold (at the
        // front) is a nullspace direction.
        let n = vectors.cols;
        let mut basis = Vec::new();
        for c in 0..n {
            if values[c].abs() > tol * scale {
                break;
            }
            let v: Vec<f64> = (0..n).map(|r| vectors[(r, c)]).collect();
            basis.push(v);
        }
        Ok(basis)
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
