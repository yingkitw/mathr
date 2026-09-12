//! Sparse matrix storage (CSR/CSC) and iterative solvers.
//!
//! Compressed storage for large systems: [`Csr`] (compressed sparse row) and
//! [`Csc`] (compressed sparse column) built from coordinate triplets or dense
//! matrices, with matrix–vector products, transposition, sparse×sparse
//! multiplication, and a conjugate-gradient solver for symmetric
//! positive-definite systems.
//!
//! Matrices are kept in canonical order (row-major, ascending column indices
//! within each row); duplicate triplets are summed on construction.

use crate::error::{MathError, Result};
use crate::matrix::Matrix;

/// A triplet `(row, col, value)` in coordinate format.
pub type Triplet = (usize, usize, f64);

/// Compressed sparse row (CSR) matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    pub rows: usize,
    pub cols: usize,
    /// Row offsets, length `rows + 1`.
    pub indptr: Vec<usize>,
    /// Column indices, length `nnz`.
    pub indices: Vec<usize>,
    /// Nonzero values, length `nnz`.
    pub values: Vec<f64>,
}

/// Compressed sparse column (CSC) matrix — the transpose-oriented layout.
#[derive(Debug, Clone, PartialEq)]
pub struct Csc {
    pub rows: usize,
    pub cols: usize,
    /// Column offsets, length `cols + 1`.
    pub indptr: Vec<usize>,
    /// Row indices, length `nnz`.
    pub indices: Vec<usize>,
    /// Nonzero values, length `nnz`.
    pub values: Vec<f64>,
}

impl Csr {
    /// Build a CSR matrix from coordinate triplets.
    /// Duplicate entries at the same (row, col) are summed.
    /// Indices within each row are sorted ascending (canonical form).
    pub fn from_triplets(rows: usize, cols: usize, triplets: &[Triplet]) -> Result<Self> {
        let mut sorted: Vec<Triplet> = triplets.to_vec();
        sorted.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged: Vec<Triplet> = Vec::with_capacity(sorted.len());
        for &(r, c, v) in &sorted {
            if r >= rows {
                return Err(MathError::InvalidArgument(format!(
                    "row index {} out of bounds (rows = {})",
                    r, rows
                )));
            }
            if c >= cols {
                return Err(MathError::InvalidArgument(format!(
                    "column index {} out of bounds (cols = {})",
                    c, cols
                )));
            }
            if let Some(last) = merged.last_mut() {
                if last.0 == r && last.1 == c {
                    last.2 += v;
                    continue;
                }
            }
            merged.push((r, c, v));
        }
        let mut indptr = vec![0usize; rows + 1];
        for &(r, _, _) in &merged {
            indptr[r + 1] += 1;
        }
        for i in 0..rows {
            indptr[i + 1] += indptr[i];
        }
        let indices = merged.iter().map(|&(_, c, _)| c).collect();
        let values = merged.iter().map(|&(_, _, v)| v).collect();
        Ok(Csr {
            rows,
            cols,
            indptr,
            indices,
            values,
        })
    }

    /// Build a CSR matrix from a dense matrix, dropping exact zeros.
    pub fn from_dense(m: &Matrix) -> Self {
        let mut indptr = vec![0usize; m.rows + 1];
        let mut indices = Vec::new();
        let mut values = Vec::new();
        for i in 0..m.rows {
            for j in 0..m.cols {
                if m[(i, j)] != 0.0 {
                    indptr[i + 1] += 1;
                    indices.push(j);
                    values.push(m[(i, j)]);
                }
            }
        }
        for i in 0..m.rows {
            indptr[i + 1] += indptr[i];
        }
        Csr {
            rows: m.rows,
            cols: m.cols,
            indptr,
            indices,
            values,
        }
    }

    /// Number of stored nonzeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Convert to a dense matrix.
    pub fn to_dense(&self) -> Matrix {
        let mut data = vec![0.0; self.rows * self.cols];
        for i in 0..self.rows {
            for k in self.indptr[i]..self.indptr[i + 1] {
                data[i * self.cols + self.indices[k]] = self.values[k];
            }
        }
        Matrix::from_row_major(self.rows, self.cols, data)
            .expect("dense conversion preserves dimensions")
    }

    /// Matrix–vector product `A · x`.
    pub fn matvec(&self, x: &[f64]) -> Result<Vec<f64>> {
        if x.len() != self.cols {
            return Err(MathError::InvalidArgument(format!(
                "matvec dimension mismatch: A is {}x{}, x has length {}",
                self.rows,
                self.cols,
                x.len()
            )));
        }
        let mut y = vec![0.0; self.rows];
        for (i, yi) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for k in self.indptr[i]..self.indptr[i + 1] {
                sum += self.values[k] * x[self.indices[k]];
            }
            *yi = sum;
        }
        Ok(y)
    }

    /// Transpose: returns the same data in CSC layout.
    pub fn transpose(&self) -> Csc {
        let mut indptr = vec![0usize; self.cols + 1];
        for &c in &self.indices {
            indptr[c + 1] += 1;
        }
        for j in 0..self.cols {
            indptr[j + 1] += indptr[j];
        }
        let mut indices = vec![0usize; self.nnz()];
        let mut values = vec![0.0; self.nnz()];
        let mut next = indptr.clone();
        for i in 0..self.rows {
            for k in self.indptr[i]..self.indptr[i + 1] {
                let c = self.indices[k];
                let pos = next[c];
                indices[pos] = i;
                values[pos] = self.values[k];
                next[c] += 1;
            }
        }
        Csc {
            rows: self.rows,
            cols: self.cols,
            indptr,
            indices,
            values,
        }
    }

    /// Sparse × sparse product `A · B` (row-wise SpGEMM).
    pub fn multiply(&self, b: &Csr) -> Result<Csr> {
        if self.cols != b.rows {
            return Err(MathError::InvalidArgument(format!(
                "cannot multiply {}x{} by {}x{}",
                self.rows,
                self.cols,
                b.rows,
                b.cols
            )));
        }
        let mut indptr = vec![0usize; self.rows + 1];
        let mut indices = Vec::new();
        let mut values = Vec::new();
        let mut scratch = vec![0.0f64; b.cols];
        let mut row_cols: Vec<usize> = Vec::new();
        for i in 0..self.rows {
            row_cols.clear();
            for k in self.indptr[i]..self.indptr[i + 1] {
                let a_col = self.indices[k];
                let a_val = self.values[k];
                for l in b.indptr[a_col]..b.indptr[a_col + 1] {
                    let c = b.indices[l];
                    if !row_cols.contains(&c) {
                        row_cols.push(c);
                        scratch[c] = a_val * b.values[l];
                    } else {
                        scratch[c] += a_val * b.values[l];
                    }
                }
            }
            row_cols.sort_unstable();
            for &c in &row_cols {
                let v = scratch[c];
                if v != 0.0 {
                    indices.push(c);
                    values.push(v);
                }
                scratch[c] = 0.0;
            }
            indptr[i + 1] = indices.len();
        }
        Ok(Csr {
            rows: self.rows,
            cols: b.cols,
            indptr,
            indices,
            values,
        })
    }
}

impl Csc {
    /// Number of stored nonzeros.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Convert back to CSR layout (scatter transpose into row-major order).
    pub fn to_csr(&self) -> Csr {
        let mut indptr = vec![0usize; self.rows + 1];
        for &r in &self.indices {
            indptr[r + 1] += 1;
        }
        for i in 0..self.rows {
            indptr[i + 1] += indptr[i];
        }
        let mut indices = vec![0usize; self.nnz()];
        let mut values = vec![0.0; self.nnz()];
        let mut next = indptr.clone();
        for j in 0..self.cols {
            for k in self.indptr[j]..self.indptr[j + 1] {
                let r = self.indices[k];
                let pos = next[r];
                indices[pos] = j;
                values[pos] = self.values[k];
                next[r] += 1;
            }
        }
        Csr {
            rows: self.rows,
            cols: self.cols,
            indptr,
            indices,
            values,
        }
    }

    /// Matrix–vector product `A · x` stored column-wise.
    pub fn matvec(&self, x: &[f64]) -> Result<Vec<f64>> {
        self.to_csr().matvec(x)
    }

    /// Convert to a dense matrix.
    pub fn to_dense(&self) -> Matrix {
        self.to_csr().to_dense()
    }
}

/// Conjugate-gradient result.
#[derive(Debug, Clone)]
pub struct CgResult {
    /// Solution vector.
    pub x: Vec<f64>,
    /// Iterations performed.
    pub iterations: usize,
    /// Final residual norm `‖b − Ax‖₂`.
    pub residual: f64,
}

/// Solve `A · x = b` with the conjugate gradient method.
///
/// `A` must be symmetric positive-definite; symmetry is verified
/// structurally (values within 1e-10 relative) and indefiniteness is
/// detected as a non-positive-curvature breakdown. Iteration stops when
/// `‖r‖₂ ≤ tol · max(1, ‖b‖₂)` or `max_iter` is reached.
pub fn conjugate_gradient(a: &Csr, b: &[f64], tol: f64, max_iter: usize) -> Result<CgResult> {
    validate_cg_input(a, b, tol)?;
    cg_inner(a, b, tol, max_iter, None)
}

/// Solve `A · x = b` with Jacobi-preconditioned conjugate gradient.
///
/// The preconditioner `M = diag(A)` removes row-scaling effects: applied to
/// `(S·A·S, S·b)` it performs exactly as many iterations as plain
/// [`conjugate_gradient`] on `(A, b)` in exact arithmetic. `A` must be
/// symmetric positive-definite with a nonzero diagonal.
pub fn conjugate_gradient_jacobi(a: &Csr, b: &[f64], tol: f64, max_iter: usize) -> Result<CgResult> {
    validate_cg_input(a, b, tol)?;
    let m_inv = diagonal_inv(a)?;
    cg_inner(a, b, tol, max_iter, Some(&m_inv))
}

fn validate_cg_input(a: &Csr, b: &[f64], tol: f64) -> Result<()> {
    if a.rows != a.cols {
        return Err(MathError::InvalidArgument(format!(
            "conjugate gradient requires a square matrix, got {}x{}",
            a.rows, a.cols
        )));
    }
    if b.len() != a.rows {
        return Err(MathError::InvalidArgument(format!(
            "right-hand side has length {}, expected {}",
            b.len(),
            a.rows
        )));
    }
    if !tol.is_finite() || tol <= 0.0 {
        return Err(MathError::InvalidArgument(format!(
            "tolerance must be positive, got {}",
            tol
        )));
    }
    if !is_symmetric(a) {
        return Err(MathError::InvalidArgument(
            "matrix is not symmetric: conjugate gradient requires a symmetric positive-definite matrix"
                .into(),
        ));
    }
    Ok(())
}

fn cg_inner(
    a: &Csr,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    m_inv: Option<&[f64]>,
) -> Result<CgResult> {
    let n = a.rows;
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let b_norm = dot(b, b).sqrt();
    let threshold = tol * b_norm.max(1.0);
    let apply_pre = |v: &[f64]| -> Vec<f64> {
        match m_inv {
            Some(minv) => v.iter().zip(minv).map(|(vi, m)| vi * m).collect(),
            None => v.to_vec(),
        }
    };
    let mut z = apply_pre(&r);
    let mut p = z.clone();
    let mut rho = dot(&r, &z);

    if dot(&r, &r).sqrt() <= threshold {
        return Ok(CgResult {
            x,
            iterations: 0,
            residual: dot(&r, &r).sqrt(),
        });
    }

    for iter in 1..=max_iter {
        if rho == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "conjugate gradient stagnation at iteration {}: residual orthogonal to preconditioned residual",
                iter
            )));
        }
        let ap = a.matvec(&p)?;
        let pap = dot(&p, &ap);
        if pap <= 0.0 {
            return Err(MathError::NotConvergent(format!(
                "conjugate gradient breakdown at iteration {}: non-positive curvature (matrix is not positive-definite)",
                iter
            )));
        }
        let alpha = rho / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm <= threshold {
            return Ok(CgResult {
                x,
                iterations: iter,
                residual: r_norm,
            });
        }
        z = apply_pre(&r);
        let rho_new = dot(&r, &z);
        let beta = rho_new / rho;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rho = rho_new;
    }

    Err(MathError::NotConvergent(format!(
        "conjugate gradient did not converge in {} iterations (residual {:.3e})",
        max_iter,
        dot(&r, &r).sqrt()
    )))
}

/// BiCGStab result.
#[derive(Debug, Clone)]
pub struct BicgstabResult {
    /// Solution vector.
    pub x: Vec<f64>,
    /// Iterations performed.
    pub iterations: usize,
    /// Final residual norm `‖b − Ax‖₂`.
    pub residual: f64,
}

/// Solve `A · x = b` with Jacobi-preconditioned BiCGStab.
///
/// Works for general (including nonsymmetric) square matrices; no symmetry
/// requirement. Breakdowns (shadow-residual orthogonality, `omega = 0`,
/// zero curvature) and iteration exhaustion return `NotConvergent`.
/// Iteration stops when `‖r‖₂ ≤ tol · max(1, ‖b‖₂)`. See also
/// [`bicgstab_ilu`] for stronger ILU(0) preconditioning.
pub fn bicgstab(a: &Csr, b: &[f64], tol: f64, max_iter: usize) -> Result<BicgstabResult> {
    validate_bicgstab_input(a, b, tol)?;
    let m_inv = diagonal_inv(a)?;
    bicgstab_apply(
        a,
        b,
        tol,
        max_iter,
        move |v: &[f64]| Ok(v.iter().zip(&m_inv).map(|(vi, m)| vi * m).collect()),
    )
}

/// Solve `A · x = b` with ILU(0)-preconditioned BiCGStab.
///
/// The preconditioner `M ≈ A` is the incomplete LU factorization with zero
/// fill-in (same sparsity pattern as `A`); for matrices whose elimination
/// introduces no fill (e.g. tridiagonal) `M = A` exactly and the solve
/// converges in one iteration. Same convergence and breakdown semantics as
/// [`bicgstab`].
pub fn bicgstab_ilu(
    a: &Csr,
    b: &[f64],
    tol: f64,
    max_iter: usize,
    ilu: &Ilu0,
) -> Result<BicgstabResult> {
    validate_bicgstab_input(a, b, tol)?;
    bicgstab_apply(a, b, tol, max_iter, |v: &[f64]| ilu.solve(v))
}

fn validate_bicgstab_input(a: &Csr, b: &[f64], tol: f64) -> Result<()> {
    if a.rows != a.cols {
        return Err(MathError::InvalidArgument(format!(
            "bicgstab requires a square matrix, got {}x{}",
            a.rows, a.cols
        )));
    }
    if b.len() != a.rows {
        return Err(MathError::InvalidArgument(format!(
            "right-hand side has length {}, expected {}",
            b.len(),
            a.rows
        )));
    }
    if !tol.is_finite() || tol <= 0.0 {
        return Err(MathError::InvalidArgument(format!(
            "tolerance must be positive, got {}",
            tol
        )));
    }
    Ok(())
}

fn bicgstab_apply<P>(a: &Csr, b: &[f64], tol: f64, max_iter: usize, apply_pre: P) -> Result<BicgstabResult>
where
    P: Fn(&[f64]) -> Result<Vec<f64>>,
{
    let n = a.rows;
    let b_norm = dot(b, b).sqrt();
    let threshold = tol * b_norm.max(1.0);
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    if b_norm <= threshold {
        return Ok(BicgstabResult {
            x,
            iterations: 0,
            residual: b_norm,
        });
    }
    let rhat = b.to_vec();
    let mut v = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut rho = 1.0;
    let mut alpha = 1.0;
    let mut omega = 1.0;

    for iter in 1..=max_iter {
        let rho_new = dot(&rhat, &r);
        if rho_new == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "bicgstab breakdown at iteration {}: shadow residual orthogonality (rho = 0)",
                iter
            )));
        }
        let beta = (rho_new / rho) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        let phat = apply_pre(&p)?;
        v = a.matvec(&phat)?;
        let denom = dot(&rhat, &v);
        if denom == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "bicgstab breakdown at iteration {}: rhat·v = 0",
                iter
            )));
        }
        alpha = rho_new / denom;
        let mut s = vec![0.0; n];
        for i in 0..n {
            x[i] += alpha * phat[i];
            s[i] = r[i] - alpha * v[i];
        }
        let s_norm = dot(&s, &s).sqrt();
        if s_norm <= threshold {
            return Ok(BicgstabResult {
                x,
                iterations: iter,
                residual: s_norm,
            });
        }
        let shat = apply_pre(&s)?;
        let t = a.matvec(&shat)?;
        let tt = dot(&t, &t);
        if tt == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "bicgstab breakdown at iteration {}: A·shat = 0",
                iter
            )));
        }
        omega = dot(&t, &s) / tt;
        if omega == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "bicgstab breakdown at iteration {}: omega = 0",
                iter
            )));
        }
        for i in 0..n {
            x[i] += omega * shat[i];
            r[i] = s[i] - omega * t[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm <= threshold {
            return Ok(BicgstabResult {
                x,
                iterations: iter,
                residual: r_norm,
            });
        }
        rho = rho_new;
    }

    Err(MathError::NotConvergent(format!(
        "bicgstab did not converge in {} iterations (residual {:.3e})",
        max_iter,
        dot(&r, &r).sqrt()
    )))
}

/// ILU(0) preconditioner: incomplete LU factorization with zero fill-in.
///
/// Stores the combined `L` (strictly lower part, implicit unit diagonal) and
/// `U` (on/above diagonal, including the diagonal) factors in the sparsity
/// pattern of the original matrix. For matrices whose Gaussian elimination
/// introduces no fill (tridiagonal, diagonal), the factorization is exact.
#[derive(Debug, Clone)]
pub struct Ilu0 {
    n: usize,
    indptr: Vec<usize>,
    indices: Vec<usize>,
    values: Vec<f64>,
}

impl Ilu0 {
    /// Factorize `A` in place over its own sparsity pattern.
    ///
    /// Errors on a structurally missing diagonal entry or a zero pivot
    /// (elimination produced `u_kk = 0`).
    pub fn factorize(a: &Csr) -> Result<Self> {
        if a.rows != a.cols {
            return Err(MathError::InvalidArgument(format!(
                "ILU(0) requires a square matrix, got {}x{}",
                a.rows, a.cols
            )));
        }
        let n = a.rows;
        let mut values = a.values.clone();
        for i in 0..n {
            let (s, e) = (a.indptr[i], a.indptr[i + 1]);
            for k_idx in s..e {
                let k = a.indices[k_idx];
                if k >= i {
                    break;
                }
                let (ks, ke) = (a.indptr[k], a.indptr[k + 1]);
                let kdiag = (ks..ke)
                    .find(|&p| a.indices[p] == k)
                    .ok_or_else(|| {
                        MathError::InvalidArgument(format!(
                            "ILU(0): missing diagonal entry in row {k}"
                        ))
                    })?;
                if values[kdiag] == 0.0 {
                    return Err(MathError::InvalidArgument(format!(
                        "ILU(0): zero pivot in row {k}"
                    )));
                }
                let lk = values[k_idx] / values[kdiag];
                values[k_idx] = lk;
                if lk == 0.0 {
                    continue;
                }
                for j_idx in (k_idx + 1)..e {
                    let j = a.indices[j_idx];
                    if let Some(kp) = (ks..ke).find(|&p| a.indices[p] == j) {
                        values[j_idx] -= lk * values[kp];
                    }
                }
            }
            let diag = (s..e).find(|&p| a.indices[p] == i).ok_or_else(|| {
                MathError::InvalidArgument(format!(
                    "ILU(0): missing diagonal entry in row {i}"
                ))
            })?;
            if values[diag] == 0.0 {
                return Err(MathError::InvalidArgument(format!(
                    "ILU(0): zero pivot in row {i}"
                )));
            }
        }
        Ok(Ilu0 {
            n,
            indptr: a.indptr.clone(),
            indices: a.indices.clone(),
            values,
        })
    }

    /// Apply the preconditioner: solve `L·U·x = r`.
    pub fn solve(&self, r: &[f64]) -> Result<Vec<f64>> {
        if r.len() != self.n {
            return Err(MathError::InvalidArgument(format!(
                "right-hand side has length {}, expected {}",
                r.len(),
                self.n
            )));
        }
        let mut x = r.to_vec();
        // Forward substitution: L·z = r (unit diagonal).
        for i in 0..self.n {
            let (s, e) = (self.indptr[i], self.indptr[i + 1]);
            let mut sum = 0.0;
            for p in s..e {
                let j = self.indices[p];
                if j >= i {
                    break;
                }
                sum += self.values[p] * x[j];
            }
            x[i] -= sum;
        }
        // Back substitution: U·y = z.
        for i in (0..self.n).rev() {
            let (s, e) = (self.indptr[i], self.indptr[i + 1]);
            let diag = (s..e)
                .find(|&p| self.indices[p] == i)
                .ok_or_else(|| {
                    MathError::InvalidArgument(format!(
                        "ILU(0): missing diagonal entry in row {i}"
                    ))
                })?;
            let mut sum = 0.0;
            for p in (diag + 1)..e {
                sum += self.values[p] * x[self.indices[p]];
            }
            x[i] = (x[i] - sum) / self.values[diag];
        }
        Ok(x)
    }

    /// Number of stored entries (equals the nnz of the source pattern).
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
}

/// Reciprocal of the diagonal, erroring on any zero (or absent) diagonal entry.
fn diagonal_inv(a: &Csr) -> Result<Vec<f64>> {
    let mut m_inv = vec![0.0; a.rows];
    for (i, mi) in m_inv.iter_mut().enumerate() {
        let mut found = false;
        for k in a.indptr[i]..a.indptr[i + 1] {
            if a.indices[k] == i {
                if a.values[k] == 0.0 {
                    return Err(MathError::InvalidArgument(format!(
                        "zero diagonal entry at row {}: Jacobi preconditioner is undefined",
                        i
                    )));
                }
                *mi = 1.0 / a.values[k];
                found = true;
                break;
            }
        }
        if !found {
            return Err(MathError::InvalidArgument(format!(
                "zero diagonal entry at row {}: Jacobi preconditioner is undefined",
                i
            )));
        }
    }
    Ok(m_inv)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Structural symmetry check: every nonzero `a_ij` must have a matching
/// `a_ji` with an equal value (1e-10 relative tolerance on the diagonal scale).
fn is_symmetric(a: &Csr) -> bool {
    for i in 0..a.rows {
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            let v = a.values[k];
            if j == i {
                continue;
            }
            let mut found = false;
            for l in a.indptr[j]..a.indptr[j + 1] {
                if a.indices[l] == i {
                    let w = a.values[l];
                    if (v - w).abs() <= 1e-10 * v.abs().max(1.0) {
                        found = true;
                    }
                    break;
                }
            }
            if !found {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12 * a.abs().max(1.0)
    }

    #[test]
    fn triplets_with_duplicates_are_summed_and_sorted() {
        //  [[4, 1],  4 at (0,0) twice -> 8; (1,0)=3 sorts before (1,1)
        //   [3, 2]]
        let a = Csr::from_triplets(2, 2, &[(1, 1, 2.0), (0, 0, 4.0), (0, 0, 4.0), (1, 0, 3.0)])
            .unwrap();
        assert_eq!(a.rows, 2);
        assert_eq!(a.cols, 2);
        assert_eq!(a.nnz(), 3);
        assert_eq!(a.indptr, vec![0, 1, 3]);
        assert_eq!(a.indices, vec![0, 0, 1]);
        assert!(close(a.values[0], 8.0));
        assert!(close(a.values[1], 3.0));
        assert!(close(a.values[2], 2.0));
    }

    #[test]
    fn out_of_bounds_triplet_is_rejected() {
        assert!(Csr::from_triplets(2, 2, &[(2, 0, 1.0)]).is_err());
        assert!(Csr::from_triplets(2, 2, &[(0, 2, 1.0)]).is_err());
    }

    #[test]
    fn dense_roundtrip_preserves_values() {
        let d = Matrix::from_rows(&[vec![1.0, 0.0, 2.0], vec![0.0, 3.0, 0.0], vec![4.0, 5.0, 6.0]])
            .unwrap();
        let s = Csr::from_dense(&d);
        assert_eq!(s.nnz(), 6);
        let back = s.to_dense();
        for i in 0..3 {
            for j in 0..3 {
                assert_abs_diff_eq!(back[(i, j)], d[(i, j)]);
            }
        }
    }

    #[test]
    fn matvec_computes_row_dot_products() {
        let d = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let y = a.matvec(&[1.0, -1.0]).unwrap();
        assert_abs_diff_eq!(y[0], -1.0);
        assert_abs_diff_eq!(y[1], -1.0);
    }

    #[test]
    fn matvec_dimension_mismatch_errors() {
        let d = Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap();
        let a = Csr::from_dense(&d);
        assert!(a.matvec(&[1.0]).is_err());
        assert!(a.matvec(&[1.0, 2.0, 3.0]).is_err());
    }

    #[test]
    fn transpose_and_csc_matvec() {
        let d = Matrix::from_rows(&[vec![1.0, 2.0, 0.0], vec![0.0, 0.0, 3.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let at = a.transpose();
        assert_eq!(at.rows, 2);
        assert_eq!(at.cols, 3);
        assert_eq!(at.nnz(), 3);
        // CSC stores the SAME matrix: matvec is still A · x.
        let y = at.matvec(&[1.0, 1.0, 1.0]).unwrap();
        assert_abs_diff_eq!(y[0], 3.0);
        assert_abs_diff_eq!(y[1], 3.0);
        // Converting back restores A exactly.
        let back = at.to_csr().to_dense();
        for i in 0..2 {
            for j in 0..3 {
                assert_abs_diff_eq!(back[(i, j)], d[(i, j)]);
            }
        }
        // Interpreting the CSC buffers as row-compressed gives Aᵀ.
        let transposed = Csr {
            rows: at.cols,
            cols: at.rows,
            indptr: at.indptr.clone(),
            indices: at.indices.clone(),
            values: at.values.clone(),
        };
        let td = transposed.to_dense();
        assert_eq!(td.rows, 3);
        assert_eq!(td.cols, 2);
        for i in 0..2 {
            for j in 0..3 {
                assert_abs_diff_eq!(td[(j, i)], d[(i, j)]);
            }
        }
    }

    #[test]
    fn sparse_multiply_matches_dense_product() {
        let da = Matrix::from_rows(&[vec![1.0, 0.0, 2.0], vec![0.0, 3.0, 0.0]]).unwrap();
        let db = Matrix::from_rows(&[vec![4.0, 0.0], vec![0.0, 5.0], vec![6.0, 0.0]]).unwrap();
        let a = Csr::from_dense(&da);
        let b = Csr::from_dense(&db);
        let c = a.multiply(&b).unwrap();
        let dc = c.to_dense();
        let expected = (&da * &db).unwrap();
        assert_eq!(dc.rows, 2);
        assert_eq!(dc.cols, 2);
        for i in 0..2 {
            for j in 0..2 {
                assert_abs_diff_eq!(dc[(i, j)], expected[(i, j)]);
            }
        }
    }

    #[test]
    fn multiply_dimension_mismatch_errors() {
        let a = Csr::from_dense(&Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap());
        let b = Csr::from_dense(&Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap());
        assert!(a.multiply(&b).is_err());
    }

    #[test]
    fn poisson_system_converges_to_exact_solution() {
        // 1-D Poisson: −u'' = f with u(0)=u(1)=0, f = sin(pi x).
        // Exact solution u(x) = sin(pi x)/pi². n = 40 interior points.
        let n = 40usize;
        let h = 1.0 / (n + 1) as f64;
        let mut triplets: Vec<Triplet> = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 2.0));
            if i > 0 {
                triplets.push((i, i - 1, -1.0));
            }
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
            }
        }
        let a = Csr::from_triplets(n, n, &triplets).unwrap();
        assert_eq!(a.nnz(), 3 * n - 2);
        let b: Vec<f64> = (1..=n)
            .map(|i| (std::f64::consts::PI * i as f64 * h).sin() * h * h)
            .collect();
        let res = conjugate_gradient(&a, &b, 1e-12, 200).unwrap();
        assert!(res.iterations < 40);
        for i in 0..n {
            let exact = (std::f64::consts::PI * (i + 1) as f64 * h).sin()
                / (std::f64::consts::PI * std::f64::consts::PI);
            // 1e-4 tolerance covers the O(h²) truncation error of the 3-point
            // stencil (~4e-6 at h = 1/41); the residual bound proves the
            // DISCRETE system is solved exactly.
            assert!(
                (res.x[i] - exact).abs() < 1e-4,
                "x[{}] = {} vs {}",
                i,
                res.x[i],
                exact
            );
        }
        assert!(res.residual < 1e-11);
        // Cross-check against a direct dense solve.
        let x_dense = a.to_dense().solve(&b).unwrap();
        for i in 0..n {
            assert_abs_diff_eq!(res.x[i], x_dense[i], epsilon = 1e-9);
        }
    }

    #[test]
    fn cg_matches_dense_solve_on_spd_system() {
        let d = Matrix::from_rows(&[vec![4.0, 1.0, 0.0], vec![1.0, 3.0, 1.0], vec![0.0, 1.0, 2.0]])
            .unwrap();
        let b = [1.0, 2.0, 3.0];
        let a = Csr::from_dense(&d);
        let cg = conjugate_gradient(&a, &b, 1e-12, 100).unwrap();
        let x_dense = d.solve(&b).unwrap();
        for i in 0..3 {
            assert_abs_diff_eq!(cg.x[i], x_dense[i], epsilon = 1e-9);
        }
        let ax = a.matvec(&cg.x).unwrap();
        for i in 0..3 {
            assert_abs_diff_eq!(ax[i], b[i], epsilon = 1e-9);
        }
    }

    #[test]
    fn cg_rejects_nonsymmetric_matrix() {
        let d = Matrix::from_rows(&[vec![2.0, 1.0], vec![0.0, 2.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let err = conjugate_gradient(&a, &[1.0, 1.0], 1e-10, 50).unwrap_err();
        assert!(err.to_string().contains("not symmetric"));
    }

    #[test]
    fn cg_breaks_down_on_indefinite_matrix() {
        // Symmetric but indefinite: eigenvalues ±1. b = [1, -1] is orthogonal
        // to the null vector, so breakdown hits at the first curvature check.
        let d = Matrix::from_rows(&[vec![0.0, 1.0], vec![1.0, 0.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let err = conjugate_gradient(&a, &[1.0, -1.0], 1e-10, 50).unwrap_err();
        assert!(matches!(err, MathError::NotConvergent(_)));
    }

    #[test]
    fn cg_validates_inputs() {
        let d = Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 2.0]]).unwrap();
        let a = Csr::from_dense(&d);
        assert!(conjugate_gradient(&a, &[1.0], 1e-10, 50).is_err());
        assert!(conjugate_gradient(&a, &[1.0, 1.0], 0.0, 50).is_err());
        let rect = Csr::from_dense(&Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap());
        assert!(conjugate_gradient(&rect, &[1.0, 1.0], 1e-10, 50).is_err());
    }

    #[test]
    fn cg_zero_rhs_returns_immediately() {
        let d = Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 3.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let res = conjugate_gradient(&a, &[0.0, 0.0], 1e-10, 50).unwrap();
        assert_eq!(res.iterations, 0);
        for v in &res.x {
            assert_abs_diff_eq!(*v, 0.0);
        }
    }

    #[test]
    fn empty_matrix_roundtrips() {
        let a = Csr::from_triplets(3, 2, &[]).unwrap();
        assert_eq!(a.nnz(), 0);
        let d = a.to_dense();
        assert_eq!(d.rows, 3);
        assert_eq!(d.cols, 2);
        for i in 0..3 {
            for j in 0..2 {
                assert_abs_diff_eq!(d[(i, j)], 0.0);
            }
        }
        let at = a.transpose();
        assert_eq!(at.nnz(), 0);
    }

    fn poisson_csr(n: usize) -> Csr {
        let mut triplets: Vec<Triplet> = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 2.0));
            if i > 0 {
                triplets.push((i, i - 1, -1.0));
            }
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
            }
        }
        Csr::from_triplets(n, n, &triplets).unwrap()
    }

    #[test]
    fn jacobi_cg_matches_dense_solve() {
        // SPD (positive leading minors 4, 39, 23), badly scaled diagonal.
        let d = Matrix::from_rows(&[vec![4.0, 1.0, 0.0], vec![1.0, 10.0, 2.0], vec![0.0, 2.0, 1.0]])
            .unwrap();
        let b = [1.0, 2.0, 3.0];
        let a = Csr::from_dense(&d);
        let res = conjugate_gradient_jacobi(&a, &b, 1e-12, 100).unwrap();
        let x_dense = d.solve(&b).unwrap();
        for i in 0..3 {
            assert_abs_diff_eq!(res.x[i], x_dense[i], epsilon = 1e-9);
        }
        assert!(res.residual < 1e-11);
    }

    #[test]
    fn jacobi_cg_is_invariant_under_diagonal_scaling() {
        // Jacobi PCG on (S·A·S, S·b) does the same work as plain CG on
        // (A, b): preconditioning cancels row scaling in exact arithmetic.
        let n = 10usize;
        let a = poisson_csr(n);
        let b: Vec<f64> = (0..n).map(|i| (i as f64).sin() + 1.0).collect();
        let plain = conjugate_gradient(&a, &b, 1e-12, 200).unwrap();
        let s: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut triplets: Vec<Triplet> = Vec::new();
        for i in 0..n {
            for k in a.indptr[i]..a.indptr[i + 1] {
                let j = a.indices[k];
                triplets.push((i, j, s[i] * a.values[k] * s[j]));
            }
        }
        let a2 = Csr::from_triplets(n, n, &triplets).unwrap();
        let b2: Vec<f64> = b.iter().zip(&s).map(|(bi, si)| bi * si).collect();
        let scaled = conjugate_gradient_jacobi(&a2, &b2, 1e-12, 200).unwrap();
        assert_eq!(
            plain.iterations, scaled.iterations,
            "iteration counts must match under diagonal scaling"
        );
        // x2 = S⁻¹ x1
        let x_dense = a2.to_dense().solve(&b2).unwrap();
        for i in 0..n {
            assert_abs_diff_eq!(scaled.x[i], x_dense[i], epsilon = 1e-9);
        }
    }

    #[test]
    fn jacobi_cg_rejects_zero_diagonal() {
        // Symmetric with a zero diagonal — Jacobi preconditioner undefined.
        let d = Matrix::from_rows(&[vec![0.0, 1.0], vec![1.0, 0.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let err = conjugate_gradient_jacobi(&a, &[1.0, 1.0], 1e-10, 50).unwrap_err();
        assert!(err.to_string().contains("zero diagonal"));
    }

    #[test]
    fn bicgstab_solves_nonsymmetric_system() {
        // Unique solution of [[2,1],[0,3]]·x = [3,4] is x = [5/6, 4/3]
        // (not [1,1] — that belongs to the SPD matrix [[2,1],[1,3]]).
        let d = Matrix::from_rows(&[vec![2.0, 1.0], vec![0.0, 3.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let res = bicgstab(&a, &[3.0, 4.0], 1e-12, 50).unwrap();
        assert_abs_diff_eq!(res.x[0], 5.0 / 6.0, epsilon = 1e-9);
        assert_abs_diff_eq!(res.x[1], 4.0 / 3.0, epsilon = 1e-9);
        let x_dense = d.solve(&[3.0, 4.0]).unwrap();
        for i in 0..2 {
            assert_abs_diff_eq!(res.x[i], x_dense[i], epsilon = 1e-9);
        }
        assert!(res.residual < 1e-11);
    }

    #[test]
    fn bicgstab_convection_diffusion_matches_dense_solve() {
        // Nonsymmetric convection-diffusion stencil:
        // sub = -h⁻² - (2h)⁻¹, diag = 2h⁻², super = -h⁻² + (2h)⁻¹
        let n = 30usize;
        let ih = (n + 1) as f64;
        let mut triplets: Vec<Triplet> = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 2.0 * ih * ih));
            if i > 0 {
                triplets.push((i, i - 1, -ih * ih - ih / 2.0));
            }
            if i + 1 < n {
                triplets.push((i, i + 1, -ih * ih + ih / 2.0));
            }
        }
        let a = Csr::from_triplets(n, n, &triplets).unwrap();
        let b: Vec<f64> = (0..n).map(|i| ((i + 1) as f64).sqrt()).collect();
        let res = bicgstab(&a, &b, 1e-12, 200).unwrap();
        let x_dense = a.to_dense().solve(&b).unwrap();
        for i in 0..n {
            assert_abs_diff_eq!(res.x[i], x_dense[i], epsilon = 1e-7);
        }
        assert!(res.residual < 1e-10);
    }

    #[test]
    fn bicgstab_breaks_down_on_singular_inconsistent_system() {
        // A = [[1,1],[1,1]] singular, b = [1,2] not in the column space:
        // the residual is driven into the null space and rhat·v hits 0.
        let d = Matrix::from_rows(&[vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let err = bicgstab(&a, &[1.0, 2.0], 1e-10, 50).unwrap_err();
        assert!(matches!(err, MathError::NotConvergent(_)));
    }

    #[test]
    fn bicgstab_rejects_zero_diagonal_and_bad_inputs() {
        let d = Matrix::from_rows(&[vec![0.0, 1.0], vec![1.0, 2.0]]).unwrap();
        let a = Csr::from_dense(&d);
        assert!(bicgstab(&a, &[1.0, 1.0], 1e-10, 50)
            .unwrap_err()
            .to_string()
            .contains("zero diagonal"));
        let good = Csr::from_dense(&Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 3.0]]).unwrap());
        assert!(bicgstab(&good, &[1.0], 1e-10, 50).is_err());
        assert!(bicgstab(&good, &[1.0, 1.0], 0.0, 50).is_err());
        let rect = Csr::from_dense(&Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap());
        assert!(bicgstab(&rect, &[1.0, 1.0], 1e-10, 50).is_err());
    }

    #[test]
    fn bicgstab_zero_rhs_returns_immediately() {
        let good = Csr::from_dense(&Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 3.0]]).unwrap());
        let res = bicgstab(&good, &[0.0, 0.0], 1e-10, 50).unwrap();
        assert_eq!(res.iterations, 0);
        for v in &res.x {
            assert_abs_diff_eq!(*v, 0.0);
        }
    }

    #[test]
    fn bicgstab_solves_spd_system_like_cg() {
        let a = poisson_csr(20);
        let b: Vec<f64> = (0..20).map(|i| (i + 1) as f64).collect();
        let cg = conjugate_gradient(&a, &b, 1e-12, 200).unwrap();
        let bg = bicgstab(&a, &b, 1e-12, 200).unwrap();
        for i in 0..20 {
            assert_abs_diff_eq!(bg.x[i], cg.x[i], epsilon = 1e-8);
        }
    }

    fn convection_diffusion_csr(n: usize) -> Csr {
        // Nonsymmetric central-difference convection-diffusion stencil.
        let ih = (n + 1) as f64;
        let mut triplets: Vec<Triplet> = Vec::new();
        for i in 0..n {
            triplets.push((i, i, 2.0 * ih * ih));
            if i > 0 {
                triplets.push((i, i - 1, -ih * ih - ih / 2.0));
            }
            if i + 1 < n {
                triplets.push((i, i + 1, -ih * ih + ih / 2.0));
            }
        }
        Csr::from_triplets(n, n, &triplets).unwrap()
    }

    #[test]
    fn ilu0_factorization_preserves_pattern() {
        let a = convection_diffusion_csr(20);
        let ilu = Ilu0::factorize(&a).unwrap();
        assert_eq!(ilu.nnz(), a.nnz());
        assert_eq!(ilu.indptr, a.indptr);
        assert_eq!(ilu.indices, a.indices);
    }

    #[test]
    fn ilu0_is_exact_for_tridiagonal() {
        // Tridiagonal elimination fills nothing: M = A exactly, so
        // solving M·x = r recovers A⁻¹r for any r.
        let a = poisson_csr(12);
        let ilu = Ilu0::factorize(&a).unwrap();
        let r: Vec<f64> = (0..12).map(|i| (i as f64).cos() + 1.0).collect();
        let x = ilu.solve(&r).unwrap();
        let ax = a.matvec(&x).unwrap();
        for i in 0..12 {
            assert_abs_diff_eq!(ax[i], r[i], epsilon = 1e-9);
        }
    }

    #[test]
    fn ilu0_errors_on_zero_pivot_and_missing_diagonal() {
        // Elimination zeroes row 1's diagonal: [[1,1],[1,1]].
        let d = Matrix::from_rows(&[vec![1.0, 1.0], vec![1.0, 1.0]]).unwrap();
        let err = Ilu0::factorize(&Csr::from_dense(&d)).unwrap_err();
        assert!(err.to_string().contains("zero pivot"), "{err}");
        // Structurally missing diagonal.
        let d2 = Matrix::from_rows(&[vec![0.0, 1.0], vec![1.0, 2.0]]).unwrap();
        let err2 = Ilu0::factorize(&Csr::from_dense(&d2)).unwrap_err();
        assert!(err2.to_string().contains("missing diagonal"), "{err2}");
        // Non-square input.
        let rect = Csr::from_dense(&Matrix::from_rows(&[vec![1.0, 2.0]]).unwrap());
        assert!(Ilu0::factorize(&rect).is_err());
        // Solve dimension check.
        let ok = Ilu0::factorize(&poisson_csr(3)).unwrap();
        assert!(ok.solve(&[1.0, 2.0]).is_err());
    }

    #[test]
    fn bicgstab_ilu_matches_dense_solve() {
        let a = convection_diffusion_csr(30);
        let b: Vec<f64> = (0..30).map(|i| ((i + 1) as f64).sqrt()).collect();
        let ilu = Ilu0::factorize(&a).unwrap();
        let res = bicgstab_ilu(&a, &b, 1e-12, 200, &ilu).unwrap();
        let x_dense = a.to_dense().solve(&b).unwrap();
        for i in 0..30 {
            assert_abs_diff_eq!(res.x[i], x_dense[i], epsilon = 1e-7);
        }
        assert!(res.residual < 1e-10);
    }

    #[test]
    fn bicgstab_ilu_converges_in_one_iteration_on_diagonal_matrix() {
        // M = A exactly: the preconditioned system is solved immediately.
        let d = Matrix::from_rows(&[vec![2.0, 0.0], vec![0.0, 3.0]]).unwrap();
        let a = Csr::from_dense(&d);
        let ilu = Ilu0::factorize(&a).unwrap();
        let res = bicgstab_ilu(&a, &[2.0, 3.0], 1e-12, 50, &ilu).unwrap();
        assert!(res.iterations <= 1, "iterations = {}", res.iterations);
    }

    #[test]
    fn bicgstab_ilu_beats_jacobi_on_convection_diffusion() {
        let a = convection_diffusion_csr(30);
        let b: Vec<f64> = (0..30).map(|i| (i as f64).sin() + 1.0).collect();
        let jacobi = bicgstab(&a, &b, 1e-12, 500).unwrap();
        let ilu = Ilu0::factorize(&a).unwrap();
        let prec = bicgstab_ilu(&a, &b, 1e-12, 500, &ilu).unwrap();
        assert!(
            prec.iterations < jacobi.iterations,
            "ILU(0) iterations = {}, Jacobi iterations = {}",
            prec.iterations,
            jacobi.iterations
        );
        for i in 0..30 {
            assert_abs_diff_eq!(prec.x[i], jacobi.x[i], epsilon = 1e-7);
        }
    }
}
