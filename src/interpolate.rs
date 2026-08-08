//! Polynomial interpolation from scratch.
//!
//! Provides Lagrange interpolation and Newton's divided-difference
//! interpolation for fitting a polynomial through a set of points.

use crate::error::{MathError, Result};

/// Lagrange interpolation: given points `(x_i, y_i)`, returns the
/// interpolated value at `x`.
pub fn lagrange_interp(points: &[(f64, f64)], x: f64) -> Result<f64> {
    if points.is_empty() {
        return Err(MathError::InvalidArgument("lagrange_interp: empty points".into()));
    }
    if points.len() == 1 {
        return Ok(points[0].1);
    }
    let n = points.len();
    let mut result = 0.0;
    for i in 0..n {
        let (xi, yi) = points[i];
        let mut term = yi;
        for j in 0..n {
            if i == j {
                continue;
            }
            let (xj, _) = points[j];
            let denom = xi - xj;
            if denom.abs() < 1e-14 {
                return Err(MathError::InvalidArgument(
                    "lagrange_interp: duplicate x values".into(),
                ));
            }
            term *= (x - xj) / denom;
        }
        result += term;
    }
    Ok(result)
}

/// Newton's divided-difference interpolation.
///
/// Builds the coefficient table once, then allows efficient evaluation
/// at multiple points via Newton's form.
pub struct NewtonInterpolator {
    coeffs: Vec<f64>,
    xs: Vec<f64>,
}

impl NewtonInterpolator {
    /// Build from a set of `(x, y)` points.
    pub fn new(points: &[(f64, f64)]) -> Result<Self> {
        if points.is_empty() {
            return Err(MathError::InvalidArgument("NewtonInterpolator: empty points".into()));
        }
        let n = points.len();
        let xs: Vec<f64> = points.iter().map(|(x, _)| *x).collect();
        let mut table = vec![vec![0.0; n]; n];
        for i in 0..n {
            table[i][0] = points[i].1;
        }
        for j in 1..n {
            for i in 0..(n - j) {
                let denom = xs[i] - xs[i + j];
                if denom.abs() < 1e-14 {
                    return Err(MathError::InvalidArgument(
                        "NewtonInterpolator: duplicate x values".into(),
                    ));
                }
                table[i][j] = (table[i][j - 1] - table[i + 1][j - 1]) / denom;
            }
        }
        let coeffs: Vec<f64> = (0..n).map(|j| table[0][j]).collect();
        Ok(Self { coeffs, xs })
    }

    /// Evaluate the interpolating polynomial at `x` using Horner-like form.
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.coeffs.len();
        let mut result = self.coeffs[n - 1];
        for i in (0..(n - 1)).rev() {
            result = result * (x - self.xs[i]) + self.coeffs[i];
        }
        result
    }

    /// Return the divided-difference coefficients.
    pub fn coefficients(&self) -> &[f64] {
        &self.coeffs
    }
}

/// Convenience: Newton interpolation at a single point.
pub fn newton_interp(points: &[(f64, f64)], x: f64) -> Result<f64> {
    let interp = NewtonInterpolator::new(points)?;
    Ok(interp.eval(x))
}

/// Linear interpolation between two points (a special case of Lagrange).
pub fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// A natural cubic spline interpolant.
///
/// Stores per-segment coefficients for each piece `S_i(x) = a_i + b_i·dx +
/// c_i·dx² + d_i·dx³` where `dx = x - x_i`. Boundary condition: second
/// derivative is zero at both ends.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    xs: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl CubicSpline {
    /// Build a natural cubic spline (S''(x_0) = S''(x_n) = 0) through the
    /// given points. The points do not need to be uniformly spaced, but
    /// they must be sorted in strictly increasing `x` order.
    pub fn new(points: &[(f64, f64)]) -> Result<Self> {
        Self::build(points, None, None)
    }

    /// Build a cubic spline with specified first derivatives at the
    /// endpoints (a "clamped" spline).
    pub fn clamped(points: &[(f64, f64)], first_deriv: f64, last_deriv: f64) -> Result<Self> {
        Self::build(points, Some(first_deriv), Some(last_deriv))
    }

    fn build(points: &[(f64, f64)], first_deriv: Option<f64>, last_deriv: Option<f64>) -> Result<Self> {
        let n = points.len();
        if n < 3 {
            return Err(MathError::InvalidArgument(
                "cubic spline: need at least 3 points".into(),
            ));
        }
        let mut xs: Vec<f64> = points.iter().map(|(x, _)| *x).collect();
        let ys: Vec<f64> = points.iter().map(|(_, y)| *y).collect();
        // Make sure x values are sorted ascending — but only error on duplicates;
        // silently sort the rest (still passes through every point).
        let mut indices: Vec<usize> = (0..n).collect();
        indices.sort_by(|&i, &j| xs[i].partial_cmp(&xs[j]).unwrap_or(std::cmp::Ordering::Equal));
        let mut sorted_xs = Vec::with_capacity(n);
        let mut sorted_ys = Vec::with_capacity(n);
        for &i in &indices {
            sorted_xs.push(xs[i]);
            sorted_ys.push(ys[i]);
        }
        for i in 1..n {
            if (sorted_xs[i] - sorted_xs[i - 1]).abs() < 1e-14 {
                return Err(MathError::InvalidArgument(
                    "cubic spline: duplicate x values".into(),
                ));
            }
        }
        xs = sorted_xs;
        let ys = sorted_ys;

        let h: Vec<f64> = (0..(n - 1)).map(|i| xs[i + 1] - xs[i]).collect();

        // Solve the tridiagonal system for second derivatives M_i.
        let mut m = vec![0.0_f64; n];

        match (first_deriv, last_deriv) {
            (None, None) => {
                // Natural boundary: M_0 = M_{n-1} = 0. Reduces to system of
                // size n-2 for unknowns M_1 .. M_{n-2}.
                let m_size = n - 2;
                if m_size == 0 {
                    // Two knots: trivial linear interpolation (every a = y_i).
                    let nseg = n - 1;
                    let mut a = vec![0.0; nseg];
                    let mut b = vec![0.0; nseg];
                    let c = vec![0.0; nseg];
                    let d = vec![0.0; nseg];
                    for i in 0..nseg {
                        a[i] = ys[i];
                        b[i] = (ys[i + 1] - ys[i]) / h[i];
                    }
                    return Ok(Self { xs, a, b, c, d });
                }
                let mut sub = vec![0.0_f64; m_size - 1];
                let mut diag = vec![0.0_f64; m_size];
                let mut sup = vec![0.0_f64; m_size - 1];
                let mut rhs = vec![0.0_f64; m_size];
                for k in 0..m_size {
                    // k corresponds to original index i = k+1 (i.e. M_{k+1}).
                    diag[k] = 2.0 * (h[k] + h[k + 1]);
                    if k + 1 < m_size {
                        sub[k] = h[k + 1]; // coeff of M_{k+2} in row k (shifted by M_{k+1}-1 = M_{k+1})
                    }
                    if k > 0 {
                        sup[k - 1] = h[k]; // coeff of M_k+1 in row k (shifted by M_{k+1})
                    }
                    rhs[k] = 6.0 * ((ys[k + 2] - ys[k + 1]) / h[k + 1] - (ys[k + 1] - ys[k]) / h[k]);
                }
                let sol = solve_tridiag(&sub, &diag, &sup, &rhs)?;
                for (k, &v) in sol.iter().enumerate() {
                    m[k + 1] = v;
                }
            }
            (Some(fp), Some(fpp)) => {
                // Clamped boundary: full tridiagonal system of size n.
                let m_size = n;
                let mut sub = vec![0.0_f64; m_size - 1];
                let mut diag = vec![0.0_f64; m_size];
                let mut sup = vec![0.0_f64; m_size - 1];
                let mut rhs = vec![0.0_f64; m_size];
                diag[0] = 2.0 * h[0];
                sup[0] = h[0];
                rhs[0] = 6.0 * ((ys[1] - ys[0]) / h[0] - fp);
                for k in 1..(m_size - 1) {
                    sub[k - 1] = h[k - 1];
                    diag[k] = 2.0 * (h[k - 1] + h[k]);
                    sup[k] = h[k];
                    rhs[k] = 6.0 * ((ys[k + 1] - ys[k]) / h[k] - (ys[k] - ys[k - 1]) / h[k - 1]);
                }
                sub[m_size - 2] = h[m_size - 2];
                diag[m_size - 1] = 2.0 * h[m_size - 2];
                rhs[m_size - 1] = 6.0 * (fpp - (ys[m_size - 1] - ys[m_size - 2]) / h[m_size - 2]);
                m = solve_tridiag(&sub, &diag, &sup, &rhs)?;
            }
            _ => {
                return Err(MathError::InvalidArgument(
                    "cubic spline: both endpoint derivatives must be given for clamped".into(),
                ));
            }
        }

        let nseg = n - 1;
        let mut a = vec![0.0; nseg];
        let mut b = vec![0.0; nseg];
        let mut c = vec![0.0; nseg];
        let mut d = vec![0.0; nseg];
        for i in 0..nseg {
            a[i] = ys[i];
            b[i] = (ys[i + 1] - ys[i]) / h[i] - h[i] * m[i] / 3.0 - h[i] * m[i + 1] / 6.0;
            c[i] = m[i] / 2.0;
            d[i] = (m[i + 1] - m[i]) / (6.0 * h[i]);
        }
        Ok(Self { xs, a, b, c, d })
    }

    /// Evaluate the spline at `x`.
    pub fn eval(&self, x: f64) -> f64 {
        let (i, dx) = self.locate(x);
        self.a[i] + self.b[i] * dx + self.c[i] * dx * dx + self.d[i] * dx * dx * dx
    }

    /// Evaluate the first derivative at `x`.
    pub fn derivative(&self, x: f64) -> f64 {
        let (i, dx) = self.locate(x);
        self.b[i] + 2.0 * self.c[i] * dx + 3.0 * self.d[i] * dx * dx
    }

    /// Evaluate the second derivative at `x`.
    pub fn second_derivative(&self, x: f64) -> f64 {
        let (i, dx) = self.locate(x);
        2.0 * self.c[i] + 6.0 * self.d[i] * dx
    }

    /// The knots (sorted x values).
    pub fn knots(&self) -> &[f64] { &self.xs }

    /// Locate the segment containing `x`. Returns `(segment_index, dx = x - x_i)`.
    fn locate(&self, x: f64) -> (usize, f64) {
        let n = self.xs.len();
        if x <= self.xs[0] {
            return (0, x - self.xs[0]);
        }
        if x >= self.xs[n - 1] {
            return (n - 2, x - self.xs[n - 2]);
        }
        // Binary search for the segment.
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.xs[mid] > x {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        (lo, x - self.xs[lo])
    }
}

/// Solve a tridiagonal system T·x = r with sub-diagonal `sub` (length m-1),
/// main diagonal `diag` (length m), and super-diagonal `sup` (length m-1),
/// where `m` is the system size.
fn solve_tridiag(sub: &[f64], diag: &[f64], sup: &[f64], rhs: &[f64]) -> Result<Vec<f64>> {
    let n = diag.len();
    if sub.len() + 1 != n || sup.len() + 1 != n || rhs.len() != n {
        return Err(MathError::InvalidArgument("solve_tridiag: length mismatch".into()));
    }
    let c = sup.to_vec();
    let mut d = rhs.to_vec();
    let mut bb = diag.to_vec();
    if bb[0].abs() < 1e-14 {
        return Err(MathError::InvalidArgument("solve_tridiag: zero pivot".into()));
    }
    for i in 1..n {
        let m = sub[i - 1] / bb[i - 1];
        bb[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / bb[n - 1];
    for i in (0..(n - 1)).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / bb[i];
    }
    Ok(x)
}

/// Chebyshev polynomial of the first kind, T_n(x), evaluated via the
/// recurrence `T_0 = 1`, `T_1 = x`, `T_{n+1} = 2x · T_n − T_{n−1}`.
pub fn chebyshev_t(n: u32, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut t_prev = 1.0_f64; // T_0
    let mut t_curr = x; // T_1
    for _ in 2..=n {
        let t_next = 2.0 * x * t_curr - t_prev;
        t_prev = t_curr;
        t_curr = t_next;
    }
    t_curr
}

/// Chebyshev nodes on [-1, 1]: `x_k = cos((2k − 1)π / (2n))`, k = 1..=n.
/// These are the optimal interpolation points that minimise the Runge
/// phenomenon.
pub fn chebyshev_nodes(n: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(n);
    let nf = n as f64;
    for k in 1..=n {
        let x = ((2 * k - 1) as f64) * std::f64::consts::PI / (2.0 * nf);
        out.push(x.cos());
    }
    out
}

/// Approximate a function on `[-1, 1]` by a truncated Chebyshev series
/// `f(x) ≈ c_0 T_0(x) + c_1 T_1(x) + ⋯ + c_n T_n(x)` with `n + 1`
/// coefficients computed by sampling the function at the Chebyshev nodes.
///
/// Returns the coefficient vector `c` whose `k`-th entry is the
/// Chebyshev coefficient for `T_k(x)`.
pub fn chebyshev_coefficients<F: Fn(f64) -> f64>(f: F, n: usize) -> Vec<f64> {
    let nodes = chebyshev_nodes(n + 1);
    let samples: Vec<f64> = nodes.iter().map(|&x| f(x)).collect();
    // Discrete Chebyshev transform (type-1 DCT):
    //   c_k = (2 − δ_{k0}) / (n+1) · Σ_{j=0}^n f(x_j) · cos(k·θ_j),
    //   θ_j = (2j + 1)π / (2(n+1))
    let mut out = vec![0.0; n + 1];
    let np1 = (n + 1) as f64;
    for k in 0..=n {
        let mut s = 0.0;
        for j in 0..=n {
            let theta = (2 * j + 1) as f64 * std::f64::consts::PI / (2.0 * np1);
            s += samples[j] * (k as f64 * theta).cos();
        }
        let factor = if k == 0 { 1.0 } else { 2.0 };
        out[k] = factor / np1 * s;
    }
    out
}

/// Evaluate a Chebyshev series `Σ c_k T_k(x)` at a point `x ∈ [-1, 1]`.
/// Uses Clenshaw's recurrence for numerical stability.
pub fn chebyshev_eval(coeffs: &[f64], x: f64) -> f64 {
    // Clenshaw: b_n = 0, b_{n+1} = 0;
    //           b_k = c_k + 2x·b_{k+1} − b_{k+2};
    // value = b_0 − x·b_1.
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }
    let mut b_next = 0.0_f64;
    let mut b_curr = 0.0_f64;
    for k in (0..n).rev() {
        let b_prev = coeffs[k] + 2.0 * x * b_curr - b_next;
        b_next = b_curr;
        b_curr = b_prev;
    }
    b_curr - x * b_next
}

/// Map Chebyshev coefficients computed on `[-1, 1]` to an interval `[a, b]`.
/// `f(x) ≈ Σ c_k T_k((2x − a − b) / (b − a))`, so evaluate at the rescaled
/// argument before passing to [`chebyshev_eval`].
pub fn chebyshev_rescale(x: f64, a: f64, b: f64) -> f64 {
    (2.0 * x - a - b) / (b - a)
}

/// Legendre polynomial of the first kind, `P_n(x)`, evaluated via the
/// three-term recurrence
/// `(n + 1) P_{n+1}(x) = (2n + 1) x P_n(x) − n P_{n−1}(x)`.
pub fn legendre_p(n: u32, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p_prev = 1.0_f64; // P_0
    let mut p_curr = x; // P_1
    for k in 1..n {
        // p_next = ((2k+1) x p_curr − k p_prev) / (k + 1)
        let kf = k as f64;
        let p_next = ((2.0 * kf + 1.0) * x * p_curr - kf * p_prev) / (kf + 1.0);
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}

/// Associated Legendre function of integer degree `n` and order `m`,
/// `P_n^m(x) = (−1)^m (1 − x²)^{m/2} · d^m P_n(x) / dx^m`,
/// with the convention `P_n^m = 0` for `m > n`.
///
/// Built from the closed-form seed `P_m^m(x)` using the upward recurrence
/// `(n − m + 1) · P_{n+1}^m = (2n + 1) · x · P_n^m − (n + m) · P_{n−1}^m`.
pub fn legendre_associated(n: u32, m: u32, x: f64) -> f64 {
    if m > n {
        return 0.0;
    }
    if m == 0 {
        return legendre_p(n, x);
    }
    let pmm = legendre_pmm(m, x);
    if n == m {
        return pmm;
    }
    let mut p_prev = pmm;
    let mut p_curr = (2.0 * m as f64 + 1.0) * x * pmm;
    if n == m + 1 {
        return p_curr;
    }
    for k in (m + 1)..n {
        let kf = k as f64;
        let mf = m as f64;
        let p_next = ((2.0 * kf + 1.0) * x * p_curr - (kf + mf) * p_prev) / (kf - mf + 1.0);
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}

/// Closed form for `P_m^m(x) = (−1)^m · (2m − 1)!! · (1 − x²)^{m/2}`.
fn legendre_pmm(m: u32, x: f64) -> f64 {
    let mut val = 1.0_f64;
    for k in 1..=m {
        val *= -((2 * k - 1) as f64);
    }
    val * (1.0 - x * x).powf(m as f64 / 2.0)
}

/// Gauss–Legendre quadrature: returns the `n` nodes and weights that
/// exactly integrate polynomials of degree `≤ 2n − 1` over `[-1, 1]`.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let n = n.max(1);
    if n == 1 {
        return (vec![0.0], vec![2.0]);
    }
    // Build the n×n symmetric Jacobi matrix in flat row-major storage.
    let mut jm = vec![0.0_f64; n * n];
    for k in 1..n {
        let kf = k as f64;
        let beta = kf / ((2.0 * kf + 1.0) * (2.0 * kf - 1.0)).sqrt();
        jm[(k - 1) * n + k] = beta;
        jm[k * n + (k - 1)] = beta;
    }
    // Compute eigendecomposition via Jacobi rotation; track V.
    let mut v = vec![0.0_f64; n * n];
    for i in 0..n { v[i * n + i] = 1.0; }
    let max_sweeps = 200usize;
    let tol = 1e-15;
    for _ in 0..max_sweeps {
        let mut p = 0usize;
        let mut q = 1;
        let mut max_off = 0.0;
        for r in 0..n {
            for c in (r + 1)..n {
                let val = jm[r * n + c].abs();
                if val > max_off { max_off = val; p = r; q = c; }
            }
        }
        if max_off < tol { break; }
        let app = jm[p * n + p];
        let aqq = jm[q * n + q];
        let apq = jm[p * n + q];
        let theta = if (app - aqq).abs() < 1e-30 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * ((2.0 * apq) / (aqq - app)).atan()
        };
        let (s, c) = theta.sin_cos();
        for i in 0..n {
            if i != p && i != q {
                let bip = jm[i * n + p];
                let biq = jm[i * n + q];
                jm[i * n + p] = c * bip - s * biq;
                jm[p * n + i] = jm[i * n + p];
                jm[i * n + q] = s * bip + c * biq;
                jm[q * n + i] = jm[i * n + q];
            }
        }
        let bpp = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        let bqq = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        jm[p * n + p] = bpp;
        jm[q * n + q] = bqq;
        jm[p * n + q] = 0.0;
        jm[q * n + p] = 0.0;
        for i in 0..n {
            let vip = v[i * n + p];
            let viq = v[i * n + q];
            v[i * n + p] = c * vip - s * viq;
            v[i * n + q] = s * vip + c * viq;
        }
    }
    let mut eig_pairs: Vec<(f64, usize)> = (0..n).map(|i| (jm[i * n + i], i)).collect();
    eig_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let nodes: Vec<f64> = eig_pairs.iter().map(|(v, _)| *v).collect();
    let mut weights = vec![0.0_f64; n];
    for (new_i, (_, old_i)) in eig_pairs.iter().enumerate() {
        let v0 = v[*old_i];
        weights[new_i] = 2.0 * v0 * v0;
    }
    (nodes, weights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn lagrange_through_points() {
        let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0)];
        // All points should be recovered exactly
        for &(x, y) in &pts {
            assert!(close(lagrange_interp(&pts, x).unwrap(), y, 1e-10));
        }
        // f(x) = x^2 + 1, so f(0.5) = 1.25
        assert!(close(lagrange_interp(&pts, 0.5).unwrap(), 1.25, 1e-10));
    }

    #[test]
    fn newton_through_points() {
        let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0), (3.0, 10.0)];
        let interp = NewtonInterpolator::new(&pts).unwrap();
        for &(x, y) in &pts {
            assert!(close(interp.eval(x), y, 1e-10));
        }
        // f(x) = x^2 + 1, so f(1.5) = 3.25
        assert!(close(interp.eval(1.5), 3.25, 1e-10));
    }

    #[test]
    fn lagrange_and_newton_agree() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 8.0), (3.0, 27.0)];
        for &x in &[0.5, 1.5, 2.5, 0.1, 2.9] {
            let l = lagrange_interp(&pts, x).unwrap();
            let n = newton_interp(&pts, x).unwrap();
            assert!(close(l, n, 1e-10), "at x={}: lagrange={} newton={}", x, l, n);
        }
    }

    #[test]
    fn single_point() {
        let pts = vec![(5.0, 42.0)];
        assert!(close(lagrange_interp(&pts, 100.0).unwrap(), 42.0, 1e-10));
    }

    #[test]
    fn lerp_basic() {
        assert!(close(lerp(0.0, 0.0, 10.0, 100.0, 5.0), 50.0, 1e-10));
        assert!(close(lerp(0.0, 0.0, 10.0, 100.0, 2.5), 25.0, 1e-10));
    }

    #[test]
    fn duplicate_x_errors() {
        let pts = vec![(1.0, 2.0), (1.0, 3.0)];
        assert!(lagrange_interp(&pts, 0.5).is_err());
        assert!(NewtonInterpolator::new(&pts).is_err());
    }

    #[test]
    fn cubic_spline_through_points() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];
        let sp = CubicSpline::new(&pts).unwrap();
        // Spline interpolates through every knot exactly.
        for &(x, y) in &pts {
            assert!(close(sp.eval(x), y, 1e-10), "x={}: got {} want {}", x, sp.eval(x), y);
        }
    }

    #[test]
    fn cubic_spline_smooth_join() {
        // C^2 continuous: first derivative is continuous across segment boundaries.
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, 1.0), (4.0, 0.0)];
        let sp = CubicSpline::new(&pts).unwrap();
        for i in 1..(pts.len() - 1) {
            let x = pts[i].0;
            let left = sp.derivative(x - 1e-7);
            let right = sp.derivative(x + 1e-7);
            assert!((left - right).abs() < 1e-3, "at {}: {} vs {}", x, left, right);
        }
    }

    #[test]
    fn cubic_spline_natural_second_deriv() {
        // Natural spline: S''(x_0) = S''(x_n) = 0.
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];
        let sp = CubicSpline::new(&pts).unwrap();
        assert!(sp.second_derivative(pts[0].0).abs() < 1e-9);
        assert!(sp.second_derivative(pts[pts.len() - 1].0).abs() < 1e-9);
    }

    #[test]
    fn cubic_spline_extrapolates_naturally() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];
        let sp = CubicSpline::new(&pts).unwrap();
        let v = sp.eval(1.5);
        // Knot at 1.5 isn't given; should be a smooth in-between value.
        assert!(v > 0.0 && v < 9.0, "got {}", v);
    }

    #[test]
    fn cubic_spline_clamped() {
        // A clamped spline with derivative constraints picks a curve that's
        // closer to the derivative-matching function.
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];
        let sp = CubicSpline::clamped(&pts, 0.0, 6.0).unwrap();
        for &(x, y) in &pts {
            assert!(close(sp.eval(x), y, 1e-10));
        }
        assert!((sp.derivative(pts[0].0) - 0.0).abs() < 1e-9);
        assert!((sp.derivative(pts[pts.len() - 1].0) - 6.0).abs() < 1e-9);
    }

    #[test]
    fn cubic_spline_reproduces_cubic() {
        // A cubic spline through four samples of a cubic function should
        // reproduce the cubic exactly (because n=4 samples of a degree-3
        // polynomial = the polynomial itself for natural spline with matching
        // values; clamped is exact when boundary derivatives match too).
        let f = |x: f64| x * x * x - 2.0 * x * x + x - 1.0;
        let fp = |x: f64| 3.0 * x * x - 4.0 * x + 1.0;
        let pts: Vec<(f64, f64)> = (0..=4).map(|i| {
            let x = i as f64;
            (x, f(x))
        }).collect();
        let sp = CubicSpline::clamped(&pts, fp(0.0), fp(4.0)).unwrap();
        for &x in &[0.5, 1.3, 2.7, 3.5] {
            assert!(close(sp.eval(x), f(x), 1e-8), "x={}: {} vs {}", x, sp.eval(x), f(x));
        }
    }

    #[test]
    fn cubic_spline_duplicate_x_errors() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (1.0, 2.0)];
        assert!(CubicSpline::new(&pts).is_err());
    }

    #[test]
    fn chebyshev_recurrence() {
        // T_0(x) = 1
        assert!(close(chebyshev_t(0, 0.5), 1.0, 1e-15));
        // T_1(x) = x
        assert!(close(chebyshev_t(1, 0.5), 0.5, 1e-15));
        // T_2(x) = 2x² − 1
        assert!(close(chebyshev_t(2, 0.5), 2.0 * 0.25 - 1.0, 1e-15));
        // T_3(x) = 4x³ − 3x
        assert!(close(chebyshev_t(3, 0.5), 4.0 * 0.125 - 3.0 * 0.5, 1e-15));
        // Chebyshev at the special point T_n(cos θ) = cos(nθ)
        let theta: f64 = 0.7;
        let x = theta.cos();
        assert!(close(chebyshev_t(5, x), (5.0 * theta).cos(), 1e-10));
    }

    #[test]
    fn chebyshev_nodes_in_range() {
        for n in [3, 5, 10, 32] {
            let nodes = chebyshev_nodes(n);
            assert_eq!(nodes.len(), n);
            for &x in &nodes {
                assert!(x >= -1.0 && x <= 1.0);
            }
        }
    }

    #[test]
    fn chebyshev_coeffs_and_eval() {
        // Approximate f(x) = x^2 by a Chebyshev series on [-1, 1].
        // T_2(x) = 2x² − 1 → x² = (T_2(x) + 1)/2 = 0.5·T_0 + 0.5·T_2.
        let coeffs = chebyshev_coefficients(|x| x * x, 4);
        assert!(close(coeffs[0], 0.5, 1e-10), "c0={}", coeffs[0]);
        assert!(close(coeffs[2], 0.5, 1e-10), "c2={}", coeffs[2]);
        assert!(close(coeffs[1], 0.0, 1e-10), "c1={}", coeffs[1]);
        assert!(close(coeffs[3], 0.0, 1e-10), "c3={}", coeffs[3]);
        // Re-evaluate and check accuracy.
        for &x in &[-0.9, -0.5, 0.0, 0.3, 0.7, 1.0] {
            assert!(close(chebyshev_eval(&coeffs, x), x * x, 1e-10),
                    "at x={}: got {}", x, chebyshev_eval(&coeffs, x));
        }
    }

    #[test]
    fn chebyshev_approximates_smooth_function() {
        // Approximate exp(x) on [-1, 1] — a series of length 8 should be plenty.
        let coeffs = chebyshev_coefficients(|x| x.exp(), 8);
        let max_err = [-0.99, -0.7, -0.4, 0.0, 0.3, 0.6, 0.99]
            .iter()
            .map(|&x| (chebyshev_eval(&coeffs, x) - x.exp()).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_err < 1e-4, "max_err={}", max_err);
    }

    #[test]
    fn chebyshev_rescale_basic() {
        // On [-1, 1], rescale is identity.
        assert!(close(chebyshev_rescale(0.0, -1.0, 1.0), 0.0, 1e-15));
        // On [0, 1]: rescale(0.5) = (1 - 0 - 1) = 0.
        assert!(close(chebyshev_rescale(0.5, 0.0, 1.0), 0.0, 1e-15));
        // On [0, 1]: rescale(1.0) = (2 - 0 - 1) = 1.
        assert!(close(chebyshev_rescale(1.0, 0.0, 1.0), 1.0, 1e-15));
        // On [0, 10]: rescale(5) = (10 - 0 - 10) / 10 = 0.
        assert!(close(chebyshev_rescale(5.0, 0.0, 10.0), 0.0, 1e-15));
    }

    #[test]
    fn legendre_recurrence() {
        // P_0 = 1, P_1 = x, P_2 = (3x² − 1)/2, P_3 = (5x³ − 3x)/2.
        assert!(close(legendre_p(0, 0.0), 1.0, 1e-15));
        assert!(close(legendre_p(1, 0.5), 0.5, 1e-15));
        assert!(close(legendre_p(2, 0.5), (3.0 * 0.25 - 1.0) / 2.0, 1e-15));
        assert!(close(legendre_p(3, 0.5), (5.0 * 0.125 - 3.0 * 0.5) / 2.0, 1e-15));
        // Endpoint values
        assert!(close(legendre_p(5, 1.0), 1.0, 1e-12));
        assert!(close(legendre_p(5, -1.0), -1.0_f64.powi(5), 1e-12));
    }

    #[test]
    fn legendre_orthogonality() {
        // ∫_{-1}^1 P_m(x) · P_n(x) dx = 2/(2n+1) if m == n else 0.
        use crate::calculus::integrate_romberg;
        for (m, n) in [(0, 0), (1, 1), (2, 2), (3, 3), (0, 2), (1, 3)] {
            let fm_n = move |x: f64| legendre_p(m, x) * legendre_p(n, x);
            let integral = integrate_romberg(fm_n, -1.0, 1.0, 10).unwrap();
            let expected = if m == n { 2.0 / (2.0 * n as f64 + 1.0) } else { 0.0 };
            let err = (integral - expected).abs();
            assert!(err < 1e-9, "m={} n={}: integral={} expected={} err={}", m, n, integral, expected, err);
        }
    }

    #[test]
    fn legendre_associated_basic() {
        // P_1^1 = -sqrt(1 - x²)
        assert!(close(legendre_associated(1, 1, 0.5), -(1.0 - 0.25_f64).sqrt(), 1e-10));
        // P_2^2 = 3 (1 - x²)
        assert!(close(legendre_associated(2, 2, 0.5), 3.0 * (1.0 - 0.25), 1e-10));
        // P_2^0 = P_2 = (3x² − 1)/2
        assert!(close(legendre_associated(2, 0, 0.5), (3.0 * 0.25 - 1.0) / 2.0, 1e-10));
        // m > n ⇒ 0
        assert!(close(legendre_associated(1, 2, 0.5), 0.0, 1e-15));
    }

    #[test]
    fn gauss_legendre_exact_for_polynomial() {
        // 2-point Gauss–Legendre exactly integrates polynomials up to degree 3.
        let (xs, ws) = gauss_legendre(2);
        // ∫_{-1}^1 x^3 dx = 0
        let integral: f64 = ws.iter().zip(xs.iter()).map(|(w, x)| w * x.powi(3)).sum();
        assert!(integral.abs() < 1e-12);
        // ∫_{-1}^1 1 dx = 2
        let integral: f64 = ws.iter().sum();
        assert!(close(integral, 2.0, 1e-10));
    }

    #[test]
    fn gauss_legendre_weights_sum() {
        // ∫_{-1}^1 1 dx = 2; weights should sum to 2.
        for n in [1, 2, 3, 4, 5, 8] {
            let (_, ws) = gauss_legendre(n);
            let s: f64 = ws.iter().sum();
            assert!(close(s, 2.0, 1e-9), "n={}: sum={}", n, s);
        }
    }

    #[test]
    fn gauss_legendre_nodes_in_range() {
        for n in [1, 2, 3, 5, 10] {
            let (xs, _) = gauss_legendre(n);
            assert_eq!(xs.len(), n);
            for &x in &xs {
                assert!(x >= -1.0 && x <= 1.0, "n={} x={}", n, x);
            }
        }
    }
}
