//! Equation solvers: Newton–Raphson, bisection, secant, Brent, and polynomial
//! root finding via the Durand–Kerner method.

use crate::error::{MathError, Result};

/// Configuration shared by all iterative solvers.
#[derive(Debug, Clone)]
pub struct SolveOptions {
    pub max_iter: usize,
    pub tol: f64,
    pub h: f64, // finite-difference step
}

impl Default for SolveOptions {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-10,
            h: 1e-6,
        }
    }
}

/// Bisection method: finds a root of `f` on `[a, b]` assuming `f(a)*f(b) < 0`.
pub fn bisect<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut lo = a;
    let mut hi = b;
    let fa = f(lo);
    let fb = f(hi);
    if fa == 0.0 {
        return Ok((lo, 0.0));
    }
    if fb == 0.0 {
        return Ok((hi, 0.0));
    }
    if fa.signum() == fb.signum() {
        return Err(MathError::InvalidArgument(format!(
            "bisect: f(a) and f(b) must have opposite signs (got f({})={}, f({})={})",
            lo, fa, hi, fb
        )));
    }
    let mut mid;
    for _ in 0..opts.max_iter {
        mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if fm.abs() < opts.tol || (hi - lo) * 0.5 < opts.tol {
            return Ok((mid, fm));
        }
        if fm.signum() == fa.signum() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Err(MathError::NotConvergent(format!(
        "bisect did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Newton–Raphson method. Requires a derivative `df`; if you don't have one,
/// use [`newton_central`] which approximates it with a centred difference.
pub fn newton<F: Fn(f64) -> f64, D: Fn(f64) -> f64>(
    f: F,
    df: D,
    x0: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut x = x0;
    let mut fx = f(x);
    for i in 0..opts.max_iter {
        if fx.abs() < opts.tol {
            return Ok((x, fx));
        }
        let dfx = df(x);
        if dfx == 0.0 {
            return Err(MathError::NotConvergent(format!(
                "newton: derivative vanished at iteration {}",
                i
            )));
        }
        let step = fx / dfx;
        x -= step;
        fx = f(x);
        if step.abs() < opts.tol {
            return Ok((x, fx));
        }
    }
    Err(MathError::NotConvergent(format!(
        "newton did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Newton method with a numerically-estimated derivative (central difference).
pub fn newton_central<F>(f: F, x0: f64, opts: SolveOptions) -> Result<(f64, f64)>
where
    F: Fn(f64) -> f64 + Clone,
{
    let h = opts.h;
    let f2 = f.clone();
    let df = move |x: f64| (f2(x + h) - f2(x - h)) / (2.0 * h);
    newton(f, df, x0, opts)
}

/// Secant method: derivative-free, but starts with two initial guesses.
pub fn secant<F: Fn(f64) -> f64>(
    f: F,
    x0: f64,
    x1: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    let mut x_prev = x0;
    let mut x_curr = x1;
    let mut f_prev = f(x_prev);
    let mut f_curr = f(x_curr);
    for _ in 0..opts.max_iter {
        if f_curr.abs() < opts.tol {
            return Ok((x_curr, f_curr));
        }
        let denom = f_curr - f_prev;
        if denom.abs() < 1e-30 {
            return Err(MathError::NotConvergent(
                "secant: zero denominator".into(),
            ));
        }
        let x_next = x_curr - f_curr * (x_curr - x_prev) / denom;
        if (x_next - x_curr).abs() < opts.tol {
            let fn_ = f(x_next);
            return Ok((x_next, fn_));
        }
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_next;
        f_curr = f(x_next);
    }
    Err(MathError::NotConvergent(format!(
        "secant did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Brent's method for a bracketed root.
///
/// Combines bisection with secant and inverse-quadratic interpolation: the
/// sign-changing bracket is always maintained, so the method cannot diverge
/// the way a pure secant iteration can, while it usually converges much
/// faster than plain bisection. The bracket endpoints may be given in either
/// order.
///
/// The state machine follows the classic Zeroin/Brent formulation: `xcur` is
/// the current best point, `xblk` the last point of opposite sign, and `xpre`
/// the previous iterate (used for the interpolation steps).
pub fn brent<F: Fn(f64) -> f64>(
    f: F,
    a: f64,
    b: f64,
    opts: SolveOptions,
) -> Result<(f64, f64)> {
    if !a.is_finite() || !b.is_finite() {
        return Err(MathError::InvalidArgument(
            "brent: bracket endpoints must be finite".into(),
        ));
    }
    let mut xpre = a;
    let mut xcur = b;
    let mut fpre = f(xpre);
    let mut fcur = f(xcur);
    if !fpre.is_finite() || !fcur.is_finite() {
        return Err(MathError::Domain(
            "brent: function must be finite at both bracket endpoints".into(),
        ));
    }
    if fpre == 0.0 {
        return Ok((xpre, fpre));
    }
    if fcur == 0.0 {
        return Ok((xcur, fcur));
    }
    if fpre.signum() == fcur.signum() {
        return Err(MathError::InvalidArgument(format!(
            "brent: f(a) and f(b) must have opposite signs (got f({})={}, f({})={})",
            xpre, fpre, xcur, fcur
        )));
    }

    // The blocking point: the last iterate whose function value has the
    // opposite sign to the current one.
    let mut xblk = xpre;
    let mut fblk = fpre;
    // The last step actually taken, and the previous one (for the
    // interpolation fallback).
    let mut spre = 0.0;
    let mut scur = 0.0;

    for _ in 0..opts.max_iter {
        // Restore the bracket whenever the current and previous values
        // straddle the root. Both step trackers are reset to the step just
        // taken: that keeps the interpolation acceptance test strict, so the
        // method falls back to bisection (which shrinks the bracket) instead
        // of stalling on a point it never crosses.
        if fpre != 0.0 && fcur != 0.0 && fpre.signum() != fcur.signum() {
            xblk = xpre;
            fblk = fpre;
            spre = xcur - xpre;
            scur = spre;
        }
        // Keep the point with the smaller |f| as the current one.
        if fblk.abs() < fcur.abs() {
            let (nxpre, nxcur, nxblk) = (xcur, xblk, xcur);
            let (nfpre, nfcur, nfblk) = (fcur, fblk, fcur);
            xpre = nxpre;
            xcur = nxcur;
            xblk = nxblk;
            fpre = nfpre;
            fcur = nfcur;
            fblk = nfblk;
        }

        let delta = (opts.tol + 4.0 * f64::EPSILON * xcur.abs()) / 2.0;
        let sbis = (xblk - xcur) / 2.0;
        if fcur == 0.0 || sbis.abs() < delta {
            return Ok((xcur, fcur));
        }

        // Try an interpolated step; fall back to bisection when it would
        // leave the bracket or is too close to the previous point.
        if spre.abs() > delta && fcur.abs() < fpre.abs() {
            let stry = if xpre == xblk {
                -fcur * (xcur - xpre) / (fcur - fpre)
            } else {
                let dpre = (fpre - fcur) / (xpre - xcur);
                let dblk = (fblk - fcur) / (xblk - xcur);
                -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre))
            };
            if 2.0 * stry.abs() < spre.abs().min(3.0 * sbis.abs() - delta) {
                spre = scur;
                scur = stry;
            } else {
                spre = sbis;
                scur = sbis;
            }
        } else {
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur;
        fpre = fcur;
        if scur.abs() > delta {
            xcur += scur;
        } else {
            // Never move less than the convergence threshold, so the
            // iteration cannot stall before it reaches it.
            xcur += if sbis > 0.0 { delta } else { -delta };
        }
        fcur = f(xcur);
        if !fcur.is_finite() {
            return Err(MathError::Domain(format!(
                "brent: function is not finite at x = {}",
                xcur
            )));
        }
    }
    Err(MathError::NotConvergent(format!(
        "brent did not converge in {} iterations",
        opts.max_iter
    )))
}

/// Configuration for the sign-scan root finder.
#[derive(Debug, Clone)]
pub struct ScanOptions {
    /// Left end of the scan interval.
    pub lo: f64,
    /// Right end of the scan interval.
    pub hi: f64,
    /// Number of sample points across the interval.
    pub samples: usize,
    /// Tolerance for accepting a root.
    pub tol: f64,
    /// Minimum separation between two distinct roots.
    pub dedup: f64,
}

impl Default for ScanOptions {
    fn default() -> Self {
        Self {
            lo: -30.0,
            hi: 30.0,
            samples: 2400,
            tol: 1e-9,
            dedup: 1e-7,
        }
    }
}

/// Find **all** real roots of `f` on an interval by scanning for sign changes.
///
/// A single Newton–Raphson run finds one root near the initial guess. This
/// samples the interval, bisects every sign change, and collects the
/// results — so `sin(x) = 0` yields `0, ±π, ±2π, …` rather than just the
/// root nearest the guess.
///
/// **Pole rejection**: a sign flip across a *pole* looks exactly like a sign
/// flip across a root, so a bisected bracket is only accepted when the
/// function actually lands near zero there. Without this, `1/(x−1) = 0`
/// would "solve" at the pole `x = 1`.
///
/// Points where the function is undefined (NaN) break the scan rather than
/// producing a spurious bracket.
///
/// Cross-port pattern from unsolve `solver.ts:286-329` (`numericSolve`).
pub fn find_all_roots<F: Fn(f64) -> f64>(f: F, opts: &ScanOptions) -> Vec<f64> {
    let n = opts.samples.max(2);
    let lo = opts.lo;
    let hi = opts.hi;
    let step = (hi - lo) / (n as f64 - 1.0);

    let mut roots: Vec<f64> = Vec::new();
    let mut prev_x = lo;
    let mut prev_y = f(lo);

    for i in 1..=n {
        let x = lo + i as f64 * step;
        let y = f(x);
        if prev_y.is_finite() && y.is_finite() {
            if prev_y == 0.0 {
                // An exact zero at a sample point is a root.
                roots.push(prev_x);
            } else if prev_y * y < 0.0 {
                // A sign change — bisect, but only accept a real root.
                if let Some(root) = bisect_bracket(&f, prev_x, x, prev_y, opts.tol) {
                    roots.push(root);
                }
            }
        }
        prev_x = x;
        prev_y = y;
    }

    dedupe_roots(roots, opts.dedup)
}

/// Bisect a sign-change bracket, rejecting poles.
///
/// A pole flips the sign just like a root does, so the result is only
/// accepted when `f` actually lands near zero at the bisection point —
/// otherwise `1/(x−p) = 0` would "solve" at the pole.
fn bisect_bracket<F: Fn(f64) -> f64>(
    f: &F,
    mut lo: f64,
    mut hi: f64,
    mut flo: f64,
    tol: f64,
) -> Option<f64> {
    // The scale of the endpoint values bounds how small |f| must be.
    let scale = flo.abs().min(f(hi).abs());
    for _ in 0..80 {
        if hi - lo <= 1e-14 {
            break;
        }
        let mid = 0.5 * (lo + hi);
        let fm = f(mid);
        if !fm.is_finite() {
            // Undefined inside the bracket — a pole, not a root.
            return None;
        }
        if flo * fm <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            flo = fm;
        }
    }
    let x = 0.5 * (lo + hi);
    let fx = f(x);
    if fx.is_finite() && fx.abs() <= tol.max(scale) {
        Some(x)
    } else {
        None
    }
}

/// Remove near-duplicate roots, keeping the first of each cluster.
fn dedupe_roots(roots: Vec<f64>, dedup: f64) -> Vec<f64> {
    let mut sorted = roots;
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<f64> = Vec::new();
    for r in sorted {
        if let Some(&last) = out.last() {
            if (r - last).abs() <= dedup {
                continue;
            }
        }
        out.push(r);
    }
    out
}

/// Find all (real) roots of a polynomial with real coefficients given by
/// `coeffs` in descending order: `coeffs[0]*x^n + ... + coeffs[n]`.
///
/// Implemented with the Durand–Kerner method, which converges quickly for
/// well-separated roots and is robust enough for moderate degrees.
pub fn polynomial_roots(coeffs: &[f64]) -> Result<Vec<(f64, f64)>> {
    let coeffs = trim_near_zero_leading(coeffs);
    let roots = durand_kerner_all_roots(&coeffs)?;

    // A repeated root is one real solution. Keep multiplicities in the complex
    // API (apart needs them), but coalesce them in this user-facing list.
    let xs = dedupe_roots(
        roots.into_iter().filter(|z| z.im == 0.0).map(|z| z.re).collect(),
        1e-6,
    );
    Ok(xs
        .into_iter()
        .map(|x| {
            let fx = polynomial_eval(&coeffs, x);
            (x, fx)
        })
        .collect())
}

/// Find all complex roots of a polynomial with real coefficients, returned
/// as `(re, im)` pairs via the Durand–Kerner method. Complex roots appear
/// together with their conjugates. Repeated roots are retained so callers can
/// recover their multiplicities.
pub fn polynomial_roots_complex(coeffs: &[f64]) -> Result<Vec<(f64, f64)>> {
    let coeffs = trim_near_zero_leading(coeffs);
    let roots = durand_kerner_all_roots(&coeffs)?;
    let mut out: Vec<(f64, f64)> = roots.iter().map(|z| (z.re, z.im)).collect();
    out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    Ok(out)
}

/// Discard leading coefficients that are negligible relative to the largest
/// coefficient. Exact leading zeros and tiny high-degree terms otherwise make
/// Durand–Kerner report spurious roots at enormous magnitudes.
fn trim_near_zero_leading(coeffs: &[f64]) -> Vec<f64> {
    if coeffs.is_empty() {
        return Vec::new();
    }
    let scale = coeffs.iter().fold(0.0_f64, |m, c| m.max(c.abs()));
    let tol = 1e-14 * scale;
    let first = coeffs
        .iter()
        .position(|c| c.abs() > tol)
        .unwrap_or(coeffs.len() - 1);
    coeffs[first..].to_vec()
}

fn snap_integer(x: f64, tol: f64) -> f64 {
    let nearest = x.round();
    if (x - nearest).abs() < tol {
        nearest
    } else {
        x
    }
}

fn polish_root(mut z: crate::complex::Complex<f64>) -> crate::complex::Complex<f64> {
    if z.im.abs() < 1e-8 * z.abs().max(1.0) {
        crate::complex::Complex::new(snap_integer(z.re, 1e-7), 0.0)
    } else {
        z.re = snap_integer(z.re, 1e-9);
        z
    }
}

/// Run the Durand–Kerner iteration and return all n complex roots of the
/// (real-coefficient, descending) polynomial.
fn durand_kerner_all_roots(coeffs: &[f64]) -> Result<Vec<crate::complex::Complex<f64>>> {
    let n = coeffs.len().saturating_sub(1);
    if n == 0 {
        return Ok(Vec::new());
    }
    // normalize so leading coefficient is 1
    let lead = coeffs[0];
    if lead == 0.0 {
        return Err(MathError::InvalidArgument("polynomial leading coefficient is 0".into()));
    }
    let a: Vec<f64> = coeffs.iter().map(|c| c / lead).collect();

    // initial guesses: complex roots spread around the unit circle
    use crate::complex::Complex;
    let mut z: Vec<Complex<f64>> = (0..n)
        .map(|k| {
            let theta = std::f64::consts::PI * (2.0 * k as f64 + 1.0) / n as f64;
            Complex::from_polar(0.4 + 0.9 * (n as f64 / 10.0), theta)
        })
        .collect();

    for _ in 0..200 {
        let mut converged = true;
        for i in 0..n {
            let mut p = Complex::new(a[0], 0.0);
            for j in 1..a.len() {
                p = p * z[i] + Complex::new(a[j], 0.0);
            }
            let mut denom = Complex::new(1.0, 0.0);
            for j in 0..n {
                if i != j {
                    denom = denom * (z[i] - z[j]);
                }
            }
            let step = if denom.abs() < 1e-30 {
                Complex::new(0.0, 0.0)
            } else {
                p / denom
            };
            z[i] = z[i] - step;
            if step.abs() > 1e-10 {
                converged = false;
            }
        }
        if converged {
            break;
        }
    }

    Ok(z.into_iter().map(polish_root).collect())
}

fn polynomial_eval(coeffs: &[f64], x: f64) -> f64 {
    let mut acc = 0.0;
    for c in coeffs {
        acc = acc * x + c;
    }
    acc
}

// ---------------------------------------------------------------------------
// Polynomial root isolation via the Vincent–Akritas–Strzebonski (VAS) method.
//
// Given a polynomial with integer coefficients, finds disjoint open intervals
// each containing exactly one real root.  Uses Descartes' rule of signs on
// Möbius-transformed polynomials with exact i128 arithmetic.
// ---------------------------------------------------------------------------

/// Count the number of sign variations (sign changes) in a coefficient
/// sequence, ignoring zeros.
fn count_sign_variations(coeffs: &[i128]) -> usize {
    let mut count = 0;
    let mut prev_sign: i8 = 0;
    for &c in coeffs {
        let sign = if c > 0 { 1 } else if c < 0 { -1 } else { 0 };
        if sign != 0 {
            if prev_sign != 0 && sign != prev_sign {
                count += 1;
            }
            prev_sign = sign;
        }
    }
    count
}

/// Compute p(x + c) for a polynomial with integer coefficients
/// (highest-degree-first ordering).
/// Uses Horner's method: process from the leading coefficient.
fn poly_translate(coeffs: &[i128], c: i128) -> Vec<i128> {
    let n = coeffs.len();
    if n == 0 {
        return Vec::new();
    }
    // Horner: r = a_n; r = r*(x+c) + a_{n-1}; ... ; r = r*(x+c) + a_0
    // r * (x + c) in highest-degree-first:
    //   r*x appends 0 at the end (shift up one degree)
    //   c*r prepends 0 at the front (same degree, padded)
    //   sum: new[0] = old[0], new[j] = old[j] + c*old[j-1], new[len] = c*old[len-1]
    let mut result = vec![coeffs[0]];
    for i in 1..n {
        let len = result.len();
        let mut new_result = vec![0i128; len + 1];
        new_result[0] = result[0];
        for j in 1..len {
            new_result[j] = result[j] + c.saturating_mul(result[j - 1]);
        }
        new_result[len] = c.saturating_mul(result[len - 1]);
        // Add coeffs[i] to the constant term (last element in highest-degree-first)
        new_result[len] += coeffs[i];
        result = new_result;
    }
    result
}

/// Compute x^n * p(1/x) — reverses the coefficient list.
/// If p(x) = a_0 + a_1*x + ... + a_n*x^n, then x^n * p(1/x) = a_n + a_{n-1}*x + ... + a_0*x^n.
fn poly_reciprocal(coeffs: &[i128]) -> Vec<i128> {
    coeffs.iter().rev().copied().collect()
}

/// Remove leading zeros (high-degree zero coefficients).
fn trim_leading_zeros(coeffs: &[i128]) -> Vec<i128> {
    let mut start = 0;
    while start < coeffs.len() - 1 && coeffs[start] == 0 {
        start += 1;
    }
    coeffs[start..].to_vec()
}

/// Compute a Cauchy bound: an integer upper bound on the absolute value of
/// all real roots of the polynomial.
fn cauchy_bound(coeffs: &[i128]) -> i128 {
    let n = coeffs.len();
    if n <= 1 {
        return 0;
    }
    let an = coeffs[0].abs();
    if an == 0 {
        return 1;
    }
    let mut max_ratio: f64 = 0.0;
    for i in 1..n {
        let ratio = coeffs[i].abs() as f64 / an as f64;
        if ratio > max_ratio {
            max_ratio = ratio;
        }
    }
    let bound = max_ratio + 1.0;
    bound.ceil() as i128
}

/// Recursively isolate positive roots of a polynomial.
///
/// `coeffs` are the polynomial coefficients (highest degree first).
/// The Möbius transformation (a, b, c, d) maps the current variable t
/// to the original x via x = (a*t + b) / (c*t + d).
/// When V(coeffs) = 1, the root is in the interval (b/d, a/c).
/// `check_zero` controls whether to detect roots at t=0 (enabled for the
/// right/shift branch to avoid duplicate detection at split points).
fn isolate_positive_roots(
    coeffs: &[i128],
    a: i128,
    b: i128,
    c: i128,
    d: i128,
    depth: usize,
    check_zero: bool,
) -> Vec<(f64, f64)> {
    if depth > 200 {
        return Vec::new(); // safety limit
    }

    let trimmed = trim_leading_zeros(coeffs);
    let n = trimmed.len();

    // Check if 0 is a root (constant term is zero).
    if check_zero && n > 1 && trimmed[n - 1] == 0 {
        let root = b as f64 / d as f64;
        let mut results = vec![(root, root)];

        // Remove the factor t: divide by t (drop the trailing zero).
        let reduced = &trimmed[..n - 1];
        let v_reduced = count_sign_variations(reduced);
        if v_reduced > 0 {
            results.extend(isolate_positive_roots(reduced, a, b, c, d, depth + 1, true));
        } else if reduced.len() > 1 && reduced[reduced.len() - 1] == 0 {
            // Another root at 0 (repeated root).
            results.extend(isolate_positive_roots(reduced, a, b, c, d, depth + 1, true));
        }
        return results;
    }

    let v = count_sign_variations(&trimmed);

    if v == 0 {
        return Vec::new();
    }

    if v == 1 {
        // Exactly one positive root in (0, ∞) for the transformed variable t.
        // This maps to x in (b/d, a/c) via the Möbius transformation.
        let lo = b as f64 / d as f64;
        let hi = if c == 0 { f64::INFINITY } else { a as f64 / c as f64 };
        let (lo, hi) = if lo < hi { (lo, hi) } else { (hi, lo) };
        return vec![(lo, hi)];
    }

    let mut results = Vec::new();

    // Split: check (0, 1) via reciprocal, and (1, ∞) via shift.
    // For (1, ∞): transform p(x) → p(x + 1), new Möbius: (a, a+b, c, c+d)
    // check_zero=true because a root at t=0 here maps to the split point
    // and is a genuine root in the right interval.
    let shifted = poly_translate(&trimmed, 1);
    results.extend(isolate_positive_roots(
        &shifted,
        a,
        a + b,
        c,
        c + d,
        depth + 1,
        true,
    ));

    // For (0, 1): map (0, 1) to (0, ∞) via x → 1/(t+1).
    // q(t) = reciprocal(p)(t+1) — first reverse coefficients, then shift by 1.
    // New Möbius: (b, a+b, d, c+d)
    // check_zero=false because a root at t=0 here maps to the split point
    // which is already handled by the right branch.
    let recip = poly_reciprocal(&trimmed);
    let shifted_recip = poly_translate(&recip, 1);
    results.extend(isolate_positive_roots(
        &shifted_recip,
        b,
        a + b,
        d,
        c + d,
        depth + 1,
        false,
    ));

    results
}

/// Isolate all real roots of a polynomial with integer coefficients.
///
/// Returns a list of disjoint open intervals `(lo, hi)`, each containing
/// exactly one real root.  Uses the Vincent–Akritas–Strzebonski method
/// with exact i128 arithmetic.
///
/// Coefficients are given highest-degree-first: `p(x) = coeffs[0]*x^n + ... + coeffs[n]`.
pub fn isolate_real_roots(coeffs: &[i64]) -> Result<Vec<(f64, f64)>> {
    if coeffs.is_empty() {
        return Err(MathError::InvalidArgument("isolate_real_roots: empty coefficients".into()));
    }
    let n = coeffs.len() - 1;
    if n == 0 {
        return Ok(Vec::new());
    }

    // Convert to i128.
    let p: Vec<i128> = coeffs.iter().map(|&x| x as i128).collect();
    let p = trim_leading_zeros(&p);

    if p.len() <= 1 {
        return Ok(Vec::new());
    }

    // Compute Cauchy bound (used for informational purposes; the VAS
    // algorithm converges without it, but it could be used for an initial
    // shift to speed up convergence).
    let _bound = cauchy_bound(&p);

    let mut all_roots = Vec::new();

    // Positive roots (including 0): search (0, ∞) with identity Möbius (1, 0, 0, 1).
    // The root-at-0 case is handled inside isolate_positive_roots.
    let pos_roots = isolate_positive_roots(&p, 1, 0, 0, 1, 0, true);
    all_roots.extend(pos_roots);

    // Negative roots: use p(-x) with identity Möbius (1, 0, 0, 1).
    // Root at 0 is already found above, so skip it here.
    let mut neg_coeffs: Vec<i128> = p.iter().enumerate().map(|(i, &c)| {
        if (p.len() - 1 - i).is_multiple_of(2) { c } else { -c }
    }).collect();
    neg_coeffs = trim_leading_zeros(&neg_coeffs);
    if neg_coeffs.len() > 1 {
        let neg_roots = isolate_positive_roots(&neg_coeffs, 1, 0, 0, 1, 0, true);
        // Negate the intervals to get negative roots, skipping root at 0.
        for (lo, hi) in neg_roots {
            if lo > 0.0 || hi > 0.0 {
                all_roots.push((-hi, -lo));
            }
        }
    }

    // Sort and return.
    all_roots.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    Ok(all_roots)
}

/// Solve a system of `n` nonlinear equations in `n` unknowns using
/// Newton's method with a finite-difference Jacobian.
///
/// `system` is a function that maps the current guess `x ∈ ℝⁿ` to
/// the residual vector `f(x) ∈ ℝⁿ`. The algorithm iterates
/// `x ← x − J(x)⁻¹ f(x)` until either `‖f(x)‖ < tol` or `max_iter` is hit.
///
/// On non-convergence or singular Jacobian, returns [`MathError::NotConvergent`].
pub fn newton_system<F: Fn(&[f64]) -> Vec<f64>>(
    system: F,
    x0: &[f64],
    opts: SolveOptions,
) -> Result<Vec<f64>> {
    let n = x0.len();
    if n == 0 {
        return Err(MathError::InvalidArgument("newton_system: empty initial guess".into()));
    }
    let h = if opts.h <= 0.0 { 1e-6 } else { opts.h };
    let mut x = x0.to_vec();
    for _ in 0..opts.max_iter {
        let f0 = system(&x);
        if f0.len() != n {
            return Err(MathError::InvalidArgument(format!(
                "newton_system: system returned {} values, expected {}",
                f0.len(), n
            )));
        }
        let norm: f64 = f0.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < opts.tol {
            return Ok(x);
        }
        // Build Jacobian by central finite differences.
        let mut jac = vec![0.0_f64; n * n];
        for j in 0..n {
            let mut xp = x.clone();
            let mut xm = x.clone();
            xp[j] += h;
            xm[j] -= h;
            let fp = system(&xp);
            let fm = system(&xm);
            for i in 0..n {
                jac[i * n + j] = (fp[i] - fm[i]) / (2.0 * h);
            }
        }
        // Solve J · delta = -f0 for delta.
        let neg_f0: Vec<f64> = f0.iter().map(|v| -v).collect();
        let m = crate::matrix::Matrix::from_row_major(n, n, jac.clone())
            .map_err(|e| MathError::InvalidArgument(format!("newton_system: {}", e)))?;
        let delta = m
            .solve(&neg_f0)
            .map_err(|_| MathError::NotConvergent("newton_system: singular Jacobian".into()))?;
        for i in 0..n {
            x[i] += delta[i];
        }
    }
    Err(MathError::NotConvergent(format!(
        "newton_system: failed to converge in {} iterations",
        opts.max_iter
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn bisect_root() {
        // x^2 - 4 = 0 on [0, 5]
        let f = |x: f64| x * x - 4.0;
        let (r, _) = bisect(f, 0.0, 5.0, SolveOptions::default()).unwrap();
        assert!((r - 2.0).abs() < 1e-6);
    }

    #[test]
    fn newton_root() {
        // cos(x) = x  -> Newton's with f=cos(x)-x, df=-sin(x)-1
        let f = |x: f64| x.cos() - x;
        let df = |x: f64| -x.sin() - 1.0;
        let (r, _) = newton(f, df, 0.5, SolveOptions::default()).unwrap();
        // known fixed point ≈ 0.7390851332151607
        assert!((r - 0.7390851332151607).abs() < 1e-8);
    }

    #[test]
    fn newton_central_root() {
        // Newton from 3 on sin(x) converges to the closest root, which is π.
        let f = |x: f64| x.sin();
        let (r, _) = newton_central(f, 3.0, SolveOptions::default()).unwrap();
        assert!((r - PI).abs() < 1e-7, "expected root at π, got {}", r);
    }

    #[test]
    fn secant_root() {
        let f = |x: f64| x.powi(3) - 2.0;
        let (r, _) = secant(f, 1.0, 2.0, SolveOptions::default()).unwrap();
        assert!((r - 2f64.cbrt()).abs() < 1e-8);
    }

    #[test]
    fn brent_root() {
        let f = |x: f64| x * x - 4.0;
        let (r, fx) = brent(f, 0.0, 5.0, SolveOptions::default()).unwrap();
        assert!((r - 2.0).abs() < 1e-10, "root = {}", r);
        assert!(fx.abs() < 1e-9, "f(root) = {}", fx);
    }

    #[test]
    fn brent_accepts_reversed_bracket_and_exact_endpoints() {
        let f = |x: f64| x * x - 4.0;
        let (r, _) = brent(f, 5.0, 0.0, SolveOptions::default()).unwrap();
        assert!((r - 2.0).abs() < 1e-10);
        let (r, fx) = brent(f, 2.0, 5.0, SolveOptions::default()).unwrap();
        assert_eq!((r, fx), (2.0, 0.0));
    }

    #[test]
    fn brent_requires_a_finite_sign_change() {
        let f = |x: f64| x * x + 1.0;
        let err = brent(f, 0.0, 2.0, SolveOptions::default()).unwrap_err();
        assert!(err.to_string().contains("opposite signs"), "got {}", err);

        let f = |x: f64| if x < 1.0 { f64::NAN } else { x - 2.0 };
        let err = brent(f, 0.0, 3.0, SolveOptions::default()).unwrap_err();
        assert!(err.to_string().contains("finite"), "got {}", err);
    }

    #[test]
    fn brent_handles_a_flat_root() {
        // x^3 has a triple root at 0 — a sign change, but a flat one, so the
        // method converges linearly (rate 2/3) like Newton instead of
        // superlinearly. The iteration budget has to allow for that.
        let f = |x: f64| x.powi(3);
        let opts = SolveOptions { max_iter: 1000, ..SolveOptions::default() };
        let (r, fx) = brent(f, -1.0, 2.0, opts).unwrap();
        assert!(r.abs() < 1e-8, "root = {}", r);
        assert!(fx.abs() < 1e-20, "f(root) = {}", fx);
    }

    #[test]
    fn brent_uses_fewer_evaluations_than_bisection() {
        use std::cell::Cell;
        let brent_evals = Cell::new(0usize);
        let bisect_evals = Cell::new(0usize);
        let bf = |x: f64| {
            brent_evals.set(brent_evals.get() + 1);
            x.cos() - x
        };
        let sf = |x: f64| {
            bisect_evals.set(bisect_evals.get() + 1);
            x.cos() - x
        };
        let opts = SolveOptions { max_iter: 100, tol: 1e-10, h: 1e-6 };
        let (br, _) = brent(bf, 0.0, 1.0, opts.clone()).unwrap();
        let (sr, _) = bisect(sf, 0.0, 1.0, opts).unwrap();
        assert!((br - sr).abs() < 1e-8);
        assert!(
            brent_evals.get() < bisect_evals.get(),
            "brent {} evals, bisection {} evals",
            brent_evals.get(),
            bisect_evals.get()
        );
    }

    #[test]
    fn polynomial_two_roots() {
        // x^2 - 5x + 6 = (x-2)(x-3)
        let r = polynomial_roots(&[1.0, -5.0, 6.0]).unwrap();
        let xs: Vec<f64> = r.iter().map(|(x, _)| *x).collect();
        assert_eq!(xs.len(), 2);
        let mut sorted = xs.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 2.0).abs() < 1e-4);
        assert!((sorted[1] - 3.0).abs() < 1e-4);
    }

    #[test]
    fn polynomial_trims_exact_and_near_zero_leading_coefficients() {
        for coeffs in [&[0.0, 1.0, -5.0, 6.0][..], &[1e-20, 1.0, -5.0, 6.0][..]] {
            let roots = polynomial_roots(coeffs).unwrap();
            let xs: Vec<f64> = roots.iter().map(|(x, _)| *x).collect();
            assert_eq!(xs, vec![2.0, 3.0], "coeffs = {:?}", coeffs);
        }
    }

    #[test]
    fn polynomial_leading_trimming_is_scale_relative() {
        // Uniformly tiny coefficients describe the same cubic as [1, 1, -5, 6];
        // the first coefficient must not be discarded merely for being small.
        let tiny = [1e-20, 1e-20, -5e-20, 6e-20];
        assert_eq!(trim_near_zero_leading(&tiny), tiny);
        assert_eq!(trim_near_zero_leading(&[0.0, 0.0]), vec![0.0]);
    }

    #[test]
    fn polynomial_repeated_real_root_is_listed_once() {
        let roots = polynomial_roots(&[1.0, -2.0, 1.0]).unwrap();
        assert_eq!(roots, vec![(1.0, 0.0)]);
    }

    #[test]
    fn polynomial_complex_roots_are_snapped() {
        let mut roots = polynomial_roots_complex(&[1.0, 0.0, 1.0]).unwrap();
        roots.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        assert_eq!(roots, vec![(0.0, -1.0), (0.0, 1.0)]);
    }

    #[test]
    fn polynomial_complex_roots_keep_multiplicity() {
        // apart relies on repeated complex roots to reconstruct multiplicity.
        let roots = polynomial_roots_complex(&[1.0, -2.0, 1.0]).unwrap();
        assert_eq!(roots, vec![(1.0, 0.0), (1.0, 0.0)]);
    }

    #[test]
    fn newton_system_linear() {
        // 2x + y = 5; x + 3y = 7  →  (x, y) = (1.6, 1.8).
        let system = |x: &[f64]| vec![2.0 * x[0] + x[1] - 5.0, x[0] + 3.0 * x[1] - 7.0];
        let sol = newton_system(system, &[0.0, 0.0], SolveOptions::default()).unwrap();
        assert!((sol[0] - 1.6).abs() < 1e-6, "x = {}", sol[0]);
        assert!((sol[1] - 1.8).abs() < 1e-6, "y = {}", sol[1]);
    }

    #[test]
    fn newton_system_nonlinear() {
        // Intersection of x² + y² = 1 and y = x² − 0.5.
        let system = |x: &[f64]| {
            vec![
                x[0] * x[0] + x[1] * x[1] - 1.0,
                x[1] - x[0] * x[0] + 0.5,
            ]
        };
        let sol = newton_system(system, &[0.9, 0.3], SolveOptions::default()).unwrap();
        let r = (sol[0] * sol[0] + sol[1] * sol[1] - 1.0).abs();
        assert!(r < 1e-6, "residual = {}", r);
        let s = (sol[1] - sol[0] * sol[0] + 0.5).abs();
        assert!(s < 1e-6, "residual = {}", s);
    }

    #[test]
    fn newton_system_3d() {
        // 3×3 linear: x + y + z = 6, 2x - y + z = 3, x + 2y - z = 2.
        // Solution: x = 1, y = 2, z = 3.
        let system = |v: &[f64]| {
            vec![
                v[0] + v[1] + v[2] - 6.0,
                2.0 * v[0] - v[1] + v[2] - 3.0,
                v[0] + 2.0 * v[1] - v[2] - 2.0,
            ]
        };
        let sol = newton_system(system, &[0.0, 0.0, 0.0], SolveOptions::default()).unwrap();
        assert!((sol[0] - 1.0).abs() < 1e-6, "x = {}", sol[0]);
        assert!((sol[1] - 2.0).abs() < 1e-6, "y = {}", sol[1]);
        assert!((sol[2] - 3.0).abs() < 1e-6, "z = {}", sol[2]);
    }

    #[test]
    fn isolate_roots_quadratic() {
        // (x-2)(x-3) = x^2 - 5x + 6
        let intervals = isolate_real_roots(&[1, -5, 6]).unwrap();
        assert_eq!(intervals.len(), 2, "expected 2 roots, got {:?}", intervals);
        // Each interval should contain exactly one root.
        for (lo, hi) in &intervals {
            assert!(
                (lo <= &2.0 && hi >= &2.0) || (lo <= &3.0 && hi >= &3.0),
                "interval ({}, {}) should contain 2 or 3", lo, hi
            );
        }
    }

    #[test]
    fn isolate_roots_cubic() {
        // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
        let intervals = isolate_real_roots(&[1, -6, 11, -6]).unwrap();
        assert_eq!(intervals.len(), 3, "expected 3 roots, got {:?}", intervals);
        let roots = [1.0, 2.0, 3.0];
        for (i, r) in roots.iter().enumerate() {
            let (lo, hi) = intervals[i];
            assert!(lo <= *r && *r <= hi, "root {} not in interval ({}, {})", r, lo, hi);
        }
    }

    #[test]
    fn isolate_roots_with_negatives() {
        // (x+2)(x-1)(x-3) = x^3 - 2x^2 - 5x + 6
        let intervals = isolate_real_roots(&[1, -2, -5, 6]).unwrap();
        assert_eq!(intervals.len(), 3, "expected 3 roots, got {:?}", intervals);
        let roots = [-2.0, 1.0, 3.0];
        for (i, r) in roots.iter().enumerate() {
            let (lo, hi) = intervals[i];
            assert!(lo <= *r && *r <= hi, "root {} not in interval ({}, {})", r, lo, hi);
        }
    }

    #[test]
    fn isolate_roots_no_real_roots() {
        // x^2 + 1 has no real roots
        let intervals = isolate_real_roots(&[1, 0, 1]).unwrap();
        assert_eq!(intervals.len(), 0, "expected 0 roots, got {:?}", intervals);
    }

    #[test]
    fn isolate_roots_root_at_zero() {
        // x * (x - 1) = x^2 - x
        let intervals = isolate_real_roots(&[1, -1, 0]).unwrap();
        assert_eq!(intervals.len(), 2, "expected 2 roots, got {:?}", intervals);
        // One root should be at 0
        assert!(intervals[0].0 <= 0.0 && 0.0 <= intervals[0].1,
            "root 0 not in interval ({}, {})", intervals[0].0, intervals[0].1);
    }

    #[test]
    fn isolate_roots_repeated() {
        // (x-1)^2 = x^2 - 2x + 1 — repeated root at 1
        let intervals = isolate_real_roots(&[1, -2, 1]).unwrap();
        // Descartes' rule counts repeated roots as one interval
        assert!(!intervals.is_empty(), "expected at least 1 root, got {:?}", intervals);
        let (lo, hi) = intervals[0];
        assert!(lo <= 1.0 && 1.0 <= hi, "root 1 not in interval ({}, {})", lo, hi);
    }

    #[test]
    fn isolate_roots_higher_degree() {
        // (x+1)(x-1)(x-2)(x+3) = x^4 + x^3 - 7x^2 - x + 6
        let intervals = isolate_real_roots(&[1, 1, -7, -1, 6]).unwrap();
        assert_eq!(intervals.len(), 4, "expected 4 roots, got {:?}", intervals);
        let roots = [-3.0, -1.0, 1.0, 2.0];
        for (i, r) in roots.iter().enumerate() {
            let (lo, hi) = intervals[i];
            assert!(lo <= *r && *r <= hi, "root {} not in interval ({}, {})", r, lo, hi);
        }
    }

    // ---- sign-scan root finder --------------------------------------------

    #[test]
    fn scan_finds_all_sine_roots() {
        // sin(x) = 0 has roots at every multiple of pi.
        let opts = ScanOptions { lo: -7.0, hi: 7.0, ..Default::default() };
        let roots = find_all_roots(|x| x.sin(), &opts);
        assert_eq!(roots.len(), 5, "got {:?}", roots);
        let want = [-2.0 * std::f64::consts::PI, -std::f64::consts::PI, 0.0, std::f64::consts::PI, 2.0 * std::f64::consts::PI];
        for (r, w) in roots.iter().zip(want.iter()) {
            assert!((r - w).abs() < 1e-6, "got {} want {}", r, w);
        }
    }

    #[test]
    fn scan_finds_single_root() {
        // exp(x) = 5 has exactly one real root.
        let opts = ScanOptions::default();
        let roots = find_all_roots(|x| x.exp() - 5.0, &opts);
        assert_eq!(roots.len(), 1, "got {:?}", roots);
        assert!((roots[0] - 5.0f64.ln()).abs() < 1e-6, "got {}", roots[0]);
    }

    #[test]
    fn scan_finds_polynomial_roots() {
        let opts = ScanOptions::default();
        let roots = find_all_roots(|x| x * x - 4.0, &opts);
        assert_eq!(roots.len(), 2, "got {:?}", roots);
        assert!((roots[0] + 2.0).abs() < 1e-6);
        assert!((roots[1] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn scan_rejects_poles() {
        // 1/(x-1) has a pole at x=1 and no root — the sign flip there
        // must NOT be reported as a root.
        let opts = ScanOptions::default();
        let roots = find_all_roots(|x| 1.0 / (x - 1.0), &opts);
        assert!(roots.is_empty(), "pole was reported as a root: {:?}", roots);
    }

    #[test]
    fn scan_rejects_tangent_poles() {
        // tan(x) has poles at +/- pi/2; the roots are at multiples of pi.
        let opts = ScanOptions { lo: -4.0, hi: 4.0, ..Default::default() };
        let roots = find_all_roots(|x| x.tan(), &opts);
        assert_eq!(roots.len(), 3, "got {:?}", roots);
        for r in &roots {
            let k = r / std::f64::consts::PI;
            assert!((k - k.round()).abs() < 1e-6, "root {} is not a multiple of pi", r);
        }
    }

    #[test]
    fn scan_finds_root_beside_a_pole() {
        // 1/(x-1) - 2 = 0 has the root x = 1.5, right next to the pole.
        let opts = ScanOptions::default();
        let roots = find_all_roots(|x| 1.0 / (x - 1.0) - 2.0, &opts);
        assert_eq!(roots.len(), 1, "got {:?}", roots);
        assert!((roots[0] - 1.5).abs() < 1e-6, "got {}", roots[0]);
    }

    #[test]
    fn scan_respects_bounds() {
        // The root at pi is outside [-1, 1].
        let opts = ScanOptions { lo: -1.0, hi: 1.0, ..Default::default() };
        let roots = find_all_roots(|x| x.sin(), &opts);
        assert_eq!(roots.len(), 1, "got {:?}", roots);
        assert!(roots[0].abs() < 1e-6);
    }

    #[test]
    fn scan_returns_sorted_roots() {
        let opts = ScanOptions::default();
        let roots = find_all_roots(|x| x.sin(), &opts);
        for w in roots.windows(2) {
            assert!(w[0] < w[1], "roots not sorted: {:?}", roots);
        }
    }

    #[test]
    fn scan_handles_nan_regions() {
        // log(x) is undefined for x <= 0 — the scan must skip that region
        // rather than produce a spurious root there. log(x) = 0 at x = 1,
        // which is a genuine root.
        let opts = ScanOptions { lo: -5.0, hi: 5.0, ..Default::default() };
        let roots = find_all_roots(|x| if x > 0.0 { x.ln() } else { f64::NAN }, &opts);
        assert_eq!(roots.len(), 1, "got {:?}", roots);
        assert!((roots[0] - 1.0).abs() < 1e-6, "got {}", roots[0]);
        // No root may be reported in the undefined region.
        assert!(roots.iter().all(|r| *r > 0.0), "got {:?}", roots);
    }

    #[test]
    fn scan_deduplicates_near_duplicates() {
        // A function that oscillates fast can bracket the same root twice.
        let opts = ScanOptions { lo: -1.0, hi: 1.0, samples: 2000, ..Default::default() };
        let roots = find_all_roots(|x| (x * 500.0).sin() * (x - 0.25), &opts);
        // The root at 0.25 must appear exactly once.
        let near = roots.iter().filter(|r| (*r - 0.25).abs() < 1e-6).count();
        assert_eq!(near, 1, "root 0.25 reported {} times: {:?}", near, roots);
    }
}
