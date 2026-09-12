//! Derivative-free function minimization.
//!
//! Complements [`crate::calculus`] (derivatives/gradients) with two
//! classic minimizers that need no derivatives:
//!
//! - **Golden-section search**: 1-D minimization of a unimodal function on
//!   a bracket `[a, b]`. Each iteration shrinks the bracket by the golden
//!   ratio `φ = (√5 − 1)/2`, so the final bracket width is `b − a` times
//!   `φ^n` — machine-precision in ~90 iterations from any reasonable start.
//! - **Nelder–Mead simplex**: N-dimensional minimization by reflecting,
//!   expanding, contracting, and shrinking a simplex of `n + 1` vertices.
//!   Robust on non-smooth and noisy objectives; converges to the exact
//!   minimum on strictly convex problems (e.g. Rosenbrock).
//!
//! # Example
//!
//! ```
//! use mathr::optim::{golden_section, nelder_mead, OptOptions};
//!
//! // 1-D: (x - 2)^2 has its minimum at x = 2
//! let (x, fx) = golden_section(|x| (x - 2.0) * (x - 2.0), 0.0, 5.0, 1e-10).unwrap();
//! assert!((x - 2.0).abs() < 1e-8);
//!
//! // N-D: Rosenbrock's banana valley, minimum (1, 1)
//! let rosen = |p: &[f64]| {
//!     let (a, b) = (p[0] - 1.0, p[0] * p[0] - p[1]);
//!     a * a + 100.0 * b * b
//! };
//! let r = nelder_mead(&rosen, &[-1.2, 1.0], &OptOptions::default()).unwrap();
//! assert!(r.converged);
//! assert!((r.x[0] - 1.0).abs() < 1e-4 && (r.x[1] - 1.0).abs() < 1e-4);
//! ```

use crate::error::{MathError, Result};

/// Tuning parameters for [`nelder_mead`].
#[derive(Debug, Clone)]
pub struct OptOptions {
    /// Maximum iterations (default 1000).
    pub max_iter: usize,
    /// Convergence tolerance on the simplex diameter (default 1e-10).
    pub tol: f64,
    /// Initial simplex step per coordinate: `start_i + step_i` (default
    /// `0.05 * max(1, |start_i|) + 0.0005`, mirroring common practice).
    pub init_step: Option<Vec<f64>>,
}

impl Default for OptOptions {
    fn default() -> Self {
        OptOptions { max_iter: 1000, tol: 1e-10, init_step: None }
    }
}

/// Outcome of a minimization run; `converged = false` when `max_iter` was
/// exhausted (the best point found is still returned).
#[derive(Debug, Clone, PartialEq)]
pub struct OptResult {
    /// Minimizer coordinates.
    pub x: Vec<f64>,
    /// Objective value at `x`.
    pub fx: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// True if the simplex diameter fell below `tol` in time.
    pub converged: bool,
}

/// Minimize a unimodal function `f` on the bracket `[a, b]` by
/// golden-section search. Returns `(x*, f(x*))`.
pub fn golden_section<F: FnMut(f64) -> f64>(mut f: F, a: f64, b: f64, tol: f64) -> Result<(f64, f64)> {
    if a >= b || !a.is_finite() || !b.is_finite() {
        return Err(MathError::InvalidArgument(format!(
            "golden_section: bad bracket [{a}, {b}]"
        )));
    }
    if tol <= 0.0 {
        return Err(MathError::InvalidArgument("golden_section: tol must be > 0".into()));
    }
    let phi = (5.0f64.sqrt() - 1.0) / 2.0;
    let mut lo = a;
    let mut hi = b;
    let mut c = hi - phi * (hi - lo);
    let mut d = lo + phi * (hi - lo);
    let mut fc = f(c);
    let mut fd = f(d);
    while hi - lo > tol {
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - phi * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + phi * (hi - lo);
            fd = f(d);
        }
    }
    let x = 0.5 * (lo + hi);
    Ok((x, f(x)))
}

/// Minimize `f: R^n -> R` from `start` with the Nelder–Mead simplex method
/// (reflection α=1, expansion γ=2, contraction ρ=0.5, shrink σ=0.5).
pub fn nelder_mead<F: Fn(&[f64]) -> f64>(
    f: F,
    start: &[f64],
    opts: &OptOptions,
) -> Result<OptResult> {
    if start.is_empty() {
        return Err(MathError::InvalidArgument("nelder_mead: empty start point".into()));
    }
    if start.iter().any(|v| !v.is_finite()) {
        return Err(MathError::InvalidArgument("nelder_mead: non-finite start point".into()));
    }
    let n = start.len();

    // initial simplex: start plus one perturbed vertex per coordinate
    let step: Vec<f64> = match &opts.init_step {
        Some(s) if s.len() == n => s.clone(),
        Some(_) => {
            return Err(MathError::InvalidArgument(
                "nelder_mead: init_step length must equal start length".into(),
            ));
        }
        None => start
            .iter()
            .map(|&s| 0.05 * s.abs().max(1.0) + 0.0005)
            .collect(),
    };
    let mut simplex: Vec<Vec<f64>> = Vec::with_capacity(n + 1);
    simplex.push(start.to_vec());
    for j in 0..n {
        let mut v = start.to_vec();
        v[j] += step[j];
        simplex.push(v);
    }
    let mut values: Vec<f64> = simplex.iter().map(|v| f(v)).collect();
    if values.iter().any(|v| !v.is_finite()) {
        return Err(MathError::Domain(
            "nelder_mead: objective produced non-finite values".into(),
        ));
    }

    let centroid = |verts: &[Vec<f64>]| -> Vec<f64> {
        let mut c = vec![0.0; n];
        for v in &verts[..n] {
            // exclude the worst vertex (last after sort)
            for (ci, cv) in c.iter_mut().enumerate() {
                *cv += v[ci];
            }
        }
        c.iter().map(|v| v / n as f64).collect()
    };

    let mut iterations = 0usize;
    let mut converged = false;

    while iterations < opts.max_iter {
        iterations += 1;
        // sort vertices by objective: best first, worst last
        let mut order: Vec<usize> = (0..=n).collect();
        order.sort_by(|&a, &b| values[a].partial_cmp(&values[b]).unwrap());
        let sorted_s: Vec<Vec<f64>> = order.iter().map(|&i| simplex[i].clone()).collect();
        let sorted_v: Vec<f64> = order.iter().map(|&i| values[i]).collect();
        simplex = sorted_s;
        values = sorted_v;

        // convergence: simplex diameter (max edge length in any coordinate)
        let mut diameter = 0.0f64;
        for v in &simplex[1..] {
            for (a, b) in v.iter().zip(simplex[0].iter()) {
                diameter = diameter.max((a - b).abs());
            }
        }
        if diameter <= opts.tol {
            converged = true;
            break;
        }

        let cen = centroid(&simplex);
        let reflect: Vec<f64> = cen
            .iter()
            .zip(simplex[n].iter())
            .map(|(c, w)| c + (c - w))
            .collect();
        let fr = f(&reflect);
        if fr < values[0] {
            // try expansion
            let expand: Vec<f64> = cen
                .iter()
                .zip(simplex[n].iter())
                .map(|(c, w)| c + 2.0 * (c - w))
                .collect();
            let fe = f(&expand);
            if fe < fr {
                simplex[n] = expand;
                values[n] = fe;
            } else {
                simplex[n] = reflect;
                values[n] = fr;
            }
        } else if fr < values[n - 1] {
            // accept reflection (better than second-worst)
            simplex[n] = reflect;
            values[n] = fr;
        } else {
            // contraction (outside if better than worst, inside otherwise)
            let outside = fr < values[n];
            let contract: Vec<f64> = cen
                .iter()
                .zip(simplex[n].iter())
                .map(|(c, w)| {
                    if outside {
                        c + 0.5 * (c - w)
                    } else {
                        c - 0.5 * (c - w)
                    }
                })
                .collect();
            let fc = f(&contract);
            if (outside && fc < fr) || (!outside && fc < values[n]) {
                simplex[n] = contract;
                values[n] = fc;
            } else {
                // shrink everything toward the best vertex
                let best = simplex[0].clone();
                for v in simplex.iter_mut().skip(1) {
                    for (vi, b) in v.iter_mut().zip(best.iter()) {
                        *vi = *b + 0.5 * (*vi - *b);
                    }
                }
                for k in 1..=n {
                    values[k] = f(&simplex[k]);
                }
            }
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(MathError::Domain(
                "nelder_mead: objective produced non-finite values".into(),
            ));
        }
    }

    // best vertex may be anywhere after the last partial sort
    let best = (0..=n)
        .enumerate()
        .min_by(|(_, a), (_, b)| values[*a].partial_cmp(&values[*b]).unwrap())
        .map(|(i, _)| i)
        .unwrap_or(0);
    Ok(OptResult {
        x: simplex[best].clone(),
        fx: values[best],
        iterations,
        converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs().max(a.abs()))
    }

    #[test]
    fn golden_section_parabola() {
        // (x - 2)^2: known minimum x = 2 (closed form)
        let (x, fx) = golden_section(|x| (x - 2.0) * (x - 2.0), 0.0, 5.0, 1e-10).unwrap();
        assert!(close(x, 2.0, 1e-8), "x={x}");
        assert!(fx <= 1e-16);
    }

    #[test]
    fn golden_section_sine() {
        // sin has its minimum at 3pi/2 on [pi, 2pi] (closed form).
        // Near a flat minimum f' = 0, comparisons between f(c) and f(d)
        // are noise-limited, so accuracy is ~1e-8, not the bracket tol.
        let (x, _) = golden_section(
            |x| x.sin(),
            std::f64::consts::PI,
            2.0 * std::f64::consts::PI,
            1e-12,
        )
        .unwrap();
        assert!(close(x, 1.5 * std::f64::consts::PI, 1e-7), "x={x}");
    }

    #[test]
    fn golden_section_bracket_shrinks_by_phi() {
        // after n iterations the bracket width is (b-a)*phi^n
        let calls = std::cell::Cell::new(0usize);
        let _ = golden_section(
            |x| {
                calls.set(calls.get() + 1);
                (x - 1.0) * (x - 1.0)
            },
            0.0,
            1.0,
            1e-3,
        )
        .unwrap();
        // (0.618)^n < 1e-3 -> n >= 14; each iteration costs 1 new evaluation
        assert!(calls.get() >= 14 && calls.get() < 40, "calls={}", calls.get());
    }

    #[test]
    fn golden_section_validation() {
        assert!(golden_section(|x| x * x, 2.0, 1.0, 1e-6).is_err());
        assert!(golden_section(|x| x * x, 0.0, f64::NAN, 1e-6).is_err());
        assert!(golden_section(|x| x * x, 0.0, 1.0, 0.0).is_err());
    }

    #[test]
    fn nelder_mead_quadratic() {
        // f = (x-3)^2 + 2(y+1)^2: minimum (3, -1), value 0
        let f = |p: &[f64]| (p[0] - 3.0) * (p[0] - 3.0) + 2.0 * (p[1] + 1.0) * (p[1] + 1.0);
        let r = nelder_mead(&f, &[0.0, 0.0], &OptOptions::default()).unwrap();
        assert!(r.converged);
        assert!(close(r.x[0], 3.0, 1e-6), "x={:?}", r.x);
        assert!(close(r.x[1], -1.0, 1e-6), "y={:?}", r.x);
        assert!(r.fx < 1e-12);
    }

    #[test]
    fn nelder_mead_rosenbrock() {
        // classic banana valley: minimum (1, 1) with value 0
        let rosen = |p: &[f64]| {
            let (a, b) = (p[0] - 1.0, p[0] * p[0] - p[1]);
            a * a + 100.0 * b * b
        };
        let r = nelder_mead(&rosen, &[-1.2, 1.0], &OptOptions::default()).unwrap();
        assert!(r.converged, "iters={}", r.iterations);
        assert!(close(r.x[0], 1.0, 1e-4), "x={:?}", r.x);
        assert!(close(r.x[1], 1.0, 1e-4), "y={:?}", r.x);
        assert!(r.fx < 1e-8, "fx={}", r.fx);
    }

    #[test]
    fn nelder_mead_higher_dimensional() {
        // f(x) = sum (xi - i)^2: minimum (1, 2, 3)
        let f = |p: &[f64]| {
            (0..p.len())
                .map(|i| (p[i] - (i + 1) as f64) * (p[i] - (i + 1) as f64))
                .sum::<f64>()
        };
        let r = nelder_mead(&f, &[0.0, 0.0, 0.0], &OptOptions::default()).unwrap();
        assert!(r.converged);
        for (i, xi) in r.x.iter().enumerate() {
            assert!(close(*xi, (i + 1) as f64, 1e-5), "x[{i}]={xi}");
        }
    }

    #[test]
    fn nelder_mead_max_iter_reports_not_converged() {
        // one-dimensional with a tiny iteration budget
        let f = |p: &[f64]| (p[0] - 1.0) * (p[0] - 1.0);
        let opts = OptOptions { max_iter: 1, ..OptOptions::default() };
        let r = nelder_mead(&f, &[0.0], &opts).unwrap();
        assert!(!r.converged);
        assert_eq!(r.iterations, 1);
    }

    #[test]
    fn nelder_mead_validation() {
        let f = |p: &[f64]| p.iter().sum();
        assert!(nelder_mead(&f, &[], &OptOptions::default()).is_err());
        assert!(nelder_mead(&f, &[f64::NAN], &OptOptions::default()).is_err());
        let opts = OptOptions { init_step: Some(vec![1.0, 1.0]), ..OptOptions::default() };
        assert!(nelder_mead(&f, &[1.0], &opts).is_err()); // step length mismatch
    }
}
