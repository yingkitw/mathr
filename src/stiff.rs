//! Implicit (A-stable) ODE solvers for stiff systems.
//!
//! Complements [`crate::ode`] (Euler, RK4, RKF45 — all explicit, all with a
//! bounded stability region) with two classic implicit methods:
//!
//! - **Backward (implicit) Euler**: `y_{n+1} = y_n + h·f(t_{n+1}, y_{n+1})`,
//!   first order, L-stable — damping and robust on the stiffest transients.
//! - **Implicit trapezoidal** (Crank–Nicolson):
//!   `y_{n+1} = y_n + h/2·(f(t_n, y_n) + f(t_{n+1}, y_{n+1}))`, second order,
//!   A-stable — accurate long-time integration with mild oscillation
//!   damping (not L-stable).
//! - **BDF2** (two-step backward differentiation formula):
//!   `(3y_{n+1} − 4y_n + y_{n−1})/(2h) = f(t_{n+1}, y_{n+1})`, second order,
//!   A-stable — the highest-order A-stable BDF. Fixed step with a single
//!   trapezoidal startup step.
//!
//! Each step solves the implicit equation with Newton's method using a
//! numeric Jacobian of `f` w.r.t. the state (central differences) and the
//! public API of [`crate::matrix::Matrix`] for the linear solves. The
//! Newton iteration starts from the previous state, which converges in a
//! couple of iterations for smooth problems.
//!
//! # Example
//!
//! ```
//! use mathr::stiff::implicit_euler_system;
//!
//! // y' = -1000 (y - 1) is stiff: explicit Euler needs h < 0.002,
//! // backward Euler is unconditionally stable.
//! let f = |_t: f64, y: &[f64]| vec![-1000.0 * (y[0] - 1.0)];
//! let y = implicit_euler_system(&f, 0.0, 1.0, &[0.0], 100).unwrap();
//! assert!((y[0] - 1.0).abs() < 1e-6);
//! ```

use crate::error::{MathError, Result};
use crate::matrix::Matrix;

/// Newton convergence tolerance (relative to state scale).
const NEWTON_TOL: f64 = 1e-10;
/// Maximum Newton iterations per time step.
const NEWTON_MAX_IT: usize = 50;

/// Central-difference Jacobian of `f(t, y)` w.r.t. `y`, row-major m×m.
fn numeric_jacobian<F: Fn(f64, &[f64]) -> Vec<f64>>(f: &F, t: f64, y: &[f64]) -> Result<Matrix> {
    let m = y.len();
    let mut jac = vec![0.0; m * m];
    for j in 0..m {
        let h = 1e-6 * y[j].abs().max(1.0);
        let mut yp = y.to_vec();
        let mut ym = y.to_vec();
        yp[j] += h;
        ym[j] -= h;
        let fp = f(t, &yp);
        let fm = f(t, &ym);
        if fp.iter().any(|v| !v.is_finite()) || fm.iter().any(|v| !v.is_finite()) {
            return Err(MathError::Domain(format!(
                "stiff solver: f produced non-finite values at t = {t}"
            )));
        }
        for i in 0..m {
            jac[i * m + j] = (fp[i] - fm[i]) / (2.0 * h);
        }
    }
    Matrix::from_row_major(m, m, jac)
}

/// Solve `G(y) = y − base − hb·f(t_new, y) = 0` by Newton's method,
/// starting from `y` (mutated in place to the solution).
fn newton_solve<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t_new: f64,
    base: &[f64],
    hb: f64,
    y: &mut [f64],
) -> Result<()> {
    let m = y.len();
    for _ in 0..NEWTON_MAX_IT {
        let fy = f(t_new, y);
        if fy.iter().any(|v| !v.is_finite()) {
            return Err(MathError::Domain(format!(
                "stiff solver: f produced non-finite values at t = {t_new}"
            )));
        }
        // G = y - base - hb*fy
        let mut g = vec![0.0; m];
        for i in 0..m {
            g[i] = y[i] - base[i] - hb * fy[i];
        }
        let scale = 1.0 + y.iter().fold(0.0f64, |mx, v| mx.max(v.abs()));
        if g.iter().fold(0.0f64, |mx, v| mx.max(v.abs())) <= NEWTON_TOL * scale {
            return Ok(());
        }
        // J_G = I - hb * J_f
        let jf = numeric_jacobian(f, t_new, y)?;
        let mut jg = vec![0.0; m * m];
        for i in 0..m {
            for j in 0..m {
                jg[i * m + j] = -hb * jf[(i, j)];
            }
            jg[i * m + i] += 1.0;
        }
        let neg_g: Vec<f64> = g.iter().map(|v| -v).collect();
        let delta = Matrix::from_row_major(m, m, jg)?.solve(&neg_g)?;
        let step_inf = delta.iter().fold(0.0f64, |mx, v| mx.max(v.abs()));
        for (yi, di) in y.iter_mut().zip(delta.iter()) {
            *yi += di;
        }
        if step_inf <= NEWTON_TOL * scale {
            // final residual check next loop costs one extra f eval; accept
            return Ok(());
        }
    }
    Err(MathError::NotConvergent(
        "stiff solver: Newton iteration did not converge".into(),
    ))
}

/// Integrate a step of backward Euler: solve
/// `y_new − y_prev − h·f(t_new, y_new) = 0`.
fn be_step<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t_old: f64,
    h: f64,
    y: &mut [f64],
) -> Result<()> {
    let base = y.to_vec();
    let t_new = t_old + h;
    newton_solve(f, t_new, &base, h, y)
}

/// Integrate a step of implicit trapezoidal: solve
/// `y_new − y_prev − h/2·(f(t_old, y_prev) + f(t_new, y_new)) = 0`.
fn trap_step<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t_old: f64,
    h: f64,
    y: &mut [f64],
) -> Result<()> {
    let f_old = f(t_old, y);
    let base: Vec<f64> = y
        .iter()
        .zip(f_old.iter())
        .map(|(a, b)| a + 0.5 * h * b)
        .collect();
    let t_new = t_old + h;
    newton_solve(f, t_new, &base, 0.5 * h, y)
}

macro_rules! integrate_impl {
    ($fname:ident, $tname:ident, $doc:expr) => {
        #[doc = $doc]
        pub fn $fname<F: Fn(f64, &[f64]) -> Vec<f64>>(
            f: &F,
            t0: f64,
            t1: f64,
            y0: &[f64],
            n: usize,
        ) -> Result<Vec<f64>> {
            let traj = $tname(f, t0, t1, y0, n)?;
            Ok(traj.last().unwrap().1.clone())
        }
    };
}

/// Full trajectory `Vec<(t, state)>` of the backward Euler method for the
/// system `dy/dt = f(t, y)` from `t0` to `t1` with `n` steps.
pub fn implicit_euler_trajectory<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n: usize,
) -> Result<Vec<(f64, Vec<f64>)>> {
    trajectory_impl(f, t0, t1, y0, n, be_step, "implicit_euler")
}

/// Full trajectory `Vec<(t, state)>` of the implicit trapezoidal method
/// (Crank–Nicolson) for `dy/dt = f(t, y)` from `t0` to `t1` with `n` steps.
pub fn trapezoid_trajectory<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n: usize,
) -> Result<Vec<(f64, Vec<f64>)>> {
    trajectory_impl(f, t0, t1, y0, n, trap_step, "trapezoid")
}

/// Full trajectory `Vec<(t, state)>` of the fixed-step BDF2 method for
/// `dy/dt = f(t, y)` from `t0` to `t1` with `n` steps.
///
/// The first step is implicit trapezoidal (second-order startup); the
/// remaining steps apply the two-step backward differentiation formula
/// `(3y_{n+1} − 4y_n + y_{n−1})/(2h) = f(t_{n+1}, y_{n+1})`, solved by
/// Newton iteration. A-stable and second-order accurate.
pub fn bdf2_trajectory<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n: usize,
) -> Result<Vec<(f64, Vec<f64>)>> {
    if n == 0 {
        return Err(MathError::InvalidArgument("bdf2: n must be > 0".into()));
    }
    if y0.is_empty() || y0.iter().any(|v| !v.is_finite()) {
        return Err(MathError::InvalidArgument(
            "bdf2: y0 must be non-empty and finite".into(),
        ));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut out = Vec::with_capacity(n + 1);
    out.push((t, y.clone()));

    // Startup: one trapezoidal step (second-order local error).
    trap_step(f, t, h, &mut y)?;
    if y.iter().any(|v| !v.is_finite()) {
        return Err(MathError::NotConvergent(format!(
            "bdf2: state became non-finite at t = {t}"
        )));
    }
    t += h;
    out.push((t, y.clone()));

    // BDF2: y − (4y_n − y_{n−1})/3 − (2h/3)·f(t_new, y) = 0.
    let mut y_prev = out[out.len() - 2].1.clone();
    for _ in 1..n {
        let y_curr = y.clone();
        let base: Vec<f64> = y
            .iter()
            .zip(&y_prev)
            .map(|(a, b)| (4.0 * a - b) / 3.0)
            .collect();
        newton_solve(f, t + h, &base, 2.0 * h / 3.0, &mut y)?;
        if y.iter().any(|v| !v.is_finite()) {
            return Err(MathError::NotConvergent(format!(
                "bdf2: state became non-finite at t = {t}"
            )));
        }
        y_prev = y_curr;
        t += h;
        out.push((t, y.clone()));
    }
    Ok(out)
}

fn trajectory_impl<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: &F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n: usize,
    step: fn(&F, f64, f64, &mut [f64]) -> Result<()>,
    name: &str,
) -> Result<Vec<(f64, Vec<f64>)>> {
    if n == 0 {
        return Err(MathError::InvalidArgument(format!("{name}: n must be > 0")));
    }
    if y0.is_empty() || y0.iter().any(|v| !v.is_finite()) {
        return Err(MathError::InvalidArgument(format!(
            "{name}: y0 must be non-empty and finite"
        )));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0.to_vec();
    let mut out = Vec::with_capacity(n + 1);
    out.push((t, y.clone()));
    for _ in 0..n {
        step(f, t, h, &mut y)?;
        if y.iter().any(|v| !v.is_finite()) {
            return Err(MathError::NotConvergent(format!(
                "{name}: state became non-finite at t = {t}"
            )));
        }
        t += h;
        out.push((t, y.clone()));
    }
    Ok(out)
}

integrate_impl!(
    implicit_euler_system,
    implicit_euler_trajectory,
    "Integrate the system `dy/dt = f(t, y)` from `t0` to `t1` with backward\n(Euler) and return the state at `t1`. A-stable (L-stable): handles stiff\ndecays with step sizes far beyond the explicit stability limit."
);
integrate_impl!(
    trapezoid_system,
    trapezoid_trajectory,
    "Integrate the system `dy/dt = f(t, y)` from `t0` to `t1` with the\nimplicit trapezoidal rule (Crank–Nicolson) and return the state at `t1`.\nA-stable and second-order accurate."
);
integrate_impl!(
    bdf2_system,
    bdf2_trajectory,
    "Integrate the system `dy/dt = f(t, y)` from `t0` to `t1` with the\nfixed-step BDF2 method (trapezoidal startup, two-step backward\ndifferentiation) and return the state at `t1`. A-stable and\nsecond-order accurate."
);

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * (1.0 + b.abs().max(a.abs()))
    }

    /// y' = -1000 (y - 1), y(0) = 0 -> y(t) = 1 - e^(-1000 t) ~ 1 at t = 1
    fn stiff_f(_t: f64, y: &[f64]) -> Vec<f64> {
        vec![-1000.0 * (y[0] - 1.0)]
    }

    #[test]
    fn backward_euler_survives_stiffness() {
        // h = 0.01 gives h*k = 10 >> 2 (explicit stability limit);
        // the BE amplification (1/(1+10)) converges monotonically to 1
        let y = implicit_euler_system(&stiff_f, 0.0, 1.0, &[0.0], 100).unwrap();
        assert!(y[0].is_finite());
        assert!(close(y[0], 1.0, 1e-6), "y={}", y[0]);
    }

    #[test]
    fn explicit_euler_blows_up_on_same_problem() {
        // contrast test: explicit Euler with the same h oscillates with
        // amplification |1 - h*k| = 9 per step: |y_100| = 9^100 ~ 1e95
        let r = crate::ode::euler(|_t, y| -1000.0 * (y - 1.0), 0.0, 1.0, 0.0, 100).unwrap();
        assert!(r.abs() > 1e10 || !r.is_finite(), "expected blow-up, got {r}");
    }

    #[test]
    fn trapezoid_survives_stiffness() {
        // amplification (1-5)/(1+5) = -2/3 per step: bounded, converges to 1
        let y = trapezoid_system(&stiff_f, 0.0, 1.0, &[0.0], 100).unwrap();
        assert!(y[0].is_finite());
        assert!(close(y[0], 1.0, 1e-3), "y={}", y[0]);
    }

    #[test]
    fn backward_euler_first_order() {
        // y' = y, y(0) = 1: error at t=1 halves when h halves (O(h))
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let e1 = (implicit_euler_system(&f, 0.0, 1.0, &[1.0], 100).unwrap()[0]
            - std::f64::consts::E)
            .abs();
        let e2 = (implicit_euler_system(&f, 0.0, 1.0, &[1.0], 200).unwrap()[0]
            - std::f64::consts::E)
            .abs();
        let ratio = e1 / e2;
        assert!(ratio > 1.8 && ratio < 2.2, "order-1 ratio={ratio}");
    }

    #[test]
    fn trapezoid_second_order() {
        // y' = y: error ratio ~ 4 when h halves (O(h^2))
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let e1 = (trapezoid_system(&f, 0.0, 1.0, &[1.0], 100).unwrap()[0]
            - std::f64::consts::E)
            .abs();
        let e2 = (trapezoid_system(&f, 0.0, 1.0, &[1.0], 200).unwrap()[0]
            - std::f64::consts::E)
            .abs();
        let ratio = e1 / e2;
        assert!(ratio > 3.5 && ratio < 4.5, "order-2 ratio={ratio}");
    }

    #[test]
    fn bdf2_second_order() {
        // y' = y: error ratio ~ 4 when h halves (O(h^2)); the single
        // trapezoidal startup step is also second-order local, so the
        // global order stays 2.
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let e1 = (bdf2_system(&f, 0.0, 1.0, &[1.0], 100).unwrap()[0] - std::f64::consts::E).abs();
        let e2 = (bdf2_system(&f, 0.0, 1.0, &[1.0], 200).unwrap()[0] - std::f64::consts::E).abs();
        let ratio = e1 / e2;
        assert!(ratio > 3.5 && ratio < 4.5, "order-2 ratio={ratio}");
    }

    #[test]
    fn bdf2_survives_stiffness() {
        // kh = 10: BDF2 amplification |r| = sqrt(0.5/11.5) ~ 0.21 per step
        // (bounded decay), startup trapezoid |r| = 2/3 -> converges to 1.
        let y = bdf2_system(&stiff_f, 0.0, 1.0, &[0.0], 100).unwrap();
        assert!(y[0].is_finite());
        assert!(close(y[0], 1.0, 1e-6), "y={}", y[0]);
    }

    #[test]
    fn bdf2_harmonic_oscillator_stays_bounded() {
        // A-stable BDF2: the principal root has |r| ~ 1 + O((h w)^2)
        // on the imaginary axis — growth over one period is negligible.
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let y = bdf2_system(&f, 0.0, 2.0 * std::f64::consts::PI, &[1.0, 0.0], 200).unwrap();
        assert!(close(y[0], 1.0, 1e-2), "y0={}", y[0]);
        assert!(close(y[1], 0.0, 1e-2), "y1={}", y[1]);
    }

    #[test]
    fn trapezoid_harmonic_oscillator_conserves_amplitude() {
        // A-stable trapezoidal rule has |amplification| = 1 for the
        // harmonic oscillator: the oscillation does not decay or blow up
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let y = trapezoid_system(&f, 0.0, 2.0 * std::f64::consts::PI, &[1.0, 0.0], 200).unwrap();
        assert!(close(y[0], 1.0, 1e-2), "y0={}", y[0]);
        assert!(close(y[1], 0.0, 1e-2), "y1={}", y[1]);
    }

    #[test]
    fn robertson_system_stays_bounded_and_conserves_mass() {
        // Robertson problem: classic stiff nonlinear system with rate
        // constants 4e-2, 1e4, 3e7. Invariants: y1+y2+y3 = 1 exactly,
        // all components non-negative, y2 stays small.
        let f = |_t: f64, y: &[f64]| {
            vec![
                -0.04 * y[0] + 1.0e4 * y[1] * y[2],
                0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1],
                3.0e7 * y[1] * y[1],
            ]
        };
        for solver in [implicit_euler_system, trapezoid_system, bdf2_system] {
            let y = solver(&f, 0.0, 10.0, &[1.0, 0.0, 0.0], 1000).unwrap();
            assert!(y.iter().all(|v| v.is_finite() && *v >= 0.0), "{:?}", y);
            let total: f64 = y.iter().sum();
            assert!(close(total, 1.0, 1e-6), "sum={total}");
            assert!(y[1] < 1.0e-3, "y2={}", y[1]);
            assert!(y[2] > y[1], "y3 should dominate y2: {:?}", y);
        }
    }

    #[test]
    fn trajectory_shape() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        let traj = implicit_euler_trajectory(&f, 0.0, 1.0, &[1.0], 10).unwrap();
        assert_eq!(traj.len(), 11);
        assert!(close(traj[0].0, 0.0, 1e-12));
        assert!(close(traj[10].0, 1.0, 1e-12));
        let traj2 = trapezoid_trajectory(&f, 0.0, 1.0, &[1.0], 10).unwrap();
        assert_eq!(traj2.len(), 11);
        let traj3 = bdf2_trajectory(&f, 0.0, 1.0, &[1.0], 10).unwrap();
        assert_eq!(traj3.len(), 11);
        assert!(close(traj3[0].0, 0.0, 1e-12));
        assert!(close(traj3[10].0, 1.0, 1e-12));
    }

    #[test]
    fn validation_errors() {
        let f = |_t: f64, y: &[f64]| vec![y[0]];
        assert!(implicit_euler_system(&f, 0.0, 1.0, &[1.0], 0).is_err());
        assert!(implicit_euler_trajectory(&f, 0.0, 1.0, &[], 10).is_err());
        assert!(trapezoid_system(&f, 0.0, 1.0, &[f64::NAN], 10).is_err());
        assert!(bdf2_system(&f, 0.0, 1.0, &[1.0], 0).is_err());
        assert!(bdf2_trajectory(&f, 0.0, 1.0, &[], 10).is_err());
        assert!(bdf2_system(&f, 0.0, 1.0, &[f64::NAN], 10).is_err());
    }
}
