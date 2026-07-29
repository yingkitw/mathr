//! Ordinary differential equation (ODE) solvers from scratch.
//!
//! Provides Euler's method, the classical 4th-order Runge–Kutta (RK4),
//! and adaptive step-size control via the Runge–Kutta–Fehlberg (RKF45) method.

use crate::error::{MathError, Result};

/// One step of Euler's method: `y_{n+1} = y_n + h * f(t_n, y_n)`.
pub fn euler_step<F: Fn(f64, f64) -> f64>(f: F, t: f64, y: f64, h: f64) -> f64 {
    y + h * f(t, y)
}

/// Integrate `dy/dt = f(t, y)` from `t0` to `t1` with initial condition `y0`
/// using Euler's method with `n` steps.
pub fn euler<F: Fn(f64, f64) -> f64>(f: F, t0: f64, t1: f64, y0: f64, n: usize) -> Result<f64> {
    if n == 0 {
        return Err(MathError::InvalidArgument("euler: n must be > 0".into()));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0;
    for _ in 0..n {
        y = euler_step(&f, t, y, h);
        t += h;
    }
    Ok(y)
}

/// One step of classical 4th-order Runge–Kutta (RK4).
pub fn rk4_step<F: Fn(f64, f64) -> f64>(f: F, t: f64, y: f64, h: f64) -> f64 {
    let k1 = f(t, y);
    let k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
    let k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
    let k4 = f(t + h, y + h * k3);
    y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

/// Integrate `dy/dt = f(t, y)` from `t0` to `t1` with initial condition `y0`
/// using RK4 with `n` steps.
pub fn rk4<F: Fn(f64, f64) -> f64>(f: F, t0: f64, t1: f64, y0: f64, n: usize) -> Result<f64> {
    if n == 0 {
        return Err(MathError::InvalidArgument("rk4: n must be > 0".into()));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0;
    for _ in 0..n {
        y = rk4_step(&f, t, y, h);
        t += h;
    }
    Ok(y)
}

/// Integrate and return the full trajectory as `Vec<(t, y)>`.
pub fn rk4_trajectory<F: Fn(f64, f64) -> f64>(
    f: F,
    t0: f64,
    t1: f64,
    y0: f64,
    n: usize,
) -> Result<Vec<(f64, f64)>> {
    if n == 0 {
        return Err(MathError::InvalidArgument("rk4_trajectory: n must be > 0".into()));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0;
    let mut out = Vec::with_capacity(n + 1);
    out.push((t, y));
    for _ in 0..n {
        y = rk4_step(&f, t, y, h);
        t += h;
        out.push((t, y));
    }
    Ok(out)
}

/// RK4 for a system of ODEs: `dy_i/dt = f_i(t, y)`.
/// `y0` is the initial state vector; returns the state at `t1`.
pub fn rk4_system<F: Fn(f64, &[f64]) -> Vec<f64>>(
    f: F,
    t0: f64,
    t1: f64,
    y0: &[f64],
    n: usize,
) -> Result<Vec<f64>> {
    if n == 0 {
        return Err(MathError::InvalidArgument("rk4_system: n must be > 0".into()));
    }
    let h = (t1 - t0) / n as f64;
    let mut t = t0;
    let mut y = y0.to_vec();
    let dim = y.len();
    for _ in 0..n {
        let k1 = f(t, &y);
        let mut y2 = vec![0.0; dim];
        for i in 0..dim {
            y2[i] = y[i] + 0.5 * h * k1[i];
        }
        let k2 = f(t + 0.5 * h, &y2);
        let mut y3 = vec![0.0; dim];
        for i in 0..dim {
            y3[i] = y[i] + 0.5 * h * k2[i];
        }
        let k3 = f(t + 0.5 * h, &y3);
        let mut y4 = vec![0.0; dim];
        for i in 0..dim {
            y4[i] = y[i] + h * k3[i];
        }
        let k4 = f(t + h, &y4);
        for i in 0..dim {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        t += h;
    }
    Ok(y)
}

/// Adaptive RK4(5) — Runge–Kutta–Fehlberg — with automatic step-size control.
/// Integrates `dy/dt = f(t, y)` from `t0` to `t1` with initial condition `y0`.
/// `tol` is the desired absolute error per step.
pub fn rkf45<F: Fn(f64, f64) -> f64>(
    f: F,
    t0: f64,
    t1: f64,
    y0: f64,
    tol: f64,
) -> Result<f64> {
    let mut t = t0;
    let mut y = y0;
    let mut h = (t1 - t0) / 100.0;
    let direction = if t1 > t0 { 1.0 } else { -1.0 };

    while (t - t1).abs() > 1e-12 {
        if (t + h - t1) * direction > 0.0 {
            h = t1 - t;
        }

        let k1 = f(t, y);
        let k2 = f(t + h / 4.0, y + h / 4.0 * k1);
        let k3 = f(t + 3.0 * h / 8.0, y + h * (3.0 * k1 + 9.0 * k2) / 32.0);
        let k4 = f(t + 12.0 * h / 13.0, y + h * (1932.0 * k1 - 7200.0 * k2 + 7296.0 * k3) / 2197.0);
        let k5 = f(t + h, y + h * (439.0 * k1 / 216.0 - 8.0 * k2 + 3680.0 * k3 / 513.0 - 845.0 * k4 / 4104.0));
        let k6 = f(
            t + h / 2.0,
            y + h * (-8.0 * k1 / 27.0 + 2.0 * k2 - 3544.0 * k3 / 2565.0 + 1859.0 * k4 / 4104.0 - 11.0 * k5 / 40.0),
        );

        // 4th-order solution
        let y4 = y + h * (25.0 * k1 / 216.0 + 1408.0 * k3 / 2565.0 + 2197.0 * k4 / 4104.0 - k5 / 5.0);
        // 5th-order solution
        let y5 = y + h * (16.0 * k1 / 135.0 + 6656.0 * k3 / 12825.0 + 28561.0 * k4 / 56430.0 - 9.0 * k5 / 50.0 + 2.0 * k6 / 55.0);

        let error = (y5 - y4).abs();
        if error < tol {
            y = y5;
            t += h;
        }

        // step-size adjustment
        if error > 0.0 {
            let s = (tol / (2.0 * error)).powf(0.25).clamp(0.1, 4.0);
            h *= s;
        }
        if h.abs() < 1e-15 {
            return Err(MathError::NotConvergent("rkf45: step size underflow".into()));
        }
    }
    Ok(y)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn euler_exponential() {
        // dy/dt = y, y(0) = 1 -> y(t) = e^t
        let f = |_t: f64, y: f64| y;
        let result = euler(f, 0.0, 1.0, 1.0, 10000).unwrap();
        assert!(close(result, std::f64::consts::E, 1e-3));
    }

    #[test]
    fn rk4_exponential() {
        let f = |_t: f64, y: f64| y;
        let result = rk4(f, 0.0, 1.0, 1.0, 100).unwrap();
        assert!(close(result, std::f64::consts::E, 1e-8));
    }

    #[test]
    fn rk4_sine() {
        // dy/dt = cos(t), y(0) = 0 -> y(t) = sin(t)
        let f = |t: f64, _y: f64| t.cos();
        let result = rk4(f, 0.0, std::f64::consts::PI, 0.0, 100).unwrap();
        assert!(close(result, 0.0, 1e-10));
    }

    #[test]
    fn rk4_trajectory_length() {
        let f = |_t: f64, y: f64| y;
        let traj = rk4_trajectory(f, 0.0, 1.0, 1.0, 10).unwrap();
        assert_eq!(traj.len(), 11);
        assert!(close(traj[0].0, 0.0, 1e-12));
        assert!(close(traj[10].0, 1.0, 1e-12));
    }

    #[test]
    fn rk4_system_harmonic_oscillator() {
        // y0 = position, y1 = velocity
        // dy0/dt = y1, dy1/dt = -y0  (omega = 1)
        // y(0) = [1, 0], solution: y0(t) = cos(t), y1(t) = -sin(t)
        let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
        let result = rk4_system(f, 0.0, std::f64::consts::PI / 2.0, &[1.0, 0.0], 100).unwrap();
        assert!(close(result[0], 0.0, 1e-8));
        assert!(close(result[1], -1.0, 1e-8));
    }

    #[test]
    fn rkf45_exponential() {
        let f = |_t: f64, y: f64| y;
        let result = rkf45(f, 0.0, 1.0, 1.0, 1e-8).unwrap();
        assert!(close(result, std::f64::consts::E, 1e-6));
    }

    #[test]
    fn rk4_vs_euler_accuracy() {
        // RK4 should be much more accurate than Euler for the same step count
        let f = |_t: f64, y: f64| y;
        let euler_result = euler(f, 0.0, 1.0, 1.0, 100).unwrap();
        let rk4_result = rk4(f, 0.0, 1.0, 1.0, 100).unwrap();
        let exact = std::f64::consts::E;
        let euler_err = (euler_result - exact).abs();
        let rk4_err = (rk4_result - exact).abs();
        assert!(rk4_err < euler_err * 1e-4);
    }
}
