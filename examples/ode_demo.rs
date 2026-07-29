use mathr::ode;

fn main() {
    // --- Euler vs RK4: dy/dt = y, y(0) = 1 (solution: e^t) ---
    let f = |_t: f64, y: f64| y;
    let t_end = 1.0;
    let y_euler = ode::euler(f, 0.0, t_end, 1.0, 100).unwrap();
    let y_rk4 = ode::rk4(f, 0.0, t_end, 1.0, 100).unwrap();
    let exact = std::f64::consts::E;
    println!("dy/dt = y, y(0) = 1, t_end = 1 (exact: e = {:.10})", exact);
    println!("  Euler (100 steps): {:.10}  error: {:.2e}", y_euler, (y_euler - exact).abs());
    println!("  RK4   (100 steps): {:.10}  error: {:.2e}", y_rk4, (y_rk4 - exact).abs());

    // --- RK4 for dy/dt = cos(t), y(0) = 0 (solution: sin(t)) ---
    let f2 = |t: f64, _y: f64| t.cos();
    let y_sin = ode::rk4(f2, 0.0, std::f64::consts::PI, 0.0, 200).unwrap();
    println!("\ndy/dt = cos(t), y(0) = 0, t_end = π (exact: 0)");
    println!("  RK4: {:.2e}", y_sin.abs());

    // --- RK4 system: harmonic oscillator ---
    // dy0/dt = y1, dy1/dt = -y0  →  y0(t) = cos(t), y1(t) = -sin(t)
    let sys = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
    let result = ode::rk4_system(sys, 0.0, std::f64::consts::PI / 2.0, &[1.0, 0.0], 100).unwrap();
    println!("\nHarmonic oscillator at t = π/2:");
    println!("  y0 = {:.10}  (exact: 0)", result[0]);
    println!("  y1 = {:.10}  (exact: -1)", result[1]);

    // --- Adaptive RKF45: dy/dt = -2y, y(0) = 1 (solution: e^{-2t}) ---
    let f3 = |_t: f64, y: f64| -2.0 * y;
    let y_rkf45 = ode::rkf45(f3, 0.0, 1.0, 1.0, 1e-8).unwrap();
    let exact3 = (-2.0_f64).exp();
    println!("\ndy/dt = -2y, y(0) = 1, t_end = 1 (exact: e⁻² = {:.10})", exact3);
    println!("  RKF45: {:.10}  error: {:.2e}", y_rkf45, (y_rkf45 - exact3).abs());
}
