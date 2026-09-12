//! Stiff ODE solvers: backward Euler, implicit trapezoidal, BDF2.
//!
//! Run with: cargo run --example stiff_demo

use mathr::stiff::{bdf2_system, implicit_euler_system, trapezoid_system};

fn main() {
    // 1. Stiff decay y' = -1000 (y - 1): explicit Euler needs h < 0.002.
    //    With h = 0.01 the implicit methods are stable, explicit is not.
    let f = |_t: f64, y: &[f64]| vec![-1000.0 * (y[0] - 1.0)];
    let be = implicit_euler_system(&f, 0.0, 1.0, &[0.0], 100).unwrap()[0];
    let trap = trapezoid_system(&f, 0.0, 1.0, &[0.0], 100).unwrap()[0];
    let bdf2 = bdf2_system(&f, 0.0, 1.0, &[0.0], 100).unwrap()[0];
    let expl = mathr::ode::euler(|_t, y| -1000.0 * (y - 1.0), 0.0, 1.0, 0.0, 100).unwrap();
    println!("y' = -1000 (y - 1), h = 0.01, t = 1  (exact ~ 1):");
    println!("  backward Euler : {be:.10}");
    println!("  trapezoidal    : {trap:.10}");
    println!("  BDF2           : {bdf2:.10}");
    println!("  explicit Euler : {expl:.3e}  (blown up)\n");

    // 2. Order check on y' = y: trapezoidal and BDF2 are 2nd order,
    //    backward Euler 1st
    let g = |_t: f64, y: &[f64]| vec![y[0]];
    let err = |v: f64| (v - std::f64::consts::E).abs();
    let be1 = err(implicit_euler_system(&g, 0.0, 1.0, &[1.0], 100).unwrap()[0]);
    let be2 = err(implicit_euler_system(&g, 0.0, 1.0, &[1.0], 200).unwrap()[0]);
    let tr1 = err(trapezoid_system(&g, 0.0, 1.0, &[1.0], 100).unwrap()[0]);
    let tr2 = err(trapezoid_system(&g, 0.0, 1.0, &[1.0], 200).unwrap()[0]);
    let b1 = err(bdf2_system(&g, 0.0, 1.0, &[1.0], 100).unwrap()[0]);
    let b2 = err(bdf2_system(&g, 0.0, 1.0, &[1.0], 200).unwrap()[0]);
    println!("y' = y error ratios when h halves (expect ~2, ~4, ~4):");
    println!("  backward Euler : {:.2}", be1 / be2);
    println!("  trapezoidal    : {:.2}", tr1 / tr2);
    println!("  BDF2           : {:.2}\n", b1 / b2);

    // 3. Robertson problem — classic stiff nonlinear system; the mass
    //    invariant y1 + y2 + y3 = 1 must hold through the stiff transient.
    let rob = |_t: f64, y: &[f64]| {
        vec![
            -0.04 * y[0] + 1.0e4 * y[1] * y[2],
            0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1],
            3.0e7 * y[1] * y[1],
        ]
    };
    for (name, y) in [
        ("backward Euler", implicit_euler_system(&rob, 0.0, 10.0, &[1.0, 0.0, 0.0], 1000).unwrap()),
        ("trapezoidal   ", trapezoid_system(&rob, 0.0, 10.0, &[1.0, 0.0, 0.0], 1000).unwrap()),
        ("BDF2          ", bdf2_system(&rob, 0.0, 10.0, &[1.0, 0.0, 0.0], 1000).unwrap()),
    ] {
        let total: f64 = y.iter().sum();
        println!("Robertson at t = 10 (h = 0.01), {name}:");
        println!("  y = [{:.6}, {:.3e}, {:.6}]  mass = {:.10}", y[0], y[1], y[2], total);
    }
}
