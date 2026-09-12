//! B-spline and cubic Hermite interpolation.
//!
//! Run with: cargo run --example bspline_demo

use mathr::bspline::cubic_bspline_interp;
use mathr::interpolate::CubicHermite;

fn main() {
    // 1. Cubic B-spline interpolation reproduces cubics exactly.
    let xs: Vec<f64> = (0..4).map(|i| i as f64).collect();
    let ys: Vec<f64> = xs.iter().map(|x| x * x * x + 2.0 * x).collect();
    let s = cubic_bspline_interp(&xs, &ys).unwrap();
    println!("y = x³ + 2x interpolated at x = 0, 1, 2, 3 (degree {}):", s.degree());
    for k in 1..6 {
        let x = k as f64 * 0.5;
        println!("  eval({x}) = {:.12}  (exact {})", s.eval(x), x * x * x + 2.0 * x);
    }

    // 2. Smooth interpolation of sin at 8 points.
    let xs: Vec<f64> = (0..8).map(|i| i as f64 * 0.5).collect();
    let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
    let s = cubic_bspline_interp(&xs, &ys).unwrap();
    let max_err: f64 = (0..64)
        .map(|k| {
            let x = k as f64 / 20.0;
            (s.eval(x) - x.sin()).abs()
        })
        .fold(0.0, f64::max);
    println!("\nsin interpolated at 8 points: max error on [0, 3.15] = {max_err:.2e}");

    // 3. Cubic Hermite with exact quadratic slopes reproduces y = x².
    let h = CubicHermite::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 2.0]).unwrap();
    println!("\ny = x² with slopes [0, 2]: hermite(0.25) = {:.12} (exact 0.0625)", h.eval(0.25));

    // 4. Same data, different slopes: slope choice shapes the curve.
    let h2 = CubicHermite::new(&[0.0, 1.0], &[0.0, 1.0], &[0.0, 0.0]).unwrap();
    println!("y endpoints [0, 1] with zero slopes: hermite(0.25) = {:.12}", h2.eval(0.25));
}
