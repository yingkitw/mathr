use mathr::interpolate;

fn main() {
    // --- Data points from f(x) = x^2 ---
    let points: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)];

    println!("Interpolation through points: {:?}", points);
    println!("  (underlying function: f(x) = x²)\n");

    // --- Lagrange interpolation ---
    let xs = [0.5, 1.5, 2.5];
    for &x in &xs {
        let y = interpolate::lagrange_interp(&points, x).unwrap();
        let exact = x * x;
        println!("  Lagrange  f({:.1}) = {:10.6}  (exact: {:10.6}, error: {:.2e})", x, y, exact, (y - exact).abs());
    }

    // --- Newton interpolation ---
    let newton = interpolate::NewtonInterpolator::new(&points).unwrap();
    for &x in &xs {
        let y = newton.eval(x);
        let exact = x * x;
        println!("  Newton    f({:.1}) = {:10.6}  (exact: {:10.6}, error: {:.2e})", x, y, exact, (y - exact).abs());
    }

    // --- Linear interpolation ---
    let lerp_y = interpolate::lerp(0.0, 10.0, 2.0, 20.0, 0.6);
    println!("\n  Linear lerp((0,10)→(2,20), x=0.6) = {:.4}", lerp_y);

    // --- Higher-degree example: f(x) = sin(x) ---
    let sin_points: Vec<(f64, f64)> = (0..=6)
        .map(|i| {
            let x = i as f64 * std::f64::consts::PI / 6.0;
            (x, x.sin())
        })
        .collect();
    println!("\nInterpolation of sin(x) through 7 equally-spaced points:");
    let newton_sin = interpolate::NewtonInterpolator::new(&sin_points).unwrap();
    for x in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5] {
        let y = newton_sin.eval(x);
        let exact = x.sin();
        println!("  Newton sin({:.2}) = {:10.6}  (exact: {:10.6}, error: {:.2e})", x, y, exact, (y - exact).abs());
    }
}
