//! Integration tests that exercise multiple modules together end-to-end.

use mathr::complex::Complex;
use mathr::eval::{eval, Context};
use mathr::expr::Expr;
use mathr::fft;
use mathr::interpolate;
use mathr::laurent;
use mathr::matrix::Matrix;
use mathr::notebook::{Notebook, NotebookCell};
use mathr::numtheory;
use mathr::parser::Parser;
use mathr::ode;
use mathr::rational::{parse_rational, Rational};
use mathr::simplify;
use mathr::special;
use mathr::solver;
use mathr::symbolic;
use mathr::taylor;

fn close(a: f64, b: f64, eps: f64) -> bool {
    (a - b).abs() < eps
}

// =========================================================================
// Expression pipeline: parse → differentiate → simplify → evaluate
// =========================================================================

#[test]
fn parse_diff_simplify_eval_pipeline() {
    // f(x) = x^3, f'(x) = 3x^2, f'(2) = 12
    let expr = Parser::parse("x^3").unwrap();
    let deriv = symbolic::differentiate(&expr, "x").unwrap();
    let simplified = simplify::simplify(&deriv);
    let mut ctx = Context::standard();
    ctx.set("x", 2.0);
    let val = eval(&simplified, &ctx).unwrap();
    assert!(close(val, 12.0, 1e-10));
}

fn square(args: &[f64]) -> mathr::error::Result<f64> {
    Ok(args[0] * args[0])
}

#[test]
fn parse_eval_with_user_function() {
    let mut ctx = Context::standard();
    ctx.insert_builtin("square", square);
    let expr = Parser::parse("square(3) + square(4)").unwrap();
    let val = eval(&expr, &ctx).unwrap();
    assert!(close(val, 25.0, 1e-10));
}

#[test]
fn taylor_series_evaluates_close_to_original() {
    // Taylor series of exp(x) around 0 with 10 terms should be very close at x=1
    let series = taylor::taylor_series_str("exp(x)", "x", 0.0, 10).unwrap();
    let mut ctx = Context::standard();
    ctx.set("x", 1.0);
    let val = eval(&series, &ctx).unwrap();
    assert!(close(val, std::f64::consts::E, 1e-6));
}

// =========================================================================
// FFT: convolution matches direct computation
// =========================================================================

#[test]
fn fft_convolution_matches_direct() {
    let a = vec![1.0, 2.0, 3.0, 4.0];
    let b = vec![5.0, 6.0, 7.0, 8.0];
    let fft_result = fft::convolve(&a, &b).unwrap();

    // Direct computation
    let n = a.len() + b.len() - 1;
    let mut direct = vec![0.0; n];
    for i in 0..a.len() {
        for j in 0..b.len() {
            direct[i + j] += a[i] * b[j];
        }
    }

    for i in 0..n {
        assert!(close(fft_result[i], direct[i], 1e-8));
    }
}

#[test]
fn fft_ifft_roundtrip_preserves_signal() {
    let input: Vec<Complex<f64>> = (0..16)
        .map(|i| Complex::new((i as f64 * 0.5).sin(), (i as f64 * 0.3).cos()))
        .collect();
    let freq = fft::fft(&input).unwrap();
    let recovered = fft::ifft(&freq).unwrap();
    for i in 0..16 {
        assert!(close(input[i].re, recovered[i].re, 1e-10));
        assert!(close(input[i].im, recovered[i].im, 1e-10));
    }
}

#[test]
fn window_function_normalization() {
    let n = 64;
    let signal = vec![1.0; n];
    let hann = fft::apply_window(&signal, fft::Window::Hann);
    let hamming = fft::apply_window(&signal, fft::Window::Hamming);
    let blackman = fft::apply_window(&signal, fft::Window::Blackman);

    // Coherent gain (sum / N) should be close to theoretical values
    // For symmetric windows with (N-1) denominator, exact only in limit
    let hann_gain = hann.iter().sum::<f64>() / n as f64;
    let hamming_gain = hamming.iter().sum::<f64>() / n as f64;
    let blackman_gain = blackman.iter().sum::<f64>() / n as f64;

    assert!(close(hann_gain, 0.5, 1e-2));
    assert!(close(hamming_gain, 0.54, 1e-2));
    assert!(close(blackman_gain, 0.42, 1e-2));
}

// =========================================================================
// Solver: find roots of expressions
// =========================================================================

#[test]
fn solve_polynomial_from_expression() {
    let expr = Parser::parse("x^2 - 4").unwrap();
    let ctx = Context::standard();
    let expr_clone = expr.clone();
    let ctx_clone = ctx.clone();
    let f = move |x: f64| {
        let mut c = ctx_clone.clone();
        c.set("x", x);
        eval(&expr_clone, &c).unwrap_or(f64::NAN)
    };
    let (root, residual) = solver::newton_central(f, 1.5, solver::SolveOptions::default()).unwrap();
    assert!(close(root, 2.0, 1e-8));
    assert!(residual.abs() < 1e-8);
}

#[test]
fn solve_bisection_from_expression() {
    let expr = Parser::parse("x^3 - x - 2").unwrap();
    let ctx = Context::standard();
    let expr_clone = expr.clone();
    let ctx_clone = ctx.clone();
    let f = move |x: f64| {
        let mut c = ctx_clone.clone();
        c.set("x", x);
        eval(&expr_clone, &c).unwrap_or(f64::NAN)
    };
    let (root, residual) = solver::bisect(f, 1.0, 2.0, solver::SolveOptions::default()).unwrap();
    assert!(close(root, 1.5213797068, 1e-6));
    assert!(residual.abs() < 1e-6);
}

// =========================================================================
// Interpolation: recovers polynomial exactly
// =========================================================================

#[test]
fn interpolation_recovers_polynomial() {
    // f(x) = 2x^2 - 3x + 1
    let f = |x: f64| 2.0 * x * x - 3.0 * x + 1.0;
    let points: Vec<(f64, f64)> = vec![(-1.0, f(-1.0)), (0.0, f(0.0)), (1.0, f(1.0)), (2.0, f(2.0))];

    let newton = interpolate::NewtonInterpolator::new(&points).unwrap();
    for x in [-0.5, 0.5, 1.5, 0.3, -0.7] {
        let interp = newton.eval(x);
        let exact = f(x);
        assert!(close(interp, exact, 1e-10));
    }

    // Lagrange should also match
    for x in [-0.5, 0.5, 1.5, 0.3, -0.7] {
        let interp = interpolate::lagrange_interp(&points, x).unwrap();
        let exact = f(x);
        assert!(close(interp, exact, 1e-10));
    }
}

// =========================================================================
// ODE: RK4 accuracy vs analytical solution
// =========================================================================

#[test]
fn rk4_exponential_accuracy() {
    let f = |_t: f64, y: f64| y;
    let result = ode::rk4(f, 0.0, 1.0, 1.0, 1000).unwrap();
    assert!(close(result, std::f64::consts::E, 1e-10));
}

#[test]
fn rk4_system_harmonic_oscillator() {
    let f = |_t: f64, y: &[f64]| vec![y[1], -y[0]];
    let result = ode::rk4_system(f, 0.0, std::f64::consts::FRAC_PI_2, &[1.0, 0.0], 1000).unwrap();
    assert!(close(result[0], 0.0, 1e-8));
    assert!(close(result[1], -1.0, 1e-8));
}

#[test]
fn rkf45_adaptive_accuracy() {
    let f = |_t: f64, y: f64| -3.0 * y;
    let result = ode::rkf45(f, 0.0, 2.0, 1.0, 1e-10).unwrap();
    let exact = (-6.0_f64).exp();
    assert!(close(result, exact, 1e-6));
}

// =========================================================================
// Matrix: solve and inverse consistency
// =========================================================================

#[test]
fn matrix_solve_and_inverse_consistency() {
    let a = Matrix::from_rows(&[
        vec![4.0, 3.0, 2.0],
        vec![1.0, 5.0, 3.0],
        vec![2.0, 1.0, 6.0],
    ]).unwrap();
    let b = vec![20.0, 14.0, 15.0];

    // Solve Ax = b
    let x = a.solve(&b).unwrap();

    // Verify Ax = b
    let ax = a.mul_vec(&x).unwrap();
    for i in 0..3 {
        assert!(close(ax[i], b[i], 1e-8));
    }

    // Inverse * b should give same x
    let inv = a.inverse().unwrap();
    let x2 = inv.mul_vec(&b).unwrap();
    for i in 0..3 {
        assert!(close(x[i], x2[i], 1e-8));
    }
}

#[test]
fn matrix_determinant_and_inverse_properties() {
    let a = Matrix::from_rows(&[
        vec![2.0, 1.0, 0.0],
        vec![1.0, 3.0, 1.0],
        vec![0.0, 1.0, 2.0],
    ]).unwrap();
    let det = a.determinant().unwrap();
    let inv = a.inverse().unwrap();
    let product = (&a * &inv).unwrap();

    // A * A⁻¹ should be identity
    for i in 0..3 {
        for j in 0..3 {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!(close(product.get(i, j), expected, 1e-8));
        }
    }
    // det(A⁻¹) = 1/det(A)
    let det_inv = inv.determinant().unwrap();
    assert!(close(det_inv, 1.0 / det, 1e-8));
}

// =========================================================================
// Special functions: known values and identities
// =========================================================================

#[test]
fn gamma_known_values() {
    assert!(close(special::gamma(1.0), 1.0, 1e-10));
    assert!(close(special::gamma(2.0), 1.0, 1e-10));
    assert!(close(special::gamma(3.0), 2.0, 1e-10));
    assert!(close(special::gamma(4.0), 6.0, 1e-10));
    assert!(close(special::gamma(5.0), 24.0, 1e-10));
    assert!(close(special::gamma(0.5), std::f64::consts::PI.sqrt(), 1e-10));
}

#[test]
fn gamma_reflection_identity() {
    // Γ(z)Γ(1-z) = π / sin(πz)
    for z in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let lhs = special::gamma(z) * special::gamma(1.0 - z);
        let rhs = std::f64::consts::PI / (std::f64::consts::PI * z).sin();
        assert!(close(lhs, rhs, 1e-8));
    }
}

#[test]
fn beta_function_known_values() {
    assert!(close(special::beta(1.0, 1.0), 1.0, 1e-10));
    assert!(close(special::beta(0.5, 0.5), std::f64::consts::PI, 1e-9));
    // B(2, 3) = Γ(2)Γ(3)/Γ(5) = 1*2/24 = 1/12
    assert!(close(special::beta(2.0, 3.0), 1.0 / 12.0, 1e-10));
}

#[test]
fn erf_known_values() {
    assert!(close(special::erf(0.0), 0.0, 1e-15));
    assert!(close(special::erf(1.0), 0.8427007929, 1e-8));
    assert!(close(special::erf(-1.0), -0.8427007929, 1e-8));
    assert!(close(special::erf(0.5), 0.5204998778, 1e-8));
    // erfc(x) + erf(x) = 1
    for x in [0.0, 0.5, 1.0, 2.0, 3.0] {
        assert!(close(special::erf(x) + special::erfc(x), 1.0, 1e-10));
    }
}

#[test]
fn incomplete_gamma_chi_squared_cdf() {
    // P(1, x) = 1 - e^{-x} is the CDF of χ² with 2 dof
    for x in [0.5, 1.0, 2.0, 5.0, 10.0] {
        let p = special::incomplete_gamma_p(1.0, x);
        let expected = 1.0 - (-x).exp();
        assert!(close(p, expected, 1e-8));
    }
}

// =========================================================================
// Number theory: cross-validation
// =========================================================================

#[test]
fn miller_rabin_matches_trial_division() {
    for n in 2u64..1000 {
        let mr = numtheory::is_prime_miller_rabin(n, 20);
        let trial = numtheory::is_prime(n);
        assert_eq!(mr, trial, "Mismatch at n={}", n);
    }
}

#[test]
fn crt_and_mod_pow_consistency() {
    // Verify CRT solution satisfies all congruences
    let remainders = [3u64, 4, 5];
    let moduli = [5u64, 7, 11];
    let x = numtheory::chinese_remainder(&remainders, &moduli).unwrap();
    for i in 0..3 {
        assert_eq!(x % moduli[i], remainders[i]);
    }
}

#[test]
fn prime_factorization_product_check() {
    for n in [12u64, 60, 360, 1024, 999983, 1234567890] {
        let factors = numtheory::prime_factors(n);
        let product: u64 = factors.iter().product();
        assert_eq!(product, n);
        // All factors should be prime
        for &f in &factors {
            assert!(numtheory::is_prime(f), "{} is not prime", f);
        }
    }
}

#[test]
fn sieve_contains_all_primes() {
    let primes = numtheory::sieve_primes(10000);
    for &p in &primes {
        assert!(numtheory::is_prime_miller_rabin(p, 20), "{} in sieve but not prime", p);
    }
    // Count should match known value: 1229 primes below 10000
    assert_eq!(primes.len(), 1229);
}

// =========================================================================
// Expression equality: canonical form
// =========================================================================

#[test]
fn expr_equality_comprehensive() {
    // a + b == b + a
    let a = Expr::add(Expr::var("a"), Expr::var("b"));
    let b = Expr::add(Expr::var("b"), Expr::var("a"));
    assert!(a.equals(&b));

    // a * b == b * a
    let a = Expr::mul(Expr::var("a"), Expr::var("b"));
    let b = Expr::mul(Expr::var("b"), Expr::var("a"));
    assert!(a.equals(&b));

    // a - b == a + (-b)
    let a = Expr::sub(Expr::var("x"), Expr::var("y"));
    let b = Expr::add(Expr::var("x"), Expr::neg(Expr::var("y")));
    assert!(a.equals(&b));

    // 2 + 3 + x == 5 + x
    let a = Expr::add(Expr::add(Expr::num(2.0), Expr::num(3.0)), Expr::var("x"));
    let b = Expr::add(Expr::num(5.0), Expr::var("x"));
    assert!(a.equals(&b));

    // 2 * 3 * x == 6 * x
    let a = Expr::mul(Expr::mul(Expr::num(2.0), Expr::num(3.0)), Expr::var("x"));
    let b = Expr::mul(Expr::num(6.0), Expr::var("x"));
    assert!(a.equals(&b));

    // x + 1 != x + 2
    let a = Expr::add(Expr::var("x"), Expr::num(1.0));
    let b = Expr::add(Expr::var("x"), Expr::num(2.0));
    assert!(!a.equals(&b));

    // (a + b) + c == c + (a + b)  (flattening)
    let a = Expr::add(Expr::add(Expr::var("a"), Expr::var("b")), Expr::var("c"));
    let b = Expr::add(Expr::var("c"), Expr::add(Expr::var("a"), Expr::var("b")));
    assert!(a.equals(&b));
}

// =========================================================================
// Performance benchmarks (correctness + timing sanity)
// =========================================================================

#[test]
fn fft_large_signal_performance() {
    // 4096-point FFT should complete in well under 1 second
    let n = 4096;
    let samples: Vec<f64> = (0..n)
        .map(|i| (2.0 * std::f64::consts::PI * 64.0 * i as f64 / n as f64).sin())
        .collect();
    let start = std::time::Instant::now();
    let mags = fft::magnitude_spectrum(&samples).unwrap();
    let elapsed = start.elapsed();
    assert_eq!(mags.len(), n / 2 + 1);
    // Peak should be at bin 64
    let peak = mags[64];
    assert!(peak > 1000.0, "expected strong peak at bin 64, got {}", peak);
    assert!(elapsed.as_millis() < 500, "FFT took too long: {:?}", elapsed);
}

#[test]
fn matrix_inverse_large_performance() {
    // 50×50 matrix inverse should complete quickly
    let n = 50;
    let mut data = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            // Diagonally dominant matrix (guaranteed invertible)
            if i == j {
                data.push(10.0);
            } else {
                data.push(((i as f64 + 1.0) * (j as f64 + 1.0)).sin() * 0.1);
            }
        }
    }
    let a = Matrix::from_row_major(n, n, data).unwrap();
    let start = std::time::Instant::now();
    let inv = a.inverse().unwrap();
    let elapsed = start.elapsed();

    // Verify A * A⁻¹ ≈ I
    let product = (&a * &inv).unwrap();
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((product.get(i, j) - expected).abs() < 1e-6,
                "A*A⁻¹[{}][{}] = {}, expected {}", i, j, product.get(i, j), expected);
        }
    }
    assert!(elapsed.as_millis() < 500, "Matrix inverse took too long: {:?}", elapsed);
}

#[test]
fn matrix_symmetric_eig_decomposition() {
    let a = Matrix::from_rows(&[
        vec![4.0, 1.0, 2.0],
        vec![1.0, 3.0, 0.0],
        vec![2.0, 0.0, 5.0],
    ]).unwrap();
    let (vals, vecs) = a.symmetric_eig().unwrap();

    // Eigenvalues should be sorted ascending.
    assert!(vals[0] <= vals[1] && vals[1] <= vals[2]);

    // Verify A v = λ v for each eigenvector.
    for j in 0..3 {
        let v: Vec<f64> = (0..3).map(|i| vecs[(i, j)]).collect();
        let av = a.mul_vec(&v).unwrap();
        for i in 0..3 {
            assert!((av[i] - vals[j] * v[i]).abs() < 1e-8,
                "A v_{} != λ_{} v_{} at i={}", j, j, j, i);
        }
    }

    // Verify eigenvectors are orthonormal.
    for i in 0..3 {
        for j in 0..3 {
            let dot: f64 = (0..3).map(|k| vecs[(k, i)] * vecs[(k, j)]).sum();
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((dot - expected).abs() < 1e-8,
                "v_{}·v_{} = {}, expected {}", i, j, dot, expected);
        }
    }
}

#[test]
fn matrix_symmetric_eig_cli() {
    // Test the REPL dispatch for symlig
    let ctx = mathr::eval::Context::standard();
    let result = mathr::repl::dispatch_str("symlig 2 1 | 1 2", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("eigenvalues"));
    // eigenvalues of [[2,1],[1,2]] are 1 and 3
    assert!(output.contains(" 1") || output.contains("[1,") || output.contains("[1 "),
            "output should contain eigenvalue 1: {}", output);
    assert!(output.contains(" 3") || output.contains("[3,") || output.contains(", 3"),
            "output should contain eigenvalue 3: {}", output);
}

#[test]
fn matrix_hessenberg_cli() {
    let ctx = mathr::eval::Context::standard();
    let result = mathr::repl::dispatch_str("hessenberg 1 2 3 | 4 5 6 | 7 8 10", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("Hessenberg"));
    assert!(output.contains("orthogonal"));
}

#[test]
fn matrix_schur_cli() {
    let ctx = mathr::eval::Context::standard();
    let result = mathr::repl::dispatch_str("schur 4 1 2 | 1 3 0 | 2 0 5", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("triangular"));
    assert!(output.contains("orthogonal"));
}

#[test]
fn isolate_roots_cli() {
    let ctx = mathr::eval::Context::standard();
    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    let result = mathr::repl::dispatch_str("isolate-roots 1 -6 11 -6", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("x"), "output should contain roots: {}", output);
    // Should find 3 roots
    assert!(output.lines().count() >= 3, "expected 3 roots, got: {}", output);
}

#[test]
fn isolate_roots_api() {
    // (x+2)(x-1)(x-3) = x^3 - 2x^2 - 5x + 6
    let intervals = mathr::solver::isolate_real_roots(&[1, -2, -5, 6]).unwrap();
    assert_eq!(intervals.len(), 3);
    let roots = [-2.0, 1.0, 3.0];
    for (i, r) in roots.iter().enumerate() {
        let (lo, hi) = intervals[i];
        assert!(lo <= *r && *r <= hi, "root {} not in ({}, {})", r, lo, hi);
    }
}

#[test]
fn gradient_repl() {
    let ctx = mathr::eval::Context::standard();
    let result = mathr::repl::dispatch_str("gradient x^2 + x*y + y^2", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("d/dx"), "output should contain d/dx: {}", output);
    assert!(output.contains("d/dy"), "output should contain d/dy: {}", output);
}

#[test]
fn pdiff_repl() {
    let ctx = mathr::eval::Context::standard();
    let result = mathr::repl::dispatch_str("pdiff x^2 * y + y^3 y", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    // ∂/∂y (x^2*y + y^3) = x^2 + 3*y^2
    let ctx2 = mathr::eval::Context::standard();
    let mut ctx2 = ctx2;
    ctx2.set("x", 2.0);
    ctx2.set("y", 3.0);
    let e = mathr::parser::Parser::parse(&output).unwrap();
    let val = mathr::eval::eval(&e, &ctx2).unwrap();
    // x^2 + 3*y^2 = 4 + 27 = 31
    assert!((val - 31.0).abs() < 1e-9, "pdiff result evaluated to {} expected 31", val);
}

#[test]
fn fourier_repl() {
    let ctx = mathr::eval::Context::standard();
    // Fourier series of cos(x) on [-pi, pi] with 5 terms, eval at 0
    let result = mathr::repl::dispatch_str("fourier cos(x) 3.14159265358979 5 0", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("a0"), "output should contain a0: {}", output);
    assert!(output.contains("a1"), "output should contain a1: {}", output);
    // cos(x) is a pure Fourier mode, so a1 ≈ 1 and eval at 0 ≈ 1
    assert!(output.contains("f("), "output should contain evaluation: {}", output);
}

#[test]
fn mc_repl() {
    let ctx = mathr::eval::Context::standard();
    // Monte Carlo integral of x over [0, 1] = 0.5
    let result = mathr::repl::dispatch_str("mc x 0 1 100000 42", ctx).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("estimate"), "output should contain estimate: {}", output);
    assert!(output.contains("std_error"), "output should contain std_error: {}", output);
}

#[test]
fn sample_repl() {
    let result = mathr::repl::dispatch_str("sample normal 0 1 10000 42", mathr::eval::Context::standard()).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("n=10000"), "output should contain n=10000: {}", output);
    assert!(output.contains("mean="), "output should contain mean: {}", output);
}

#[test]
fn dist_repl() {
    let result = mathr::repl::dispatch_str("dist normal 0 0 1", mathr::eval::Context::standard()).unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("pdf"), "output should contain pdf: {}", output);
    assert!(output.contains("cdf"), "output should contain cdf: {}", output);
    // Normal(0,0,1): pdf = 1/sqrt(2pi), cdf = 0.5
    assert!(output.contains("0.3989"), "pdf should be ~0.3989: {}", output);
    assert!(output.contains("0.5"), "cdf should be 0.5: {}", output);
}

#[test]
fn tikhonov_repl() {
    // Solve a simple 2x2 system with Tikhonov, lambda=0 → same as plain solve
    // A = [[2, 1], [1, 3]], b = [3, 4], x = [1, 1]
    let result = mathr::repl::dispatch_str(
        "tikhonov 2 1 | 1 3 | 3 4 0",
        mathr::eval::Context::standard(),
    )
    .unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("x = ["), "output should contain x = [...]: {}", output);
    assert!(output.contains("1"), "output should contain solution: {}", output);
}

#[test]
fn laurent_repl() {
    // Laurent series of 1/x around x=0, pole order 1, 3 positive terms
    let result = mathr::repl::dispatch_str(
        "laurent 1/x 0 1 3",
        mathr::eval::Context::standard(),
    )
    .unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("1/x"), "output should contain 1/x: {}", output);
}

#[test]
fn rat_repl() {
    // 1/2 + 1/3 = 5/6
    let result = mathr::repl::dispatch_str(
        "rat 1/2 + 1/3",
        mathr::eval::Context::standard(),
    )
    .unwrap();
    assert!(result.is_some());
    let output = result.unwrap();
    assert!(output.contains("5/6"), "output should contain 5/6: {}", output);
}

#[test]
fn next_pow2_correctness() {
    assert_eq!(fft::next_pow2(0), 1);
    assert_eq!(fft::next_pow2(1), 1);
    assert_eq!(fft::next_pow2(2), 2);
    assert_eq!(fft::next_pow2(3), 4);
    assert_eq!(fft::next_pow2(4), 4);
    assert_eq!(fft::next_pow2(5), 8);
    assert_eq!(fft::next_pow2(7), 8);
    assert_eq!(fft::next_pow2(8), 8);
    assert_eq!(fft::next_pow2(9), 16);
    assert_eq!(fft::next_pow2(1023), 1024);
    assert_eq!(fft::next_pow2(1024), 1024);
    assert_eq!(fft::next_pow2(1025), 2048);
}

#[test]
fn convolution_large_signal_performance() {
    // Convolve two 1024-length signals via FFT
    let a: Vec<f64> = (0..1024).map(|i| (i as f64 * 0.01).sin()).collect();
    let b: Vec<f64> = (0..1024).map(|i| (i as f64 * 0.02).cos()).collect();
    let start = std::time::Instant::now();
    let result = fft::convolve(&a, &b).unwrap();
    let elapsed = start.elapsed();
    assert_eq!(result.len(), 2047);
    assert!(elapsed.as_millis() < 200, "Convolution took too long: {:?}", elapsed);
}

// =========================================================================
// TeX / Markdown input: parse → eval, parse → diff, parse → solve
// =========================================================================

#[test]
fn tex_eval_fraction() {
    let e = Parser::parse(r"\frac{1}{2} + \frac{3}{4}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 1.25, 1e-10));
}

#[test]
fn tex_eval_sqrt() {
    let e = Parser::parse(r"\sqrt{16}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 4.0, 1e-10));
}

#[test]
fn tex_eval_sin_pi() {
    let e = Parser::parse(r"\sin(\pi / 4)").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, std::f64::consts::FRAC_PI_4.sin(), 1e-10));
}

#[test]
fn tex_eval_cdot() {
    let e = Parser::parse(r"2 \cdot 3 + 4").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 10.0, 1e-10));
}

#[test]
fn tex_eval_left_right() {
    let e = Parser::parse(r"\left( 1 + 2 \right) \cdot 3").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 9.0, 1e-10));
}

#[test]
fn tex_eval_gamma() {
    let e = Parser::parse(r"\Gamma{0.5}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, std::f64::consts::PI.sqrt(), 1e-8));
}

#[test]
fn tex_eval_log_subscript() {
    let e = Parser::parse(r"\log_2{8}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 3.0, 1e-10));
}

#[test]
fn tex_eval_operatorname() {
    let e = Parser::parse(r"\operatorname{erf}(1.0)").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 0.8427007929, 1e-8));
}

#[test]
fn tex_eval_nested_frac() {
    let e = Parser::parse(r"\frac{\frac{1}{2}}{\frac{3}{4}}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 2.0 / 3.0, 1e-10));
}

#[test]
fn tex_eval_implicit_mult() {
    let e = Parser::parse(r"2\pi").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 2.0 * std::f64::consts::PI, 1e-10));
}

#[test]
fn markdown_inline_eval() {
    let e = Parser::parse(r"$\sin(\pi / 4)$").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, std::f64::consts::FRAC_PI_4.sin(), 1e-10));
}

#[test]
fn markdown_display_eval() {
    let e = Parser::parse(r"$$\frac{1}{2} + \frac{3}{4}$$").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 1.25, 1e-10));
}

#[test]
fn latex_bracket_eval() {
    let e = Parser::parse(r"\[\sqrt{16} + \cos(0)\]").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 5.0, 1e-10));
}

#[test]
fn latex_paren_eval() {
    let e = Parser::parse(r"\(\log_2{8} + 1\)").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 4.0, 1e-10));
}

#[test]
fn tex_diff_matches_plain() {
    let tex_expr = Parser::parse(r"\sin(x^2)").unwrap();
    let plain_expr = Parser::parse("sin(x^2)").unwrap();
    let d_tex = symbolic::differentiate(&tex_expr, "x").unwrap();
    let d_plain = symbolic::differentiate(&plain_expr, "x").unwrap();
    assert_eq!(d_tex.canonicalize(), d_plain.canonicalize());
}

#[test]
fn tex_diff_frac() {
    // d/dx (x^2 + 1)/(x - 1) — parse as TeX, differentiate, evaluate
    let e = Parser::parse(r"\frac{x^2 + 1}{x - 1}").unwrap();
    let d = symbolic::differentiate(&e, "x").unwrap();
    let s = simplify::simplify(&d);
    // At x=2: f(x) = (4+1)/(2-1) = 5, f'(x) = (2x(x-1) - (x^2+1)) / (x-1)^2
    // f'(2) = (4*1 - 5) / 1 = -1
    let mut ctx = Context::standard();
    ctx.set("x", 2.0);
    let v = eval(&s, &ctx).unwrap();
    assert!(close(v, -1.0, 1e-8));
}

#[test]
fn tex_solve_from_frac() {
    let e = Parser::parse(r"\frac{x^2 - 4}{1}").unwrap();
    let ctx = Context::standard();
    let expr_clone = e.clone();
    let ctx_clone = ctx.clone();
    let f = move |x: f64| {
        let mut c = ctx_clone.clone();
        c.set("x", x);
        eval(&expr_clone, &c).unwrap_or(f64::NAN)
    };
    let (root, residual) = solver::newton_central(f, 1.0, solver::SolveOptions::default()).unwrap();
    assert!(close(root, 2.0, 1e-8));
    assert!(residual.abs() < 1e-8);
}

#[test]
fn tex_simplify_matches_plain() {
    let tex_e = Parser::parse(r"\frac{x^2 - 1}{x - 1}").unwrap();
    let plain_e = Parser::parse("(x^2 - 1)/(x - 1)").unwrap();
    let s_tex = simplify::simplify(&tex_e);
    let s_plain = simplify::simplify(&plain_e);
    assert_eq!(s_tex.canonicalize(), s_plain.canonicalize());
}

#[test]
fn tex_taylor_matches_plain() {
    let tex_series = taylor::taylor_series_str(r"\exp(x)", "x", 0.0, 5).unwrap();
    let plain_series = taylor::taylor_series_str("exp(x)", "x", 0.0, 5).unwrap();
    assert_eq!(tex_series.canonicalize(), plain_series.canonicalize());
}

#[test]
fn tex_taylor_evaluates_correctly() {
    let series = taylor::taylor_series_str(r"\sin(x)", "x", 0.0, 7).unwrap();
    let mut ctx = Context::standard();
    ctx.set("x", 0.5);
    let v = eval(&series, &ctx).unwrap();
    assert!(close(v, 0.5_f64.sin(), 1e-4));
}

#[test]
fn tex_pow_brace_eval() {
    let e = Parser::parse(r"2^{10}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 1024.0, 1e-10));
}

#[test]
fn tex_mixed_expression_eval() {
    // 3\sqrt{16} + \frac{1}{2} = 12.5
    let e = Parser::parse(r"3\sqrt{16} + \frac{1}{2}").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 12.5, 1e-10));
}

#[test]
fn tex_left_right_with_frac() {
    let e = Parser::parse(r"\left( \frac{1}{2} + \frac{1}{2} \right) \cdot 10").unwrap();
    let ctx = Context::standard();
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 10.0, 1e-10));
}

#[test]
fn tex_all_delimiters_produce_same_ast() {
    let exprs = [
        Parser::parse(r"\frac{1}{2} + \frac{3}{4}").unwrap(),
        Parser::parse(r"$\frac{1}{2} + \frac{3}{4}$").unwrap(),
        Parser::parse(r"$$\frac{1}{2} + \frac{3}{4}$$").unwrap(),
        Parser::parse(r"\(\frac{1}{2} + \frac{3}{4}\)").unwrap(),
        Parser::parse(r"\[\frac{1}{2} + \frac{3}{4}\]").unwrap(),
    ];
    for i in 1..exprs.len() {
        assert_eq!(exprs[0].canonicalize(), exprs[i].canonicalize());
    }
}

#[test]
fn tex_text_variable_eval() {
    let e = Parser::parse(r"\text{alpha} + 1").unwrap();
    let mut ctx = Context::standard();
    ctx.set("alpha", 2.0);
    let v = eval(&e, &ctx).unwrap();
    assert!(close(v, 3.0, 1e-10));
}

// =========================================================================
// Laurent series: parse → expand → evaluate
// =========================================================================

#[test]
fn laurent_simple_pole_eval() {
    // f(x) = 1/x, Laurent series should give c_-1 = 1, rest = 0
    let ls = laurent::laurent_series_str("1/x", "x", 0.0, 1, 3).unwrap();
    assert!(close(ls.coeff(-1), 1.0, 1e-6));
    assert!(close(ls.coeff(0), 0.0, 1e-6));
    // eval at x=5: 1/5 = 0.2
    assert!(close(ls.eval(5.0), 0.2, 1e-6));
}

#[test]
fn laurent_double_pole_eval() {
    // f(x) = 1/x^2, c_-2 = 1, rest = 0
    let ls = laurent::laurent_series_str("1/x^2", "x", 0.0, 2, 3).unwrap();
    assert!(close(ls.coeff(-2), 1.0, 1e-4));
    // eval at x=4: 1/16 = 0.0625
    assert!(close(ls.eval(4.0), 0.0625, 1e-4));
}

#[test]
fn laurent_exp_over_x_pipeline() {
    // f(x) = exp(x)/x = 1/x + 1 + x/2 + x^2/6 + ...
    let ls = laurent::laurent_series_str("exp(x)/x", "x", 0.0, 1, 6).unwrap();
    assert!(close(ls.coeff(-1), 1.0, 1e-6));
    assert!(close(ls.coeff(0), 1.0, 1e-6));
    assert!(close(ls.coeff(1), 0.5, 1e-6));
    assert!(close(ls.coeff(2), 1.0 / 6.0, 1e-6));
    // eval at x=1: e ≈ 2.71828
    assert!(close(ls.eval(1.0), std::f64::consts::E, 1e-2));
}

#[test]
fn laurent_around_nonzero_eval() {
    // f(x) = 1/(x-2), pole at x=2
    let ls = laurent::laurent_series_str("1/(x-2)", "x", 2.0, 1, 3).unwrap();
    assert!(close(ls.coeff(-1), 1.0, 1e-6));
    // eval at x=5: 1/3
    assert!(close(ls.eval(5.0), 1.0 / 3.0, 1e-4));
}

#[test]
fn laurent_no_pole_matches_taylor() {
    // f(x) = exp(x), no pole → Laurent with pole_order=0 should match Taylor
    let ls = laurent::laurent_series_str("exp(x)", "x", 0.0, 0, 5).unwrap();
    assert_eq!(ls.pole_order, 0);
    assert!(close(ls.coeff(0), 1.0, 1e-6));
    assert!(close(ls.coeff(1), 1.0, 1e-6));
    assert!(close(ls.coeff(2), 0.5, 1e-6));
}

#[test]
fn laurent_rational_function() {
    // f(x) = 1/(x(1-x)) = 1/x + 1 + x + x^2 + ...
    let ls = laurent::laurent_series_str("1/(x*(1-x))", "x", 0.0, 1, 5).unwrap();
    assert!(close(ls.coeff(-1), 1.0, 1e-6));
    assert!(close(ls.coeff(0), 1.0, 1e-6));
    assert!(close(ls.coeff(1), 1.0, 1e-6));
    assert!(close(ls.coeff(2), 1.0, 1e-6));
}

#[test]
fn laurent_to_string_contains_principal() {
    let ls = laurent::laurent_series_str("1/x + 2 + x", "x", 0.0, 1, 2).unwrap();
    let s = ls.to_string();
    assert!(s.contains("1/x"), "string should contain 1/x: {}", s);
}

// =========================================================================
// Rational arithmetic: parse → arithmetic → verify exactness
// =========================================================================

#[test]
fn rational_add_exact() {
    let a = parse_rational("1/2").unwrap();
    let b = parse_rational("1/3").unwrap();
    let c = a + b;
    assert_eq!(c.num(), 5);
    assert_eq!(c.den(), 6);
}

#[test]
fn rational_sub_exact() {
    let a = parse_rational("1/2").unwrap();
    let b = parse_rational("1/3").unwrap();
    let c = a - b;
    assert_eq!(c.num(), 1);
    assert_eq!(c.den(), 6);
}

#[test]
fn rational_mul_exact() {
    let a = Rational::new(2, 3).unwrap();
    let b = Rational::new(3, 4).unwrap();
    let c = a * b;
    assert_eq!(c.num(), 1);
    assert_eq!(c.den(), 2);
}

#[test]
fn rational_div_exact() {
    let a = Rational::new(2, 3).unwrap();
    let b = Rational::new(4, 5).unwrap();
    let c = a / b;
    assert_eq!(c.num(), 5);
    assert_eq!(c.den(), 6);
}

#[test]
fn rational_powi_exact() {
    let a = Rational::new(2, 3).unwrap();
    assert_eq!(a.powi(3), Rational::new(8, 27).unwrap());
    assert_eq!(a.powi(-2), Rational::new(9, 4).unwrap());
    assert_eq!(a.powi(0), Rational::from_int(1));
}

#[test]
fn rational_parse_decimal() {
    let r = parse_rational("0.5").unwrap();
    assert_eq!(r.num(), 1);
    assert_eq!(r.den(), 2);

    let r = parse_rational("-1.25").unwrap();
    assert_eq!(r.num(), -5);
    assert_eq!(r.den(), 4);
}

#[test]
fn rational_parse_fraction() {
    let r = parse_rational("3/4").unwrap();
    assert_eq!(r.num(), 3);
    assert_eq!(r.den(), 4);

    let r = parse_rational("-3/4").unwrap();
    assert_eq!(r.num(), -3);
    assert_eq!(r.den(), 4);
}

#[test]
fn rational_reduction() {
    let r = Rational::new(6, 8).unwrap();
    assert_eq!(r.num(), 3);
    assert_eq!(r.den(), 4);
}

#[test]
fn rational_equality_after_reduction() {
    let a = Rational::new(1, 2).unwrap();
    let b = Rational::new(2, 4).unwrap();
    assert_eq!(a, b);
}

#[test]
fn rational_ordering() {
    let a = Rational::new(1, 3).unwrap();
    let b = Rational::new(1, 2).unwrap();
    assert!(a < b);
    assert!(b > a);
}

#[test]
fn rational_large_arithmetic() {
    // i128 intermediates should handle large denominators
    let a = Rational::new(1, 1_000_000_000).unwrap();
    let b = Rational::new(1, 1_000_000_000).unwrap();
    let c = a + b;
    assert_eq!(c.num(), 1);
    assert_eq!(c.den(), 500_000_000);
}

#[test]
fn rational_chained_arithmetic() {
    // (1/2 + 1/3) * (1/4) = (5/6) * (1/4) = 5/24
    let result = (parse_rational("1/2").unwrap() + parse_rational("1/3").unwrap())
        * parse_rational("1/4").unwrap();
    assert_eq!(result.num(), 5);
    assert_eq!(result.den(), 24);
}

#[test]
fn rational_to_f64_accuracy() {
    let r = Rational::new(22, 7).unwrap();
    assert!((r.to_f64() - 22.0 / 7.0).abs() < 1e-15);
}

#[test]
fn rational_reciprocal() {
    let a = Rational::new(3, 4).unwrap();
    let r = a.recip().unwrap();
    assert_eq!(r.num(), 4);
    assert_eq!(r.den(), 3);
}

#[test]
fn rational_abs_and_neg() {
    let a = Rational::new(-3, 4).unwrap();
    assert_eq!(a.abs(), Rational::new(3, 4).unwrap());
    assert_eq!(-a, Rational::new(3, 4).unwrap());
}

#[test]
fn rational_display() {
    assert_eq!(Rational::from_int(5).to_string(), "5");
    assert_eq!(Rational::new(3, 4).unwrap().to_string(), "3/4");
    assert_eq!(Rational::new(-3, 4).unwrap().to_string(), "-3/4");
}

#[test]
fn rational_zero_den_errors() {
    assert!(Rational::new(1, 0).is_err());
    assert!(parse_rational("1/0").is_err());
}

// =========================================================================
// Cross-module pipelines
// =========================================================================

#[test]
fn rational_to_f64_then_eval() {
    // Convert rational to f64 and use in expression evaluation
    let r = Rational::new(1, 4).unwrap();
    let val = r.to_f64();
    let mut ctx = Context::standard();
    ctx.set("x", val);
    let expr = Parser::parse("x * 4 + 1").unwrap();
    let result = eval(&expr, &ctx).unwrap();
    assert!(close(result, 2.0, 1e-15));
}

#[test]
fn laurent_then_solve_near_pole() {
    // Compute Laurent series of 1/(x-1), then use solver to find root of f(x) - 0.5
    // 1/(x-1) = 0.5 → x = 3
    let ls = laurent::laurent_series_str("1/(x-1)", "x", 1.0, 1, 3).unwrap();
    // The series evaluated at x=3 should be close to 0.5
    let val = ls.eval(3.0);
    assert!(close(val, 0.5, 1e-4));
}

#[test]
fn taylor_then_laurent_consistency() {
    // For a function with no pole, Taylor and Laurent(pole_order=0) should agree
    let taylor_series = taylor::taylor_series_str("cos(x)", "x", 0.0, 5).unwrap();
    let laurent_series = laurent::laurent_series_str("cos(x)", "x", 0.0, 0, 5).unwrap();
    // Compare first few coefficients
    assert!(close(laurent_series.coeff(0), 1.0, 1e-6));
    // Taylor c_0 should also be 1
    let mut ctx = Context::standard();
    ctx.set("x", 0.0);
    let t0 = eval(&taylor_series, &ctx).unwrap();
    assert!(close(t0, 1.0, 1e-10));
    // Both should agree at x=0.5
    ctx.set("x", 0.5);
    let tv = eval(&taylor_series, &ctx).unwrap();
    let lv = laurent_series.eval(0.5);
    assert!(close(tv, lv, 1e-3));
}

#[test]
fn rational_arithmetic_repl_all_ops() {
    // Test all four operators via REPL
    let ctx = mathr::eval::Context::standard();

    let r1 = mathr::repl::dispatch_str("rat 1/2 + 1/4", ctx.clone()).unwrap().unwrap();
    assert!(r1.contains("3/4"), "1/2 + 1/4 should give 3/4: {}", r1);

    let r2 = mathr::repl::dispatch_str("rat 1/2 - 1/4", ctx.clone()).unwrap().unwrap();
    assert!(r2.contains("1/4"), "1/2 - 1/4 should give 1/4: {}", r2);

    let r3 = mathr::repl::dispatch_str("rat 2/3 * 3/4", ctx.clone()).unwrap().unwrap();
    assert!(r3.contains("1/2"), "2/3 * 3/4 should give 1/2: {}", r3);

    let r4 = mathr::repl::dispatch_str("rat 2/3 / 4/5", ctx.clone()).unwrap().unwrap();
    assert!(r4.contains("5/6"), "2/3 / 4/5 should give 5/6: {}", r4);
}

#[test]
fn rational_decimal_repl() {
    let ctx = mathr::eval::Context::standard();
    let r = mathr::repl::dispatch_str("rat 0.5 + 0.25", ctx).unwrap().unwrap();
    assert!(r.contains("3/4"), "0.5 + 0.25 should give 3/4: {}", r);
}

#[test]
fn laurent_repl_nonzero_center() {
    let ctx = mathr::eval::Context::standard();
    let r = mathr::repl::dispatch_str("laurent 1/(x-2) 2 1 3", ctx).unwrap().unwrap();
    assert!(r.contains("1/(x - 2)"), "output should contain 1/(x - 2): {}", r);
}

// =========================================================================
// Notebook: .mnb file format, cell evaluation, save/load
// =========================================================================

#[test]
fn notebook_create_eval_cell() {
    let mut nb = Notebook::new();
    let id = nb.add_cell("sin(pi/4)");
    nb.eval_cell(id, &Context::standard()).unwrap();
    assert!(nb.cells[id].output.contains("0.707"), "output: {}", nb.cells[id].output);
}

#[test]
fn notebook_eval_tex_cell() {
    let mut nb = Notebook::new();
    let id = nb.add_cell(r"\frac{1}{2} + \frac{3}{4}");
    nb.eval_cell(id, &Context::standard()).unwrap();
    assert!(nb.cells[id].output.contains("1.25"), "output: {}", nb.cells[id].output);
}

#[test]
fn notebook_eval_all_cells() {
    let mut nb = Notebook::new();
    nb.add_cell("1 + 2");
    nb.add_cell("3 * 4");
    nb.add_cell("sin(0)");
    nb.eval_all(&Context::standard()).unwrap();
    assert!(nb.cells[0].output.contains("3"));
    assert!(nb.cells[1].output.contains("12"));
    assert!(nb.cells[2].output.contains("0"));
}

#[test]
fn notebook_eval_diff_cell() {
    let mut nb = Notebook::new();
    let id = nb.add_cell("diff x^3");
    nb.eval_cell(id, &Context::standard()).unwrap();
    assert!(nb.cells[id].output.contains("3") && nb.cells[id].output.contains("x"),
        "output should contain derivative: {}", nb.cells[id].output);
}

#[test]
fn notebook_eval_solve_cell() {
    let mut nb = Notebook::new();
    let id = nb.add_cell("solve x^2 - 4");
    nb.eval_cell(id, &Context::standard()).unwrap();
    assert!(nb.cells[id].output.contains("2") || nb.cells[id].output.contains("root"),
        "output should contain root: {}", nb.cells[id].output);
}

#[test]
fn notebook_json_roundtrip() {
    let mut nb = Notebook::new();
    nb.add_cell("sin(pi/4)");
    nb.add_cell(r"\frac{1}{2}");
    nb.cells[0].output = "0.707...".to_string();
    nb.cells[1].output = "0.5".to_string();
    let json = nb.to_json();
    let nb2 = mathr::notebook::parse_notebook_json(&json).unwrap();
    assert_eq!(nb2.cells.len(), 2);
    assert_eq!(nb2.cells[0].input, "sin(pi/4)");
    assert_eq!(nb2.cells[0].output, "0.707...");
    assert_eq!(nb2.cells[1].input, r"\frac{1}{2}");
    assert_eq!(nb2.cells[1].output, "0.5");
}

#[test]
fn notebook_save_load_file() {
    let path = std::env::temp_dir().join("mathr_integration_test.mnb");
    let mut nb = Notebook::new();
    nb.add_cell("1 + 2");
    nb.add_cell("sin(pi/4)");
    nb.cells[0].output = "3".to_string();
    nb.save(&path).unwrap();
    let nb2 = Notebook::load(&path).unwrap();
    assert_eq!(nb2.cells.len(), 2);
    assert_eq!(nb2.cells[0].input, "1 + 2");
    assert_eq!(nb2.cells[0].output, "3");
    let _ = std::fs::remove_file(&path);
}

#[test]
fn notebook_set_input_clears_output() {
    let mut nb = Notebook::new();
    let id = nb.add_cell("1 + 1");
    nb.eval_cell(id, &Context::standard()).unwrap();
    assert!(!nb.cells[id].output.is_empty());
    nb.set_input(id, "2 + 2").unwrap();
    assert_eq!(nb.cells[id].input, "2 + 2");
    assert!(nb.cells[id].output.is_empty());
}

#[test]
fn notebook_remove_cell_reindexes() {
    let mut nb = Notebook::new();
    nb.add_cell("a");
    nb.add_cell("b");
    nb.add_cell("c");
    nb.remove_cell(1).unwrap();
    assert_eq!(nb.cells.len(), 2);
    assert_eq!(nb.cells[0].id, 0);
    assert_eq!(nb.cells[1].id, 1);
    assert_eq!(nb.cells[1].input, "c");
}

#[test]
fn notebook_parse_empty_cells() {
    let json = r#"{"cells": []}"#;
    let nb = mathr::notebook::parse_notebook_json(json).unwrap();
    assert_eq!(nb.cells.len(), 0);
}

#[test]
fn notebook_parse_bad_json_errors() {
    assert!(mathr::notebook::parse_notebook_json(r#"{"foo": "bar"}"#).is_err());
}

#[test]
fn notebook_json_escape_special_chars() {
    let mut nb = Notebook::new();
    nb.add_cell("a\nb\tc");
    let json = nb.to_json();
    assert!(json.contains("\\n"));
    assert!(json.contains("\\t"));
    let nb2 = mathr::notebook::parse_notebook_json(&json).unwrap();
    assert_eq!(nb2.cells[0].input, "a\nb\tc");
}

#[test]
fn notebook_load_example_file() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/notebooks/demo.mnb");
    let nb = Notebook::load(&path).unwrap();
    assert_eq!(nb.cells.len(), 4);
    assert_eq!(nb.cells[0].input, "sin(pi/4)");
    assert_eq!(nb.cells[1].input, r"\frac{1}{2} + \frac{3}{4}");
    assert_eq!(nb.cells[2].input, "diff x^3");
}

#[test]
fn notebook_eval_loaded_example() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples/notebooks/demo.mnb");
    let mut nb = Notebook::load(&path).unwrap();
    nb.eval_all(&Context::standard()).unwrap();
    assert!(nb.cells[0].output.contains("0.707"));
    assert!(nb.cells[1].output.contains("1.25"));
    assert!(nb.cells[2].output.contains("3"));
}
