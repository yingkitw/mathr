//! Integration tests that exercise multiple modules together end-to-end.

use mathr::complex::Complex;
use mathr::eval::{eval, Context};
use mathr::expr::Expr;
use mathr::fft;
use mathr::interpolate;
use mathr::matrix::Matrix;
use mathr::numtheory;
use mathr::parser::Parser;
use mathr::ode;
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
    ctx.funcs.insert("square".into(), mathr::eval::Func::Builtin(square));
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
