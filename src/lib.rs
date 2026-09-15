//! # mathr
//!
//! A pure-Rust mathematics library and CLI that bundles symbolic and
//! numerical computation:
//!
//! - **Expression AST, parser and evaluator** ([`expr`], [`parser`], [`eval`])
//!   with a standard library of elementary functions (`sin`, `cos`, `exp`,
//!   `log`, ...) and constants (`pi`, `e`, ...).
//! - **FFT from scratch** ([`fft`]): Cooley–Tukey radix-2 (forward, inverse,
//!   2D, real-input) plus magnitude / power spectra, convolution,
//!   cross-correlation, and window functions (Hann, Hamming, Blackman).
//! - **Numerical calculus** ([`calculus`]): high-order finite-difference
//!   derivatives, trapezoidal and Simpson's integrators, adaptive quadrature,
//!   partial derivatives and gradients.
//! - **Symbolic algebra** ([`symbolic`], [`simplify`]): differentiation by
//!   the usual calculus rules, plus a simplifier that flattens identities and
//!   folds constants.  Expression equality via canonical form ([`expr`]).
//! - **Equation solving** ([`solver`]): bisection, Newton–Raphson (with
//!   numeric-derivative fallback), secant, and Durand–Kerner for polynomial
//!   roots.
//! - **Plotting** ([`plot`]): PNG output via [`plotters`], with single-,
//!   multi- and scatter-plot variants.
//! - **Matrix operations** ([`matrix`]): arithmetic, determinant, inverse,
//!   linear system solving, trace, LU decomposition with partial pivoting,
//!   rank estimation, Cholesky decomposition, SVD, power iteration, and
//!   symmetric eigenvalue decomposition via the QR algorithm, Hessenberg
//!   decomposition, and real Schur decomposition.
//! - **Sparse matrices** ([`sparse`]): CSR/CSC compressed storage with
//!   matrix–vector products, transposition, sparse×sparse multiplication,
//!   Jacobi-preconditioned conjugate gradient for symmetric positive-definite
//!   systems, and BiCGStab for general nonsymmetric systems.
//! - **Statistics** ([`stats`]): mean, median, variance, stddev, quartiles,
//!   correlation, linear regression.
//! - **Curve fitting** ([`curvefit`]): Levenberg–Marquardt nonlinear least
//!   squares with numeric Jacobians, parameter standard errors.
//! - **Logistic regression** ([`logit`]): binary classification via IRLS with
//!   Wald standard errors and stable log-likelihood.
//! - **Distributions & tests** ([`dists`]): Student-t, chi-squared, F,
//!   binomial, Poisson, uniform pdf/cdf, normal quantile, regularized
//!   incomplete beta, and classical hypothesis tests (t-tests,
//!   chi-squared goodness-of-fit, one-way ANOVA).
//! - **Number theory** ([`numtheory`]): GCD, LCM, primality, factorization,
//!   binomial coefficients, Fibonacci, sieve, Euler's totient, Miller–Rabin,
//!   Chinese Remainder Theorem, modular exponentiation, Jacobi symbol,
//!   continued fractions, linear Diophantine solver.
//! - **ODE solvers** ([`ode`]): Euler, RK4, RK4 systems, adaptive RKF45.
//! - **Stiff ODE solvers** ([`stiff`]): backward Euler and implicit
//!   trapezoidal (Crank–Nicolson) for stiff systems — A-stable, per-step
//!   Newton with numeric Jacobians.
//! - **Taylor series** ([`taylor`]): symbolic Taylor expansion around a point.
//! - **Laurent series** ([`laurent`]): expansion around poles with negative powers.
//! - **Rational arithmetic** ([`rational`]): exact `Rational` type with GCD reduction, arithmetic, parsing.
//! - **Notebook** ([`notebook`]): `.mnb` file format with cells of TeX/math expressions and results.
//! - **Web notebook server** ([`server`]): minimal HTTP server serving a Jupyter-like web UI.
//! - **Interpolation** ([`interpolate`]): Lagrange, Newton divided-difference,
//!   linear interpolation, and natural / clamped cubic splines.
//! - **Monotone interpolation** ([`pchip`]): piecewise cubic Hermite
//!   (Fritsch–Carlson) that never overshoots the local data range.
//! - **B-spline interpolation** ([`bspline`]): Cox–de Boor basis functions,
//!   de Boor curve evaluation, and clamped interpolating splines via knot
//!   averaging (with automatic degree reduction for few points).
//! - **Optimization** ([`optim`]): derivative-free minimization — golden-
//!   section search in 1-D and Nelder–Mead simplex in N dimensions.
//! - **Special functions** ([`special`]): Gamma, log-Gamma, Beta, erf, erfc,
//!   sinc, incomplete gamma P, Bessel functions `J_0`, `J_1`, `J_n`.
//!
//! The CLI (`mathr "<expr or command>"`) is a thin wrapper around the same
//! library functions. It accepts plain math expressions, LaTeX/TeX input,
//! and command keywords like `diff`, `solve`, `int`, `gcd`, etc.

pub mod calculus;
pub mod ceval;
pub mod complex;
pub mod curvefit;
pub mod error;
pub mod eval;
pub mod expr;
pub mod fft;
pub mod interpolate;
pub mod logit;
pub mod pchip;
pub mod stiff;
pub mod matrix;
pub mod numtheory;
pub mod optim;
pub mod bigint;
pub mod bspline;
pub mod autodiff;
pub mod ode;
pub mod parser;
pub mod plot;
pub mod repl;
pub mod simplify;
pub mod solver;
pub mod sparse;
pub mod special;
pub mod dists;
pub mod stats;
pub mod symbolic;
pub mod taylor;
pub mod laurent;
pub mod rational;
pub mod fastmath;
pub mod notebook;
pub mod server;
pub mod mathml;
pub mod serialize;
pub mod interval;
pub mod bigdec;
pub mod limit;
pub mod poly;
pub mod apart;
pub mod qsolve;
pub mod sumprod;

pub use error::{MathError, Result};

/// Re-exports of the most common types for downstream `use mathr::*;`.
pub mod prelude {
    pub use crate::complex::Complex;
    pub use crate::error::{MathError, Result};
    pub use crate::eval::{eval, eval_str, Context, Func};
    pub use crate::expr::Expr;
    pub use crate::parser::Parser;
    pub use crate::simplify::simplify;
    pub use crate::matrix::{Cholesky, EigenPair, Lu, Matrix, PowerIterOptions, Svd};
    pub use crate::interpolate::{
        chebyshev_coefficients, chebyshev_eval, chebyshev_nodes, chebyshev_rescale, chebyshev_t,
        gauss_legendre, lagrange_interp, legendre_associated, legendre_p, lerp, newton_interp,
        CubicSpline, NewtonInterpolator,
    };
    pub use crate::numtheory::{
        binomial, chinese_remainder, continued_fraction, continued_fraction_value, diophantine,
        discrete_log, euler_totient, extended_gcd, factorial, fibonacci, gcd, is_prime,
        is_prime_miller_rabin, jacobi_symbol, lcm, mod_inverse, mod_pow, prime_factors,
        sieve_primes,
    };
    pub use crate::autodiff::{
        derivative as ad_derivative, eval as dual_eval, gradient as ad_gradient,
        jacobian as ad_jacobian, Dual,
    };
    pub use crate::mathml::{from_mathml, to_mathml, to_mathml_doc};
    pub use crate::ceval::{eval_complex, eval_complex_str};
    pub use crate::qsolve::solve_symbolic;
    pub use crate::sumprod::{product, summation};
    pub use crate::serialize::{from_json, from_rpn, from_sexpr, to_json, to_rpn, to_sexpr};
    pub use crate::interval::{Interval, eval_interval};
    pub use crate::bigdec::{eval_decimal, BigDecimal};
    pub use crate::limit::{limit, LimitValue};
    pub use crate::poly::expand;
    pub use crate::apart::apart;
    pub use crate::curvefit::{curve_fit, LmFit, LmOptions};
    pub use crate::logit::{logistic_regression, predict_proba, LogitFit, LogitOptions as LogitOpts};
    pub use crate::stiff::{bdf2_system, bdf2_trajectory, implicit_euler_system, implicit_euler_trajectory, trapezoid_system, trapezoid_trajectory};
    pub use crate::pchip::Pchip;
    pub use crate::bspline::{basis_function, cubic_bspline_interp, BSpline};
    pub use crate::interpolate::CubicHermite;
    pub use crate::optim::{golden_section, nelder_mead, OptOptions, OptResult};
    pub use crate::ode::{euler, rk4, rk4_system, rkf45};
    pub use crate::solver::{
        bisect, isolate_real_roots, newton_central, newton_system, polynomial_roots, secant, SolveOptions,
    };
    pub use crate::special::{bessel_j0, bessel_j1, bessel_jn, beta, erfc, erf, gamma, log_gamma, sinc};
    pub use crate::calculus::{
        derivative, fourier_eval, fourier_series, integrate_adaptive, integrate_romberg,
        integrate_simpson, integrate_trap, monte_carlo_integrate_1d, monte_carlo_integrate_nd,
        partial, second_derivative, FourierSeries,
    };
    pub use crate::stats::{
        correlation, cumulants, exp_cdf, exp_pdf, linear_regression, mean, median, moments,
        normal_cdf, normal_pdf, Rng, stddev, variance, Summary,
    };
    pub use crate::dists::{
        anova_oneway, beta_inc, binomial_cdf, binomial_pmf, bootstrap_ci, chi2_cdf, chi2_pdf,
        chi2_ppf, chi_square_gof, chi_square_uniform, f_cdf, f_pdf, f_ppf, kruskal_wallis,
        mann_whitney_u, normal_ppf, poisson_cdf, poisson_pmf, spearman_corr, student_t_cdf,
        student_t_pdf, student_t_ppf, t_test_one, t_test_paired, t_test_two, uniform_cdf,
        uniform_pdf, uniform_ppf, wilcoxon_signed_rank, TestResult,
    };
    pub use crate::symbolic::{differentiate, gradient, integrate};
    pub use crate::taylor::taylor_series;
    pub use crate::laurent::{laurent_series, laurent_series_str, LaurentSeries};
    pub use crate::rational::{parse_rational, Rational};
    pub use crate::notebook::{CellType, Notebook, NotebookCell};
    pub use crate::server::NotebookServer;
    pub use crate::fastmath::{fast_cos, fast_exp, fast_log, fast_pow, fast_sin, fast_sqrt, fast_tan, ChebyshevApprox};
    pub use crate::sparse::{bicgstab, bicgstab_ilu, conjugate_gradient, conjugate_gradient_jacobi, CgResult, Csc, Csr, Ilu0};
}