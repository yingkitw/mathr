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
//!   linear system solving, trace.
//! - **Statistics** ([`stats`]): mean, median, variance, stddev, quartiles,
//!   correlation, linear regression.
//! - **Number theory** ([`numtheory`]): GCD, LCM, primality, factorization,
//!   binomial coefficients, Fibonacci, sieve, Euler's totient, Miller–Rabin,
//!   Chinese Remainder Theorem, modular exponentiation.
//! - **ODE solvers** ([`ode`]): Euler, RK4, RK4 systems, adaptive RKF45.
//! - **Taylor series** ([`taylor`]): symbolic Taylor expansion around a point.
//! - **Interpolation** ([`interpolate`]): Lagrange, Newton divided-difference,
//!   linear interpolation.
//! - **Special functions** ([`special`]): Gamma, log-Gamma, Beta, erf, erfc,
//!   sinc, incomplete gamma P.
//!
//! The CLI (`mathr "<expr or command>"`) is a thin wrapper around the same
//! library functions. It accepts plain math expressions, LaTeX/TeX input,
//! and command keywords like `diff`, `solve`, `int`, `gcd`, etc.

pub mod calculus;
pub mod complex;
pub mod error;
pub mod eval;
pub mod expr;
pub mod fft;
pub mod interpolate;
pub mod matrix;
pub mod numtheory;
pub mod ode;
pub mod parser;
pub mod plot;
pub mod repl;
pub mod simplify;
pub mod solver;
pub mod special;
pub mod stats;
pub mod symbolic;
pub mod taylor;

pub use error::{MathError, Result};

/// Re-exports of the most common types for downstream `use mathr::*;`.
pub mod prelude {
    pub use crate::complex::Complex;
    pub use crate::error::{MathError, Result};
    pub use crate::eval::{eval, eval_str, Context, Func};
    pub use crate::expr::Expr;
    pub use crate::parser::Parser;
    pub use crate::simplify::simplify;
    pub use crate::matrix::Matrix;
    pub use crate::interpolate::{lagrange_interp, newton_interp, NewtonInterpolator};
    pub use crate::numtheory::{binomial, chinese_remainder, factorial, fibonacci, gcd, is_prime, is_prime_miller_rabin, lcm, mod_pow, prime_factors, sieve_primes};
    pub use crate::ode::{euler, rk4, rk4_system, rkf45};
    pub use crate::solver::{bisect, newton_central, polynomial_roots, secant, SolveOptions};
    pub use crate::special::{beta, erfc, erf, gamma, log_gamma, sinc};
    pub use crate::stats::{correlation, linear_regression, mean, median, stddev, variance, Summary};
    pub use crate::symbolic::differentiate;
    pub use crate::taylor::taylor_series;
}