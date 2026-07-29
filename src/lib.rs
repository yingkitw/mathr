//! # maths
//!
//! A pure-Rust mathematics library and CLI that bundles symbolic and
//! numerical computation:
//!
//! - **Expression AST, parser and evaluator** ([`expr`], [`parser`], [`eval`])
//!   with a standard library of elementary functions (`sin`, `cos`, `exp`,
//!   `log`, ...) and constants (`pi`, `e`, ...).
//! - **FFT from scratch** ([`fft`]): Cooley–Tukey radix-2 (forward, inverse,
//!   2D, real-input) plus magnitude / power spectra.
//! - **Numerical calculus** ([`calculus`]): high-order finite-difference
//!   derivatives, trapezoidal and Simpson's integrators, adaptive quadrature,
//!   partial derivatives and gradients.
//! - **Symbolic algebra** ([`symbolic`], [`simplify`]): differentiation by
//!   the usual calculus rules, plus a simplifier that flattens identities and
//!   folds constants.
//! - **Equation solving** ([`solver`]): bisection, Newton–Raphson (with
//!   numeric-derivative fallback), secant, and Durand–Kerner for polynomial
//!   roots.
//! - **Plotting** ([`plot`]): PNG output via [`plotters`], with single-,
//!   multi- and scatter-plot variants.
//!
//! The CLI (`maths eval|diff|integrate|solve|plot|fft|repl`) is a thin
//! wrapper around the same library functions.

pub mod calculus;
pub mod complex;
pub mod error;
pub mod eval;
pub mod expr;
pub mod fft;
pub mod parser;
pub mod plot;
pub mod repl;
pub mod simplify;
pub mod solver;
pub mod symbolic;

pub use error::{MathError, Result};

/// Re-exports of the most common types for downstream `use maths::*;`.
pub mod prelude {
    pub use crate::complex::Complex;
    pub use crate::error::{MathError, Result};
    pub use crate::eval::{eval, eval_str, Context, Func};
    pub use crate::expr::Expr;
    pub use crate::parser::Parser;
    pub use crate::simplify::simplify;
    pub use crate::solver::{bisect, newton_central, polynomial_roots, secant, SolveOptions};
    pub use crate::symbolic::differentiate;
}