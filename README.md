# mathr — Pure-Rust Math Library & CLI

[![crates.io](https://img.shields.io/crates/v/mathr.svg)](https://crates.io/crates/mathr)
[![docs.rs](https://docs.rs/mathr/badge.svg)](https://docs.rs/mathr)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

A pure-Rust mathematics library and command-line tool for **symbolic and numerical computation** — built from scratch with zero external math dependencies.

## Features

- **Symbolic algebra** — differentiation (product, quotient, chain rules) and **symbolic integration** (polynomial, exponential, trigonometric, and inverse-trigonometric primitives) with automatic simplification
- **Numerical calculus** — high-order finite-difference derivatives, trapezoidal / Simpson's / adaptive quadrature, **Romberg with Richardson extrapolation**, partial derivatives and gradients
- **FFT** — Cooley–Tukey radix-2 from scratch: forward, inverse, 2D, real-input, magnitude / power spectra, convolution, cross-correlation, window functions (Hann, Hamming, Blackman, Rectangular)
- **Equation solving** — bisection, Newton–Raphson, secant, Durand–Kerner polynomial root finding, and **Newton's method for nonlinear systems**
- **Matrix operations** — arithmetic, determinant, inverse, linear-system solver, trace, transpose, **rank** estimation, **LU decomposition** with partial pivoting, **Cholesky decomposition** `A = L·Lᵀ` for SPD matrices, **SVD** `A = U · Σ · Vᵀ`, **dominant eigenvalue/eigenvector** via power iteration
- **Statistics** — mean, median, variance, standard deviation, quartiles, IQR, correlation, linear regression
- **Number theory** — GCD, LCM, primality (trial + Miller–Rabin), factorization, sieve of Eratosthenes, binomial coefficients, factorial, Fibonacci, Euler's totient, Chinese Remainder Theorem, modular exponentiation, **Jacobi symbol**, **continued fractions**, **linear Diophantine solver**, **discrete logarithm** (baby-step giant-step)
- **ODE solvers** — Euler, RK4, RK4 systems, adaptive RKF45
- **Taylor series** — symbolic expansion around any point
- **Interpolation** — Lagrange, Newton divided-difference, linear, **natural / clamped cubic spline**, **Chebyshev polynomials** and series approximation, **Legendre polynomials** and associated, **Gauss–Legendre quadrature**
- **Special functions** — Gamma, log-Gamma, Beta, erf, erfc, sinc, incomplete gamma P, **Bessel functions** `J_0`, `J_1`, `J_n`
- **Complex numbers** — generic `Complex<T>` with arithmetic, polar conversion, powers
- **Plotting** — PNG output via `plotters` (line, multi-series, scatter)
- **LaTeX / TeX input** — parse `\frac`, `\sqrt`, `\sin`, `\pi`, `\left(\right)`, `^{...}`, `\Gamma`, `\log_2`, and more; supports `$...$`, `$$...$$`, `\[...\]`, `\(...\)` delimiters
- **Interactive REPL** — rustyline-powered with history, bracket matching, variable/function bindings
- **Expression parser** — recursive-descent, implicit multiplication, scientific notation, 30+ built-in functions

## Quick Start

```bash
cargo install mathr
```

Or build from source:

```bash
git clone https://github.com/yingkitw/mathr.git
cd mathr
cargo build --release
```

## CLI Usage

Just pass a string — `mathr` figures out what to do:

```bash
mathr "sin(pi/4) + 2^3"               # evaluate an expression
mathr "gamma(0.5)"                    # special functions
mathr "diff x^3 + 2*x^2"              # symbolic derivative
mathr "integrate sin(x)"              # symbolic integration → -cos(x)
mathr "simplify 2*x + 3*x + 0"        # algebraic simplification
mathr "solve x^2 - 4"                 # find roots (Newton–Raphson)
mathr "int sin(x) 0 pi"               # numerical integral over [0, π]
mathr "romberg sin(x) 0 3.14159"      # Romberg-integral of sin on [0, π]
mathr "taylor exp(x) 0 5"             # Taylor series (5 terms around 0)
mathr "poly-roots 1 -5 6"             # polynomial roots (Durand–Kerner)
mathr "fft 1 0 -1 0 1 0 -1 0"         # FFT magnitude spectrum
mathr "conv 1 2 3 x 1 1"              # convolution (x separates signals)
mathr "stats 1 2 3 4 5 6 7 8"         # descriptive statistics

# Matrix operations (rows separated by `|`)
mathr "lu 1 2 3 | 4 -6 0 | -2 7 2"                    # LU decomposition
mathr "cholesky 4 12 -16 | 12 37 -43 | -16 -43 98"    # Cholesky decomposition
mathr "eig 2 1 | 1 2"                                 # dominant eigenpair
mathr "svd 1 2 | 3 4 | 5 6"                           # SVD (rectangular OK)

# Interpolation
mathr "spline 0 0 1 1 2 4 3 9 1.5"    # cubic spline at x=1.5
mathr "chebyshev 5 0.3"                # Chebyshev T_5(0.3)
mathr "legendre 5 0.3"                 # Legendre P_5(0.3)

# Number theory
mathr "gcd 48 36"                                         # GCD
mathr "is-prime 97"                                       # primality test
mathr "factor 360"                                        # prime factorization
mathr "fib 50"                                            # Fibonacci
mathr "mr-prime 2305843009213693951"                     # Miller–Rabin
mathr "jacobi 5 7"                                        # Jacobi symbol
mathr "cf 22 7"                                           # continued fraction
mathr "diophantine 3 5 7"                                 # linear Diophantine
mathr "dlog 2 27 101"                                     # discrete log

# Special functions
mathr "bessel_j(2, 5)"                                   # Bessel J_2(5)

# Plot to PNG
mathr "plot sin(x) -6.28 6.28 wave.png"

echo "sin(pi/2)" | mathr       # read from stdin
mathr                         # interactive REPL
```

## LaTeX / TeX Input

`mathr` accepts LaTeX math formulas — with or without Markdown delimiters:

```bash
mathr "\frac{1}{2} + \frac{3}{4}"           # → 1.25
mathr "$\sin(\pi / 4)$"                     # → 0.7071...  (inline $...$)
mathr "$$\sqrt{16} + \cos(\pi)$$"           # → 3          (display $$...$$)
mathr "\[\frac{x^2 - 4}{1}\]"               # → evaluate with \[...\]
mathr "\(\log_2{8}\)"                       # → 3          (\(...\) inline)
mathr "2 \cdot 3 + 4"                       # → 10
mathr "\left( 1 + 2 \right) \cdot 3"        # → 9
mathr "\Gamma{0.5}"                         # → 1.7724...  (√π)
mathr "diff \sin(x^2)"                      # → cos(x^2)*2*x
mathr "solve $\frac{x^2 - 4}{1}$"           # → root ≈ 2
```

**Supported delimiters**: `$...$`, `$$...$$`, `\[...\]`, `\(...\)`, or raw TeX with no delimiters.

**Supported TeX commands**: `\frac`, `\sqrt`, `\pi`, `\tau`, `\infty`, `\cdot`, `\times`, `\left(`, `\right)`, `^{...}`, `\sin`, `\cos`, `\tan`, `\arcsin`, `\arccos`, `\arctan`, `\sinh`, `\cosh`, `\tanh`, `\exp`, `\ln`, `\log`, `\log_2`, `\log_{10}`, `\Gamma`, `\operatorname{...}`, `\text{...}`.

## REPL

```
mathr> sin(pi/4)
0.7071067812
mathr> let x = 3
mathr> x^2 + 1
10
mathr> fn f(x) = x^2 + 2*x + 1
mathr> f(5)
36
mathr> diff sin(x^2)
2*x*cos(x^2)
mathr> integrate sin(x)
-cos(x)
mathr> taylor exp(x) 0 5
1 + x + 0.5*x^2 + 0.1666666667*x^3 + 0.0416666667*x^4
mathr> romberg sin(x) 0 3.14159
2
mathr> svd 1 2 | 3 4 | 5 6
σ = [9.525518, 0.514301]
mathr> quit
```

## Crate API

```rust
use mathr::prelude::*;

// Evaluate an expression
let val = eval_str("sin(pi/4) + 2^3", &[])?;

// Symbolic differentiation and integration
let expr = Parser::parse("x^3 + 2*x")?;
let deriv = differentiate(&expr, "x")?;
let integr = integrate(&Parser::parse("x^2")?, "x")?;   // x³/3

// FFT magnitude spectrum
let samples = vec![1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
let mags = mathr::fft::magnitude_spectrum(&samples)?;

// Matrix operations
let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]])?;
let det = m.determinant()?;
let inv = m.inverse()?;
let lu = m.lu()?;
let chol = m.cholesky()?;
let svd = m.svd()?;                              // A = U Σ Vᵀ
let eigen = m.power_iteration(PowerIterOptions::default())?;

// Descriptive statistics
let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
let s = mathr::stats::summary(&data)?;

// Number theory
let primes = mathr::numtheory::sieve_primes(100);
let fib50 = mathr::numtheory::fibonacci(50);
let is_prime = mathr::numtheory::is_prime_miller_rabin(2305843009213693951, 20);
let cf = mathr::numtheory::continued_fraction(22, 7)?;       // [3; 7]
let j  = mathr::numtheory::jacobi_symbol(5, 7);                // -1
let (x, y) = mathr::numtheory::diophantine(3, 5, 7)?;          // (14, -7)
let dl = mathr::numtheory::discrete_log(2, 27, 101);            // Some(7)

// ODE: solve y' = y, y(0) = 1, on [0, 1]
let y = mathr::ode::rk4(|_t, y| y, 0.0, 1.0, 1.0, 100)?;

// Taylor series expansion
let series = mathr::taylor::taylor_series_str("exp(x)", "x", 0.0, 5)?;

// Interpolation
let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0)];
let y = lagrange_interp(&pts, 0.5)?;
let sp = CubicSpline::new(&[(0.0, 0.0), (1.0, 1.0), (2.0, 4.0), (3.0, 9.0)])?;
let v = sp.eval(1.5);
// Chebyshev series
let coeffs = chebyshev_coefficients(|x| x.sin(), 8);
let y_eval = chebyshev_eval(&coeffs, 0.5);
// Gauss–Legendre quadrature
let (nodes, weights) = gauss_legendre(8);

// Special functions
let g = mathr::special::gamma(0.5);
let j0 = mathr::special::bessel_j0(5.0);
let j5 = mathr::special::bessel_jn(5, 2.0);

// Numerical integration
let integral = mathr::calculus::integrate_romberg(|x| x.exp(), 0.0, 1.0, 10)?;

// Newton's method for nonlinear systems
let system = |x: &[f64]| vec![
    2.0 * x[0] + x[1] - 5.0,
    x[0] + 3.0 * x[1] - 7.0,
];
let sol = mathr::solver::newton_system(system, &[0.0, 0.0], SolveOptions::default())?;
```

## Modules

| Module | Description |
|--------|-------------|
| `expr` | Expression AST (`Expr`) with canonicalization and equality |
| `parser` | Recursive-descent parser with LaTeX/TeX support |
| `eval` | Tree-walking evaluator with `Context` (variables, functions) |
| `simplify` | Constant folding and algebraic identity simplification |
| `symbolic` | Symbolic differentiation and integration |
| `calculus` | Numerical derivatives, quadrature, gradients, Romberg |
| `solver` | Bisection, Newton, secant, polynomial roots, **Newton for systems** |
| `fft` | Cooley–Tukey FFT, convolution, cross-correlation, windows |
| `complex` | Generic complex number type |
| `matrix` | Matrix arithmetic, determinant, inverse, solve, **LU**, **Cholesky**, **SVD**, **power iteration**, rank |
| `stats` | Descriptive statistics, correlation, regression |
| `numtheory` | GCD, LCM, primality, factorization, sieve, CRT, totient, **Jacobi**, **Diophantine**, **continued fractions**, **discrete log** |
| `ode` | Euler, RK4, RK4 systems, adaptive RKF45 |
| `taylor` | Symbolic Taylor series expansion |
| `interpolate` | Lagrange, Newton, linear, **cubic spline**, **Chebyshev**, **Legendre**, **Gauss–Legendre** |
| `special` | Gamma, Beta, erf, erfc, sinc, incomplete gamma, **Bessel J_0/J_1/J_n** |
| `plot` | PNG plotting via `plotters` |

## Dependencies

| Crate | Purpose |
|-------|---------|
| `clap` | CLI argument parsing |
| `anyhow` | Error handling in binary |
| `thiserror` | Error types in library |
| `rustyline` | REPL line editing with history |
| `plotters` | PNG plot rendering |
| `num-traits` | Numeric trait bounds |
| `approx` (dev) | Float comparison in tests |

## License

Apache-2.0