# mathr — Rust Math Library, CLI Calculator & Web Notebook

[![crates.io](https://img.shields.io/crates/v/mathr.svg)](https://crates.io/crates/mathr)
[![docs.rs](https://docs.rs/mathr/badge.svg)](https://docs.rs/mathr)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)

**mathr** is a pure-Rust mathematics library and command-line calculator for symbolic and numerical computation — built from scratch with zero external math dependencies. It includes a Jupyter-like web notebook with KaTeX math rendering, step-by-step solving, and exact fraction arithmetic.

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [CLI Usage](#cli-usage)
- [LaTeX / TeX Input](#latex--tex-input)
- [REPL](#repl)
- [Web Notebook](#web-notebook)
- [Crate API](#crate-api)
- [Modules](#modules)
- [Dependencies](#dependencies)
- [License](#license)

## Features

### Symbolic Computation
- **Symbolic differentiation** — product, quotient, chain rules; partial derivatives and gradients
- **Symbolic integration** — polynomial, exponential, trigonometric, and inverse-trigonometric primitives
- **Algebraic simplification** — constant folding and algebraic identity simplification
- **Taylor series** — symbolic expansion around any point
- **Laurent series** — expansion around poles with negative powers
- **Exact rational arithmetic** — `Rational` type with GCD reduction; expressions with integer fractions evaluate exactly (e.g. `\frac{1}{2} + \frac{3}{4}` → `5/4`)

### Numerical Computation
- **Numerical calculus** — high-order finite-difference derivatives, trapezoidal / Simpson's / adaptive quadrature, **Romberg with Richardson extrapolation**
- **FFT** — Cooley–Tukey radix-2: forward, inverse, 2D, real-input, magnitude / power spectra, convolution, cross-correlation, window functions (Hann, Hamming, Blackman, Rectangular)
- **Equation solving** — bisection, Newton–Raphson, secant, Durand–Kerner polynomial roots, **Newton's method for nonlinear systems**, VAS root isolation
- **ODE solvers** — Euler, RK4, RK4 systems, adaptive RKF45
- **Monte Carlo integration** — 1-D and N-D with reproducible LCG and standard error
- **Fourier series** — numerical coefficient computation via Simpson's rule

### Linear Algebra
- **Matrix operations** — arithmetic, determinant, inverse, linear-system solver, trace, transpose, rank
- **LU decomposition** with partial pivoting
- **Cholesky decomposition** `A = L·Lᵀ` for SPD matrices
- **SVD** `A = U · Σ · Vᵀ` via one-sided Jacobi rotations
- **Eigenvalue solvers** — power iteration, symmetric QR algorithm (Householder + Wilkinson shift)
- **Hessenberg** and **real Schur decomposition**
- **Tikhonov regularisation** for ill-conditioned and rectangular systems

### Number Theory
- GCD, LCM, primality (trial + **Miller–Rabin**), factorization, sieve of Eratosthenes
- Binomial coefficients, factorial, Fibonacci (fast doubling), Euler's totient
- **Chinese Remainder Theorem**, modular exponentiation, modular inverse
- **Jacobi symbol**, **continued fractions**, **linear Diophantine solver**, **discrete logarithm** (baby-step giant-step)
- **Big integers** — arbitrary-precision primality, factorization (Pollard's rho), factorial, Fibonacci, binomial, modular exponentiation, totient via `num-bigint`. REPL commands `fact`, `fib`, `binom` auto-upgrade to BigInt on overflow (no separate "big" command needed for these)
- **Automatic differentiation** — dual numbers for exact forward-mode AD; derivatives, gradients, and Jacobians of arbitrary compositions

### Interpolation & Special Functions
- **Interpolation** — Lagrange, Newton, linear, **cubic spline**, **Chebyshev** polynomials and series, **Legendre** polynomials, **Gauss–Legendre quadrature**
- **Special functions** — Gamma, log-Gamma, Beta, erf, erfc, sinc, incomplete gamma P, **Bessel functions** `J_0`, `J_1`, `J_n`
- **Fast math** — Chebyshev-based approximations of `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `pow` with argument reduction (~1e-12 accuracy)

### Input & Output
- **LaTeX / TeX input** — parse `\frac`, `\sqrt`, `\sin`, `\pi`, `\left(\right)`, `^{...}`, `\Gamma`, `\log_2`, and more; supports `$...$`, `$$...$$`, `\[...\]`, `\(...\)` delimiters
- **Interactive REPL** — rustyline-powered with history, variable/function bindings
- **Expression parser** — recursive-descent, implicit multiplication, scientific notation, 30+ built-in functions
- **PNG plotting** — line, multi-series, scatter via `plotters`
- **Web notebook** — Jupyter-like UI with KaTeX math rendering, step-by-step solving, exact fraction arithmetic, `.mnb` file format, cell-based evaluation

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
mathr "det 1 2 | 3 4"                                # matrix determinant

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

# Web notebook (Jupyter-like UI)
mathr notebook                 # open web UI at http://127.0.0.1:3000
mathr notebook examples/notebooks/demo.mnb   # load a notebook file
mathr notebook examples/notebooks/demo.mnb 8080  # custom port
```

## LaTeX / TeX Input

`mathr` accepts LaTeX math formulas — with or without Markdown delimiters:

```bash
mathr "\frac{1}{2} + \frac{3}{4}"           # → 5/4 (exact fraction)
mathr "$\sin(\pi / 4)$"                     # → 0.7071...  (inline $...$)
mathr "$$\sqrt{16} + \cos(\pi)$$"           # → 3          (display $$...$$)
mathr "\[\frac{x^2 - 4}{1}\]"               # → evaluate with \[...\]
mathr "\(\log_2{8}\)"                       # → 3          (\(...\) inline)
mathr "2 \cdot 3 + 4"                       # → 10
mathr "\left( 1 + 2 \right) \cdot 3"        # → 9
mathr "\Gamma{0.5}"                         # → 1.7724...  (√π)
mathr "diff \sin(x^2)"                      # → cos(x^2)*2*x
mathr "solve $\frac{x^2 - 4}{1}$"           # → root ≈ 2
mathr "\frac{x^2 - 1}{x - 1}"              # → (x^2 - 1)/(x - 1) (simplify on unbound vars)
```

**Supported delimiters**: `$...$`, `$$...$$`, `\[...\]`, `\(...\)`, or raw TeX with no delimiters.

**Supported TeX commands**: `\frac`, `\sqrt`, `\pi`, `\tau`, `\infty`, `\cdot`, `\times`, `\left(`, `\right)`, `^{...}`, `\sin`, `\cos`, `\tan`, `\arcsin`, `\arccos`, `\arctan`, `\sinh`, `\cosh`, `\tanh`, `\exp`, `\ln`, `\log`, `\log_2`, `\log_{10}`, `\Gamma`, `\operatorname{...}`, `\text{...}`.

## REPL

The interactive REPL supports variable bindings, function definitions, and all commands:

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
mathr> det 1 2 | 3 4
det = -2
mathr> fast sin 1.5
fast sin(1.5) = 0.997495  (exact: 0.997495, err: 0e0)
mathr> big prime 1000000007
1000000007 is prime
mathr> big fact 25
25! = 15511210043330985984000000
mathr> big fib 100
F_100 = 354224848179261915075
mathr> fact 25
15511210043330985984000000
mathr> fib 100
354224848179261915075
mathr> binom 100 50
100891344545564193334812497256
mathr> factor 360
2^3 · 3^2 · 5
mathr> 5!
120
mathr> 20!
2432902008176640000
mathr> |-5|
5
mathr> 7 mod 3
1
mathr> gcd(12, 8)
4
mathr> C(5, 2)
10
mathr> \binom{5}{2}
10
mathr> ad sin(x^2) at x=1.5
f(x) = 0.997495,  f'(x) = -0.313312
mathr> ad grad x^2 + y^3 with x=2,y=3
∇f = [∂f/∂x = 4, ∂f/∂y = 27]
mathr> big modpow 2 100 1000000007
2^100 ≡ 97637128 (mod 1000000007)
mathr> quit
```

## Web Notebook

Launch a Jupyter-like web notebook with KaTeX math rendering:

```bash
mathr notebook                        # open web UI at http://127.0.0.1:3000
mathr notebook examples/notebooks/demo.mnb      # load a notebook file
mathr notebook examples/notebooks/demo.mnb 8080 # custom port
```

### Example Notebooks

| File | Topics |
|------|--------|
| `examples/notebooks/demo.mnb` | General overview |
| `examples/notebooks/calculus.mnb` | Differentiation, integration, Taylor series, gradients |
| `examples/notebooks/fractions.mnb` | Exact fraction arithmetic with `\frac` and `rat` |
| `examples/notebooks/solving.mnb` | Root finding, polynomial roots, simplification |
| `examples/notebooks/linear_algebra.mnb` | LU, Cholesky, SVD, eigenvalues, FFT, stats |
| `examples/notebooks/number_theory.mnb` | GCD, primality, factorization, Diophantine, discrete log |
| `examples/notebooks/special_functions.mnb` | Gamma, erf, Bessel, sinc |
| `examples/notebooks/latex_demo.mnb` | LaTeX/TeX input with `\sin`, `\frac`, `\Gamma` |
| `examples/notebooks/series_interp.mnb` | Taylor/Laurent series, splines, Chebyshev, Legendre |
| `examples/notebooks/notebook_features.mnb` | Shared context, text cells, inline plots, Markdown, execution counters |

### Notebook Features

- **Shared context across cells** — variables and functions defined in one cell (via `let`/`fn`) are available in subsequent cells, just like Jupyter
- **Cell types** — Math cells (evaluated with KaTeX rendering) and Text cells (Markdown-rendered documentation)
- **Inline plots** — `plot` commands render PNG images directly in the notebook (no file management)
- **Cell management** — add, delete, duplicate, move up/down, and toggle cell type
- **Execution status** — each cell shows running/done/error status with `In [n]:` execution counters
- **Context panel** — collapsible panel showing all bound variables and user functions
- **Reset & Run All** — resets the shared context and re-evaluates all cells in order
- **Markdown rendering** — text cells render Markdown (headings, lists, code, blockquotes) via marked.js
- **KaTeX math rendering** — input expressions and output results rendered as math notation
- **Step-by-step solving** — shows intermediate steps for `diff`, `solve`, `taylor`, `integrate`, `simplify`, `rat`, `laurent`
- **Exact fraction arithmetic** — `\frac{1}{2} + \frac{3}{4}` evaluates to `5/4`, not `1.25`
- **Live input preview** — each cell shows a rendered math/Markdown preview as you type
- **Save / load** — `.mnb` JSON file format with cells of TeX/math input and evaluated output
- **Keyboard shortcuts** — Shift/Cmd/Ctrl+Enter to run a cell, Alt+Enter to run and add a new cell

### `.mnb` File Format

```json
{
  "cells": [
    { "id": 0, "input": "let x = 5", "output": "x = 5", "cell_type": "math" },
    { "id": 1, "input": "x * 3", "output": "= 15", "cell_type": "math" },
    { "id": 2, "input": "# My Notes", "output": "# My Notes", "cell_type": "text" }
  ]
}
```

### REST API

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/` | Serve web UI HTML |
| `POST` | `/api/eval` | Evaluate expression (updates shared context); returns `{input, output, steps}` |
| `GET` | `/api/notebook` | Get current notebook as JSON |
| `POST` | `/api/notebook` | Replace notebook state (auto-saves to file) |
| `POST` | `/api/save` | Save notebook to file |
| `POST` | `/api/reset` | Reset the shared evaluation context |
| `GET` | `/api/context` | Get current variables and user functions |

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
| `calculus` | Numerical derivatives, quadrature, gradients, Romberg, Monte Carlo, Fourier series |
| `solver` | Bisection, Newton, secant, polynomial roots, **Newton for systems**, VAS root isolation |
| `fft` | Cooley–Tukey FFT, convolution, cross-correlation, windows |
| `complex` | Generic complex number type |
| `matrix` | Matrix arithmetic, determinant, inverse, solve, **LU**, **Cholesky**, **SVD**, **eigenvalues**, rank |
| `stats` | Descriptive statistics, correlation, regression, stochastic primitives |
| `numtheory` | GCD, LCM, primality, factorization, sieve, CRT, totient, **Jacobi**, **Diophantine**, **continued fractions**, **discrete log** |
| `bigint` | Arbitrary-precision integers: primality (Miller–Rabin), factorization (Pollard's rho), GCD, LCM, factorial, Fibonacci, binomial, mod_pow, totient |
| `autodiff` | Automatic differentiation via dual numbers: `derivative`, `gradient`, `jacobian`, `Dual` type with full arithmetic |
| `ode` | Euler, RK4, RK4 systems, adaptive RKF45 |
| `taylor` | Symbolic Taylor series expansion |
| `laurent` | Laurent series expansion around poles |
| `interpolate` | Lagrange, Newton, linear, **cubic spline**, **Chebyshev**, **Legendre**, **Gauss–Legendre** |
| `special` | Gamma, Beta, erf, erfc, sinc, incomplete gamma, **Bessel J_0/J_1/J_n** |
| `fastmath` | Chebyshev-based fast approximations of `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `pow` |
| `rational` | Exact rational arithmetic (`Rational` type), `eval_rational` for exact AST evaluation |
| `notebook` | `.mnb` notebook format, JSON cells with TeX/math input + output, **cell types** (math/text), **cell reordering**, **shared context** |
| `server` | Minimal HTTP server for web notebook UI (KaTeX rendering, step-by-step solving) |
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
| `num-bigint` | Arbitrary-precision integers |
| `num-integer` | Integer trait methods (mod_floor, is_even) |
| `approx` (dev) | Float comparison in tests |

## License

Apache-2.0

## Links

- [crates.io](https://crates.io/crates/mathr)
- [docs.rs](https://docs.rs/mathr)
- [GitHub](https://github.com/yingkitw/mathr)
- [Report an issue](https://github.com/yingkitw/mathr/issues)