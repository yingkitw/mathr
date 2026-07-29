# mathr

A pure-Rust mathematics library and CLI for symbolic and numerical computation — FFT from scratch, calculus, equation solving, plotting, matrix operations, statistics, number theory, ODE solvers, Taylor series, interpolation, and special functions.

## Quick Start

```bash
cargo build --release
```

### CLI Usage

Just give it a string — `mathr` figures out what to do:

```bash
mathr "sin(pi/4) + 2^3"          # evaluate an expression
mathr "gamma(0.5)"               # special functions as expressions
mathr "diff x^3 + 2*x^2"         # symbolic derivative
mathr "simplify 2*x + 3*x + 0"   # simplify
mathr "solve x^2 - 4"            # find roots
mathr "int sin(x) 0 pi"          # numerical integral
mathr "taylor exp(x) 0 5"        # Taylor series
mathr "poly-roots 1 -5 6"        # polynomial roots
mathr "fft 1 0 -1 0 1 0 -1 0"    # FFT magnitude spectrum
mathr "conv 1 2 3 x 1 1"         # convolution (x separates signals)
mathr "stats 1 2 3 4 5 6 7 8"    # descriptive statistics
mathr "gcd 48 36"                # number theory
mathr "is-prime 97"
mathr "factor 360"
mathr "fib 50"
mathr "mr-prime 2305843009213693951"
mathr "plot sin(x) -6.28 6.28 wave.png"
echo "sin(pi/2)" | mathr         # read from stdin
mathr                            # interactive REPL
```

### TeX / LaTeX Input

`mathr` accepts LaTeX math formulas — with or without Markdown delimiters:

```bash
mathr "\frac{1}{2} + \frac{3}{4}"           # → 1.25
mathr "$\sin(\pi / 4)$"                     # → 0.7071...  (inline math)
mathr "$$\sqrt{16} + \cos(\pi)$$"           # → 3          (display math)
mathr "\[\frac{x^2 - 4}{1}\]"               # → solve with \[...\]
mathr "\(\log_2{8}\)"                       # → 3          (\(...\) inline)
mathr "2 \cdot 3 + 4"                       # → 10
mathr "\left( 1 + 2 \right) \cdot 3"        # → 9
mathr "\Gamma{0.5}"                         # → 1.7724...
mathr "diff \sin(x^2)"                      # → cos(x^2)*2*x
mathr "solve $\frac{x^2 - 4}{1}$"           # → root ≈ 2
```

Supported delimiters: `$...$`, `$$...$$`, `\[...\]`, `\(...\)`, or raw TeX
with no delimiters at all.

Supported TeX commands: `\frac`, `\sqrt`, `\pi`, `\tau`, `\infty`, `\cdot`,
`\times`, `\left(`, `\right)`, `^{...}`, `\sin`, `\cos`, `\tan`, `\arcsin`,
`\arccos`, `\arctan`, `\sinh`, `\cosh`, `\tanh`, `\exp`, `\ln`, `\log`,
`\log_2`, `\log_{10}`, `\Gamma`, `\operatorname{...}`, `\text{...}`.

### REPL

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
mathr> taylor exp(x) 0 5
1 + x + 0.5*x^2 + 0.1666666667*x^3 + 0.0416666667*x^4
mathr> quit
```

## Features

- **Expression parsing & evaluation**: Recursive-descent parser, implicit multiplication, scientific notation, constants (`pi`, `e`, `tau`, `inf`), 27+ built-in functions
- **Symbolic differentiation**: Product, quotient, chain rules with simplification
- **Algebraic simplification**: Constant folding, identity flattening
- **Numerical calculus**: 5-point stencil derivatives, trapezoidal/Simpson's/adaptive quadrature, partial derivatives, gradients
- **Equation solving**: Bisection, Newton–Raphson, secant, Durand–Kerner polynomial roots
- **FFT from scratch**: Cooley–Tukey radix-2 (forward, inverse, 2D, real-input), magnitude/power spectra, convolution, cross-correlation, window functions (Hann, Hamming, Blackman)
- **Complex numbers**: Generic `Complex<T>` with arithmetic, `abs`, `arg`, `exp`, `from_polar`
- **Plotting**: PNG output via `plotters` (line, multi-series, scatter)
- **Matrix operations**: Arithmetic, determinant, inverse, linear system solve, trace, transpose
- **Statistics**: Mean, median, variance, stddev, quartiles, IQR, correlation, linear regression
- **Number theory**: GCD, LCM, extended GCD, modular inverse, primality, factorization, sieve, binomial, factorial, Fibonacci, Euler's totient, Miller–Rabin, Chinese Remainder Theorem, modular exponentiation
- **ODE solvers**: Euler, RK4, RK4 systems, adaptive RKF45
- **Taylor series**: Symbolic Taylor expansion around any point
- **Interpolation**: Lagrange, Newton divided-difference, linear
- **Special functions**: Gamma, log-Gamma, Beta, erf, erfc, sinc, incomplete gamma P
- **Expression equality**: Canonical form comparison with commutative sorting and constant folding
- **Interactive REPL**: rustyline-powered with history, bracket matching, variable/function bindings

## Crate API

```rust
use mathr::prelude::*;

// Evaluate
let val = eval_str("sin(pi/4) + 2^3", &[])?;

// Symbolic differentiation
let expr = Parser::parse("x^3 + 2*x")?;
let deriv = differentiate(&expr, "x")?;

// FFT
let samples = vec![1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
let mags = mathr::fft::magnitude_spectrum(&samples)?;

// Matrix
let m = Matrix::from_rows(&[vec![1.0, 2.0], vec![3.0, 4.0]])?;
let det = m.determinant()?;
let inv = m.inverse()?;

// Statistics
let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
let s = mathr::stats::summary(&data)?;

// Number theory
let primes = mathr::numtheory::sieve_primes(100);
let fib50 = mathr::numtheory::fibonacci(50);

// ODE
let y = mathr::ode::rk4(|_t, y| y, 0.0, 1.0, 1.0, 100)?;

// Taylor series
let series = mathr::taylor::taylor_series_str("exp(x)", "x", 0.0, 5)?;

// Interpolation
let pts = vec![(0.0, 1.0), (1.0, 2.0), (2.0, 5.0)];
let y = lagrange_interp(&pts, 0.5)?;

// Special functions
let g = mathr::special::gamma(0.5); // √π
let e = mathr::special::erf(1.0);

// Miller–Rabin primality
let is_prime = mathr::numtheory::is_prime_miller_rabin(2305843009213693951, 20);
```

## Dependencies

| Crate | Purpose |
|-------|---------|
| `clap` | CLI argument parsing |
| `anyhow` | Error handling in binary |
| `thiserror` | Error types in library |
| `rustyline` | REPL line editing |
| `plotters` | PNG plotting |
| `num-traits` | Numeric traits |
| `approx` (dev) | Float comparison in tests |

## License

Apache-2.0
