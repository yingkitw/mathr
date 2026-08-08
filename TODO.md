# TODO

## Done

### Core
- [x] Expression AST, parser, evaluator with Context
- [x] Expression equality checking (canonical form with commutative sorting, constant folding)

### Symbolic algebra
- [x] Symbolic differentiation (product, quotient, chain rules)
- [x] **Multi-variable symbolic differentiation** (partial derivatives and gradients)
- [x] Algebraic simplification (constant folding, identities)
- [x] **Symbolic integration** for polynomial, exponential, trigonometric, and inverse-trigonometric primitives

### Numerical calculus
- [x] Numerical derivatives (high-order finite difference) and gradients
- [x] Trapezoidal, Simpson's, and adaptive quadrature
- [x] **Romberg with Richardson extrapolation**
- [x] **Fourier series** (numerical coefficient computation via Simpson's rule, evaluation)
- [x] **Monte Carlo integration** (1-D and N-D, reproducible LCG, standard error)
- [x] **Stochastic primitives** (Rng, uniform/normal/exponential sampling, normal/exp PDF/CDF)
- [x] **Moment / cumulant helpers** (skewness, excess kurtosis, cumulants up to order 4)
- [x] **Hilbert-matrix-aware solvers** (Tikhonov regularisation, Hilbert matrix construction, least-squares for rectangular systems)

### Equation solving
- [x] Bisection, Newton–Raphson, secant
- [x] Durand–Kerner polynomial root finding
- [x] **Polynomial root isolation** (VAS method with i128 exact arithmetic)
- [x] **Newton's method for nonlinear systems** (central-difference Jacobian)

### FFT
- [x] Cooley–Tukey radix-2 from scratch: forward, inverse, 2D, real-input
- [x] Magnitude / power spectra
- [x] Convolution and cross-correlation
- [x] Window functions (Hann, Hamming, Blackman, Rectangular)

### Matrix
- [x] Arithmetic (add, sub, mul, scalar), transpose, trace
- [x] Gaussian elimination: determinant, inverse, linear solve
- [x] **Rank** estimation
- [x] **LU decomposition** with partial pivoting
- [x] **Cholesky decomposition** `A = L · Lᵀ` for symmetric positive-definite matrices
- [x] **SVD** `A = U · Σ · Vᵀ` via one-sided Jacobi rotations
- [x] **Power iteration** for the dominant eigenvalue/eigenvector
- [x] **QR algorithm** for full symmetric eigenvalue decomposition (Householder tridiagonalisation + Wilkinson-shift QR iteration)
- [x] **Hessenberg decomposition** `A = Q·H·Qᵀ` via Householder reflections
- [x] **Real Schur decomposition** `A = Q·T·Qᵀ` via shifted QR on Hessenberg form
- [x] **`det` REPL command** for matrix determinant

### Statistics
- [x] Mean, median, variance, standard deviation
- [x] Quartiles, IQR
- [x] Pearson correlation
- [x] Linear regression

### Number theory
- [x] GCD, LCM, extended GCD, modular inverse, modular exponentiation
- [x] Primality (trial division + Miller–Rabin, deterministic for n < 3.3e24)
- [x] Prime factorization, sieve of Eratosthenes
- [x] Binomial coefficients, factorial, Fibonacci (fast doubling)
- [x] Euler's totient
- [x] **Jacobi symbol**
- [x] **Continued fractions** (rational and real-valued approximants)
- [x] **Linear Diophantine solver**
- [x] **Discrete logarithm** (baby-step giant-step)
- [x] Chinese Remainder Theorem

### ODE
- [x] Euler, RK4, RK4 systems
- [x] Adaptive RKF45

### Interpolation
- [x] Lagrange, Newton divided-difference, linear
- [x] **Natural / clamped cubic spline** (Thomas algorithm for second derivatives)
- [x] **Chebyshev polynomials** `T_n(x)`, Chebyshev nodes, series approximation, Clenshaw evaluation
- [x] **Legendre polynomials** `P_n(x)`, associated `P_n^m(x)`, **Gauss–Legendre quadrature**

### Special functions
- [x] Gamma, log-Gamma (Lanczos)
- [x] Beta function
- [x] erf, erfc (via incomplete gamma)
- [x] sinc, incomplete gamma P
- [x] **Bessel functions** `J_0(x)`, `J_1(x)`, integer-order `J_n(x)` (Maclaurin series + asymptotic + forward recurrence)

### Taylor series
- [x] Symbolic expansion around any point
- [x] **Laurent series** (expansion around poles, negative powers, principal + analytic parts)

### Rational arithmetic
- [x] **Rational number type** (exact arithmetic, GCD reduction, i128 intermediate, parsing, REPL)

### Web notebook
- [x] **Math notebook** (`.mnb` file format, JSON cells with TeX/math input + output)
- [x] **Web notebook server** (minimal HTTP server, Jupyter-like UI with KaTeX rendering, cell eval, save/load)
- [x] **Step-by-step solving** (`dispatch_steps` — shows intermediate steps for diff, solve, taylor, integrate, simplify, rat, laurent)
- [x] **Exact rational evaluation in notebook** (fraction expressions evaluated as `Rational`, returning exact fractions instead of decimals)
- [x] **KaTeX math rendering** (input preview + output rendering, plain-to-LaTeX converter, decimal-to-fraction display)
- [x] **Fallback to simplify on unbound variables** (expressions with variables that can't be evaluated are simplified instead of erroring)
- [x] **`det` command in notebook/REPL** (matrix determinant via `det <rows>`)
- [x] **Shared context across cells** (`let`/`fn` bindings persist across cells, `/api/reset` + `/api/context` endpoints, `dispatch_with_ctx` API)
- [x] **Cell types** (Math cells evaluated with KaTeX, Text cells for documentation — `CellType` enum, JSON `cell_type` field)
- [x] **Cell management** (move up/down, duplicate, toggle type, execution status indicators, context panel, Reset & Run All)
- [x] **Inline plots** (`plot` commands render PNG images directly in the notebook via base64, `plot_function_to_bytes`/`plot_multi_to_bytes`/`plot_scatter_to_bytes`)
- [x] **Markdown text cells** (text cells render Markdown via marked.js — headings, lists, code, blockquotes)
- [x] **Execution counters** (`In [n]:` indicators like Jupyter, Alt+Enter to run + add cell)

### Other
- [x] Complex number type with arithmetic, polar conversion, powers
- [x] PNG plotting (line, multi-series, scatter) via `plotters`
- [x] LaTeX / TeX input (`\frac`, `\sqrt`, `\sin`, `\pi`, `\left(\right)`, `^{...}`, `\Gamma`, `\log_2`, …; `$...$`, `$$...$$`, `\[...\]`, `\(...\)`)
- [x] Interactive REPL (rustyline-powered with history)
- [x] CLI subcommands and REPL dispatch for all features
- [x] **406 inline unit tests** + 156 integration tests — all passing
- [x] AGENTS.md, README.md, ARCHITECTURE.md, SPEC.md

### Fast math
- [x] **Chebyshev-based fast math library** — `ChebyshevApprox` struct, `fast_sin`/`fast_cos`/`fast_tan`/`fast_exp`/`fast_log`/`fast_sqrt`/`fast_pow` with argument reduction, `fast` REPL command

### Big integers
- [x] **Big integer support** — `bigint` module with arbitrary-precision primality (Miller–Rabin), factorization (trial division + Pollard's rho), GCD, LCM, factorial, Fibonacci (fast doubling), binomial, modular exponentiation, totient. REPL commands `fact`/`fib`/`binom` auto-upgrade to BigInt on u64 overflow (SymPy-style). `big` command for explicit big-integer ops on inputs > u64::MAX.

### Automatic differentiation
- [x] **Dual numbers** — `autodiff` module with `Dual` type (value + derivative), full arithmetic operator overloads, elementary functions (sin, cos, tan, exp, ln, log, sqrt, powf, pow_dual, asin, acos, atan, sinh, cosh, tanh, abs), `eval` for `Expr` AST, `derivative`, `gradient`, `jacobian`, `ad` REPL command. Naming follows AD literature conventions.

## Brainstorming

### High Priority
(none currently)

### Medium Priority
(none currently — big integer support completed)

### Low Priority
- [ ] Arbitrary-precision arithmetic (BigDecimal)
- [ ] Interval arithmetic for rigorous bounds
- [x] **Automatic differentiation (dual numbers)** — `autodiff` module with `Dual` type, `derivative`, `gradient`, `jacobian`, `ad` REPL command
- [ ] Expression serialization (S-expressions, JSON, RPN)
- [ ] GPU-accelerated FFT (via `wgpu`)
- [ ] 3D plotting (surface plots, contour plots)
- [ ] Animated plot output (GIF/WebM)