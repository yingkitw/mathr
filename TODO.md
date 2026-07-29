# TODO

## Done

- [x] Expression AST, parser, evaluator with Context
- [x] Symbolic differentiation (product, quotient, chain rules)
- [x] Algebraic simplification (constant folding, identities)
- [x] Numerical calculus (derivatives, trapezoidal/Simpson's/adaptive quadrature, gradients)
- [x] Equation solvers (bisection, Newton–Raphson, secant, Durand–Kerner)
- [x] FFT from scratch (Cooley–Tukey radix-2, inverse, 2D, real-input, spectra)
- [x] Complex number type from scratch
- [x] PNG plotting (single, multi, scatter)
- [x] Interactive REPL (rustyline)
- [x] Matrix operations (arithmetic, determinant, inverse, solve, trace, transpose)
- [x] Descriptive statistics (mean, median, variance, stddev, quartiles, correlation, regression)
- [x] Number theory (GCD, LCM, extended GCD, mod inverse, primality, factorization, sieve, binomial, factorial, Fibonacci, Euler's totient)
- [x] ODE solvers (Euler, RK4, RK4 systems, adaptive RKF45)
- [x] Taylor series expansion (symbolic, around any point)
- [x] Polynomial interpolation (Lagrange, Newton divided differences, linear)
- [x] FFT convolution and cross-correlation
- [x] Window functions (Hann, Hamming, Blackman, Rectangular)
- [x] Miller–Rabin primality test (deterministic for n < 3.3e24)
- [x] Modular exponentiation, Chinese Remainder Theorem
- [x] Special functions (Gamma, log-Gamma, Beta, erf, erfc, sinc, incomplete gamma P)
- [x] Expression equality checking (canonical form with commutative sorting, constant folding)
- [x] CLI subcommands for all features
- [x] REPL `taylor` command
- [x] 129 inline unit tests — all passing
- [x] AGENTS.md, README.md, ARCHITECTURE.md, SPEC.md

## Brainstorming

### High Priority
- [ ] Eigenvalue/eigenvector computation (power iteration, QR algorithm)
- [ ] LU decomposition for matrix operations
- [ ] SVD (singular value decomposition)
- [ ] Cubic spline interpolation
- [ ] Chebyshev polynomial approximation
- [ ] Symbolic integration for simple cases (polynomial, exponential)

### Medium Priority
- [ ] Multi-variable symbolic differentiation (partial derivatives)
- [ ] Series expansion (Laurent series, Fourier series)
- [ ] Continued fractions
- [ ] Rational number type (exact arithmetic)
- [ ] Big integer support for number theory (arbitrary precision)
- [ ] Discrete logarithm
- [ ] Bessel functions, Legendre polynomials
- [ ] Newton's method for systems of nonlinear equations
- [ ] Romberg integration
- [ ] Monte Carlo integration

### Low Priority
- [ ] Arbitrary-precision arithmetic (BigDecimal)
- [ ] Interval arithmetic for rigorous bounds
- [ ] Automatic differentiation (dual numbers)
- [ ] Expression serialization (S-expressions, JSON)
- [ ] LaTeX output for expressions
- [ ] WebAssembly build for browser use
- [ ] Python bindings (PyO3)
- [ ] GPU-accelerated FFT (via wgpu)
- [ ] 3D plotting (surface plots, contour plots)
- [ ] Animated plot output (GIF/WebM)
