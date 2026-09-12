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
- [x] **Condition number / nullspace** — `Matrix::condition_number` (2-norm, σ_max/σ_min from SVD, inf for singular) and `Matrix::nullspace` (orthonormal basis of {x : Ax = 0} from the eigendecomposition of AᵀA, works for wide matrices). `cond <rows>` and `null <rows>` REPL commands.
- [x] **More special functions** — `digamma`, `trigamma`, `polygamma(m, x)` (via Hurwitz zeta), `harmonic(n)`, `zeta(s)` (Euler–Maclaurin, exact trivial zeros/pole), `hurwitz(s, a)`, elliptic integrals `elliptic_k(k)` (Carlson RF), `elliptic_e(k)` (AGM identity), `elliptic_f(φ, k)` (RF), `elliptic_e_inc(φ, k)` (Simpson quadrature). All registered as expression functions; validated against reference digits (Apéry, π²/6, Legendre K/E values).

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
- [x] **868 inline unit tests** + 226 integration tests — all passing
- [x] AGENTS.md, README.md, ARCHITECTURE.md, SPEC.md

### Fast math
- [x] **Chebyshev-based fast math library** — `ChebyshevApprox` struct, `fast_sin`/`fast_cos`/`fast_tan`/`fast_exp`/`fast_log`/`fast_sqrt`/`fast_pow` with argument reduction, `fast` REPL command

### Big integers
- [x] **Big integer support** — `bigint` module with arbitrary-precision primality (Miller–Rabin), factorization (trial division + Pollard's rho), GCD, LCM, factorial, Fibonacci (fast doubling), binomial, modular exponentiation, totient. REPL commands `fact`/`fib`/`binom` auto-upgrade to BigInt on u64 overflow (SymPy-style). `big` command for explicit big-integer ops on inputs > u64::MAX.

### Automatic differentiation
- [x] **Dual numbers** — `autodiff` module with `Dual` type (value + derivative), full arithmetic operator overloads, elementary functions (sin, cos, tan, exp, ln, log, sqrt, powf, pow_dual, asin, acos, atan, sinh, cosh, tanh, abs), `eval` for `Expr` AST, `derivative`, `gradient`, `jacobian`, `ad` REPL command. Naming follows AD literature conventions.

### Serialization
- [x] **Expression serialization** — `serialize` module with three interchangeable formats: S-expressions (`to_sexpr`/`from_sexpr`), JSON (`to_json`/`from_json`), and RPN (`to_rpn`/`from_rpn`). All round-trip via `Expr::equals`. Hand-rolled JSON parser (no serde dep). `serialize <fmt> <expr>` and `serialize <fmt> import <text>` REPL commands.

### Interval arithmetic
- [x] **Interval arithmetic** — `interval` module with `Interval` type (closed `[lo, hi]` bounds), full arithmetic (+, -, *, /, neg), `sqr`/`powi` with zero-crossing awareness, elementary functions (sin, cos, tan, exp, ln, log10, log2, sqrt, abs, asin, acos, atan, sinh, cosh, tanh, cbrt) with monotonicity and extrema tracking, set operations (intersect, hull, contains, overlaps), and `eval_interval` for rigorous bounds over the `Expr` AST. `interval <expr> with <var>=[lo,hi],...` REPL command. Pure Rust, no deps.

### Matrix decompositions
- [x] **QR decomposition** — Householder reflections for `A = Q·R` where Q is orthogonal and R is upper-triangular. Handles rectangular m×n matrices (tall and wide). Least-squares solve for overdetermined systems. `qr <rows>` REPL command.

### Arbitrary-precision decimals
- [x] **Arbitrary-precision decimal arithmetic** — `bigdec` module on the `bigdecimal` crate: precision-controlled `add`/`sub`/`mul`/`div` (BigInt long division sized from operand digit counts), `sqrt` (Newton), `exp`/`ln` (argument reduction + Taylor/atanh series), `pi` (Machin), `e`, `log2`/`log10`, `sin`/`cos`/`tan` (π/2 quadrant reduction), `atan`/`asin`/`acos`, `sinh`/`cosh`/`tanh`, `pow`, `eval_decimal`/`eval_decimal_rounded` over the `Expr` AST with variable bindings. All transcendentals compute at `prec + 10` guard digits and correctly round. `dec <expr> [prec <n>] [with <var>=<val>,...]` REPL command. Parser-folded `pi`/`e`/`tau` f64 constants are restored to full precision.

### Symbolic limits
- [x] **Symbolic limits** — `limit` module: `lim x→a f(x)` for finite points and ±∞. Three-stage strategy: direct substitution on the simplified expression, L'Hôpital's rule for `0/0`/`∞/∞` quotients (recursive, max 6 applications), numeric probing with geometric steps from both sides. Classifies finite limits (integer/12-sig snapping), poles (`±∞`), opposite-side divergence, and oscillation (`does not exist`). `limit <expr> [<var>] <point>` REPL command with variable inference; step-by-step output via `dispatch_steps` (notebook UI).

### Polynomial algebra
- [x] **Polynomial expansion/collect** — `poly` module: `expand` distributes products and non-negative integer powers into a collected sum of monomials (multivariate, `f64` coefficients, like-term collection, descending-degree + lex ordering, `Sub` joins for negative terms). Non-polynomial parts (functions, symbolic powers, variable denominators) stay intact with expanded children; safety caps on exponent (64) and term count (20k) prevent blowup. `to_poly`/`poly_to_expr` exported as the foundation for partial fractions. `expand <expr>` REPL command + step-by-step.
- [x] **Partial fraction decomposition** — `apart` module: `apart <expr> [var]` decomposes univariate rational functions via polynomial long division + numeric factorization of the denominator (Durand–Kerner complex roots, conjugate pairs → irreducible quadratics, root clustering for multiplicities) + square linear system for the coefficients (solved by Gaussian elimination). Sign-aware output joins (`-1/(x + 1) + 1/x`), integer/12-sig coefficient snapping, step-by-step in the notebook.

### Probability & statistics
- [x] **Probability distributions & hypothesis tests** — `dists` module: Student-t, chi-squared, F, binomial, Poisson, uniform (pdf/cdf/pmf), `normal_ppf` (Acklam + Halley refinement), regularized incomplete beta `beta_inc` (continued fraction). Hypothesis tests returning statistic/df/p-value: one-sample, Welch two-sample, and paired t-tests, chi-squared goodness-of-fit (+uniform), one-way ANOVA. `dist` REPL command extended (t/chi2/f/binom/poisson/uniform); new `ttest`/`chitest`/`anova` commands with `|` group separator.
- [x] **Quantile (inverse CDF) functions** — `student_t_ppf`, `chi2_ppf`, `f_ppf`, `uniform_ppf` via generic bracketed bisection to machine precision (`invert_cdf`); `qtile <dist> <p> <params...>` REPL command for critical values (`normal` accepts optional `[mu sigma]`). Validated against closed forms (Cauchy tan, chi²(2) = −2ln(1−p), F(1,ν) = t²) and standard t/chi² tables.
- [x] **Nonlinear least squares (Levenberg–Marquardt)** — `curvefit` module: `curve_fit` with numeric central-difference Jacobian, Marquardt diagonal scaling, damping sub-loop (λ adaptive), convergence by SSE-drop/step size, parameter standard errors from σ²(JᵀJ)⁻¹. `fit <expr> [in <var>] [with p=v,...] x1 y1 ...` REPL command (params inferred or guessed); `curvefit_demo` example. Validated against closed-form `linear_regression` + Michaelis–Menten/exponential recovery; sine-frequency local-minimum behavior documented.
- [x] **Nonparametric tests & bootstrap CIs** — Mann–Whitney U (tie-corrected variance, continuity correction), Wilcoxon signed-rank (zeros dropped, tie-corrected), Kruskal–Wallis H (chi-squared, tie-corrected), Spearman rank correlation (Pearson on average ranks), percentile bootstrap confidence intervals (seeded, reproducible). REPL commands `mwu`/`wilcoxon`/`kw`/`boot`. Hand-computed rank/tie test values validated; `TestResult.df` now `Option<f64>`.
- [x] **Logistic regression (IRLS)** — `logit` module: `logistic_regression` (Newton/IRLS on the log-likelihood, working response z = Xβ + (y−p)/w, Marquardt-free), Wald standard errors from (XᵀWX)⁻¹, overflow-safe sigmoid and log-likelihood, `predict_proba`. 1e-10 ridge keeps separation well-defined (divergence reported honestly). `logit <y...> with <x1...> [| <x2...>]` REPL command; `logit_demo` example. Validated against closed-form empirical-logit solutions (balanced two-point: b1 = 2ln2, b0 = −3ln2, SE = √3) and zero-score equations at the MLE.
- [x] **Stiff ODE solvers** — `stiff` module: backward Euler (L-stable) + implicit trapezoidal/Crank–Nicolson (A-stable, 2nd order) for systems, per-step Newton with central-difference numeric Jacobians and `Matrix::solve`. `*_system` (final state) and `*_trajectory` variants. Validated by stiff-decay stability contrast (explicit Euler amplifies 9×/step to 1e95; implicit methods exact to 1e-10), observed convergence orders (2.01 / 4.00), harmonic-oscillator amplitude conservation, and Robertson-problem mass conservation to 1e-10. `stiff_demo` example. BDF2 added 2026-09-12 (see Low Priority Done).
- [x] **Monotone interpolation (PCHIP)** — `pchip` module: Fritsch–Carlson piecewise cubic Hermite (`Pchip::new`/`eval`/`slopes`) with non-uniform-grid weighted-harmonic-mean slopes, zero slope at local extrema, endpoint clipping, clamped extrapolation, n=2 degenerate-to-line case. Never overshoots the local knot range (THE property, tested densely), knot-exact, linear data exact, sin accuracy ~O(h³) verified by refinement ratio. `pchip x1 y1 ... x_at` REPL command mirroring `spline`.
- [x] **Optimization (derivative-free minimization)** — `optim` module: `golden_section` (1-D unimodal bracket, φ-shrink, FnMut for Expr closures) and `nelder_mead` (N-D simplex with reflection/expansion/contraction/shrink, `OptResult { x, fx, iterations, converged }`). `minimize <expr> [var] a b` REPL command over expression objectives. Validated on closed-form parabola/sine minima, quadratic bowls up to 3-D, and Rosenbrock's banana valley (converged to 1e-4 from the classic start).
- [x] **B-spline / Hermite interpolation** — `bspline` module: Cox–de Boor basis functions, general `BSpline` (de Boor evaluation, affine parameter map, clamped endpoints), and `cubic_bspline_interp` (de Boor knot averaging → nonsingular collocation, automatic degree reduction for n < 4). Plus `CubicHermite` in `interpolate.rs` (piecewise cubic Hermite with user-supplied slopes). `bspline x1 y1 ... x_at` and `hermite x1 y1 d1 ... x_at` REPL commands; `bspline_demo` example. Validated by exact cubic/quadratic/linear reproduction, sin interpolation (8.7e-4 max error at 8 points), knot/slope exactness.
- [x] **Sparse matrices (CSR/CSC)** — `sparse` module: `Csr`/`Csc` compressed storage from coordinate triplets (duplicates summed, canonical order) or dense matrices, `matvec`, `transpose` (scatter layout conversion), row-wise sparse×sparse `multiply` (SpGEMM), and `conjugate_gradient` (structural symmetry check, non-positive-curvature breakdown detection, relative-residual stopping). `spcg <rows> | <b...>` REPL command; `sparse_demo` example; validated against dense `Matrix::solve` and the 1-D Poisson system.
- [x] **Preconditioned CG + BiCGStab** — `conjugate_gradient_jacobi` (diagonal preconditioning; provably iteration-invariant under diagonal scaling) and `bicgstab` (general nonsymmetric systems, Jacobi-preconditioned, four explicit breakdown guards). REPL: `spcg jacobi ...` flag + `spbicg <rows> | <b...>`; validated against dense solves, convection-diffusion stencils, and singular-inconsistent breakdown. ILU(0) added 2026-09-12 (see Low Priority Done).

## Brainstorming

### High Priority
(none currently — optimization completed this cycle)

### Medium Priority
(none currently)

### Low Priority
- [x] **ILU(0)-preconditioned Krylov solvers** — `Ilu0::factorize` (zero-fill LU over A's sparsity pattern, row-wise IKJ) + triangular solves; `bicgstab_ilu` preconditioned BiCGStab via a generic preconditioner-closure core. Exact (1 iteration) on no-fill matrices like tridiagonals; `spbicg ilu` REPL flag.
- [x] **Stiff ODE solvers (implicit BDF)** — `bdf2_system`/`bdf2_trajectory` in `stiff.rs`: fixed-step BDF2 (A-stable, 2nd order) with one trapezoidal startup step, reusing the shared Newton core; validated by order-2 convergence ratio, stiff-decay stability, and Robertson mass conservation
- [x] **B-spline / Hermite interpolation** — `bspline` module (basis functions, de Boor eval, cubic interpolation via knot averaging) + `CubicHermite`; `bspline`/`hermite` REPL commands
- [x] **Sparse matrices (CSR/CSC)** — compressed storage for large systems: `sparse` module with `Csr`/`Csc`, matvec, transpose, SpGEMM, `conjugate_gradient`, `spcg` REPL command
- [x] **Preconditioned conjugate gradient** — Jacobi (diagonal) preconditioning: `conjugate_gradient_jacobi`, `spcg jacobi` REPL flag
- [x] **BiCGStab for general sparse systems** — nonsymmetric solve without symmetry requirement: `bicgstab`, `spbicg` REPL command
- [x] **Arbitrary-precision arithmetic (BigDecimal)** — `bigdec` module: `pi`/`e`/`sqrt`/`exp`/`ln`/trig to N significant digits, `eval_decimal` over `Expr` AST, `dec` REPL command
- [x] **Automatic differentiation (dual numbers)** — `autodiff` module with `Dual` type, `derivative`, `gradient`, `jacobian`, `ad` REPL command
- [x] **MathML support** — W3C Presentation MathML export (`to_mathml`, `to_mathml_doc`) and import (`from_mathml`). `mathml` REPL command for both directions. Supports `<mn>`, `<mi>`, `<mo>`, `<mrow>`, `<mfrac>`, `<msup>`, `<msub>`, `<msqrt>`, `<mroot>`, `<mfenced>`, `<mstyle>`, `<mtext>`.
- [x] **Expression serialization** — S-expressions, JSON, RPN (`serialize` module, `serialize` REPL command)
- [x] **Interval arithmetic** — rigorous bounds via `Interval` type + `eval_interval` over Expr AST (`interval` module, `interval` REPL command)
- [ ] GPU-accelerated FFT (via `wgpu`)
- [ ] 3D plotting (surface plots, contour plots)
- [ ] Animated plot output (GIF/WebM)