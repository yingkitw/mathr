# TODO

## Done

### Core
- [x] Expression AST, parser, evaluator with Context
- [x] Expression equality checking (canonical form with commutative sorting, constant folding)

### Symbolic algebra
- [x] Symbolic differentiation (product, quotient, chain rules)
- [x] **Multi-variable symbolic differentiation** (partial derivatives and gradients)
- [x] **Implicit differentiation** (`idiff` / `implicit_diff`: `dy/dx = −F_x/F_y` for `F(x,y)=0`)
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
- [x] **Mobile/tablet adaptive web UI** (CSS breakpoints at 900px/600px: wrapping toolbar, full-width cell actions, 16px inputs to avoid iOS zoom, safe-area insets)

### Other
- [x] Complex number type with arithmetic, polar conversion, powers
- [x] PNG plotting (line, multi-series, scatter) via `plotters`
- [x] LaTeX / TeX input (`\frac`, `\sqrt`, `\sin`, `\pi`, `\left(\right)`, `^{...}`, `\Gamma`, `\log_2`, …; `$...$`, `$$...$$`, `\[...\]`, `\(...\)`)
- [x] Interactive REPL (rustyline-powered with history)
- [x] CLI subcommands and REPL dispatch for all features
- [x] **Modularize repl.rs** — the 4,962-line monolith is now a `repl/` directory: `mod.rs` (~2,190 lines: the REPL loop, `dispatch_str`/`dispatch_steps`/`dispatch_with_ctx`, `dispatch_inner` routing, shared parsing helpers, `format_value`, `HELP`, and all 77 dispatch tests) plus five domain modules holding the 78 `do_*` handlers verbatim: `calc.rs` (13 calculus handlers), `solve.rs` (39 solving/numerics/linear-algebra/signal handlers, incl. `parse_matrix`), `stats.rs` (13 statistics handlers + `format_test_result`), `numtheory.rs` (5 handlers), `io.rs` (6 export/import handlers). Pure code motion — the public API (`dispatch_str`, `dispatch_steps`, `dispatch_with_ctx`, `run`, `format_value`, `guess_var`) is unchanged, so all call sites and tests needed no edits; child modules see the parent's private helpers through normal descendant visibility. Caught during the move: `split_assignments` initially landed in no file, and `do_sumprod`'s doc comment landed in io.rs — both fixed. (2026-10-31)
- [x] **Modularize the remaining 2,000+ line files (`matrix.rs`, `mcp.rs`)** — `matrix.rs` (2,747 lines, one `impl Matrix` block) is now a `matrix/` directory: `mod.rs` (~1,530: the struct, constructors, core ops, det/solve/inverse/rank, `Display`, operator impls, `SwapWithin`, and all ~250 tests) plus `decomp.rs` (LU/Cholesky/SVD/QR/Tikhonov/Hilbert as additional inherent `impl Matrix` blocks + the `Lu`/`Cholesky`/`Svd`/`Qr` result types) and `eigen.rs` (Hessenberg/Schur/symmetric-eig/power-iteration/condition/nullspace + `givens` + `EigenPair`/`PowerIterOptions`). `mcp.rs` (2,194 lines) is now `mcp/`: `mod.rs` (JSON-RPC protocol, stdio loop, tests), `json.rs` (the hand-rolled JSON layer), `tools.rs` (the 16 tool schemas + execution). Rust allows inherent `impl` blocks in any crate module, so both public APIs — `crate::matrix::*` method calls and `crate::mcp::*` items — are unchanged with **zero call-site edits**. The parent-private-field rule (descendants see ancestor privates) keeps `Matrix.data` accessible from the children. (2026-10-31)
- [x] **1302 inline unit tests** + 360 integration tests — all passing; `cargo clippy --all-targets` exits 0; doc-sync enforced by `tests/docs_sync.rs` (AGENTS.md tree, ARCHITECTURE.md mentions, TODO.md counts)
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

- [x] **Summation & product (Σ/∏)** — `sumprod` module: `summation`/`product` over inclusive integer bounds in three tiers — O(1) closed forms (var-independent, arithmetic series, Faulhaber k²/k³ via i128-exact formulas, geometric `c^var`), exact rational term loop (`eval_rational` accumulation cross-checked against the f64 pass because i64-backed Rationals silently overflow on denominator-heavy sums like Σ1/k past k≈43), and f64 fallback. `sum`/`prod <expr> [<var>] <a> <b>` REPL commands (minimize-style trailing bounds); empty ranges = 0 / 1; 1M-term cap. (2026-09-15, from the qalculate competitive pass)
- [x] **Doc-sync enforcement** — `tests/docs_sync.rs` makes doc drift a test failure: AGENTS.md project structure must list every `pub mod` in `lib.rs`, and the TODO.md Done test counts must match actual `#[test]` occurrences in `src/` and `tests/` (counted as attribute-led lines; the claim must match what `cargo test` reports, including the checker itself). Wired into AGENTS.md Steps 2, 6, and the audit scorecard (2026-09-15).

## Brainstorming

### Codebase Scorecard (2026-09-15, /30 per AGENTS.md Step 7)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 920 unit + 236 integration + 13 doc-tests green; examples compile |
| Numerical fidelity | 5 | Validated vs reference values throughout (see MEMORY.md) |
| Maintainability | 2 → 4 | Was: `cargo clippy` hard-failed (8 approx_constant errors once all targets compiled) + ~60 warnings. Fixed: all errors, gate exits 0, mechanical warnings auto-fixed; remaining ~40 stylistic warnings in frozen mature modules documented as accepted debt. Bonus: `identical if blocks` lint exposed a real silent-wrong-answer parsing bug in `romberg` (fixed + 3 regression tests) |
| Docs alignment | 3 | TODO.md test counts were stale (fixed 2026-09-15); AGENTS.md structure missing 7 newer modules (fixed 2026-09-15); SPEC/README otherwise aligned |
| Wiring & ergonomics | 5 | 40 modules declared in lib.rs, prelude + REPL dispatch complete |
| Footprint | 4 | 9 runtime deps, all justified |

Next weakest: Docs alignment → now enforced: `tests/docs_sync.rs` (2026-09-15) makes module-list/test-count drift a `cargo test` failure.

### Codebase Scorecard — second pass (2026-09-15, /30)
| Dimension | Score | Delta vs baseline |
|---|---|---|
| Correctness & coverage | 5 | — (937 unit + 247 integration + 13 doc-tests; +19 tests since baseline) |
| Numerical fidelity | 5 | — (new modules validated via closed forms: Euler identity, asin(2), qsolve back-substitution) |
| Maintainability | 4 | +2 (clippy gate exits 0; 53 stylistic warnings in frozen modules = documented debt per Mature Module Policy) |
| Docs alignment | 5 | +2 (AGENTS tree + ARCHITECTURE mentions + TODO counts all enforced by `docs_sync`; SPEC/README updated per feature) |
| Wiring & ergonomics | 5 | — (44 modules, prelude exports, dispatch + steps complete) |
| Footprint | 4 | — (9 deps unchanged; ceval/qsolve added zero new deps) |

**Total: 28/30** (baseline 19/30). Remaining 2 points: dependency footprint at its practical ceiling, clippy style debt frozen by policy. Next improvement lever: competitive intelligence to seed features users feel (the 3 speculative Low items stay parked).

### Codebase Scorecard — third pass (2026-10-03, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 953 unit + 259 integration + 13 doc-tests green; idiff + Div Display regression covered |
| Numerical fidelity | 5 | Circle `dy/dx = −x/y` eval checks; coefficient cancel preserves value |
| Maintainability | 4 | Clippy gate still exits 0; mature-module Display fix was one-line + regression test |
| Docs alignment | 5 | docs_sync green; SPEC/README/ARCHITECTURE/TODO/MEMORY updated for `idiff` + responsive web UI |
| Wiring & ergonomics | 5 | `idiff` in REPL, `dispatch_steps`, prelude, webui command list, calculus_demo |
| Footprint | 4 | No new deps |

### Codebase Scorecard — fourth pass (2026-10-13, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 1014 unit (+61) + 259 integration + 13 doc-tests + 4 docs_sync green; +21 new like-term/coefficient tests; 1 apart test relaxed for canonicalize display flip (output unchanged mathematically) |
| Numerical fidelity | 5 | Like-term collection is purely structural; no numeric drift |
| Maintainability | 4 | simplify.rs: doc-comment drift fixed; new helper `split_coefficient`/`rebuild_term` keep `collect` short |
| Docs alignment | 5 | MEMORY.md + AGENTS.md unchanged in scope; TODO + count updated; docs_sync all green |
| Wiring & ergonomics | 5 | No new REPL command (pure library improvement to `simplify()`) |
| Footprint | 4 | No new deps |

**Total: 28/30**. Remaining 2 points: clippy style debt in mature modules frozen by policy. Next lever: High-priority continued-fraction rational snapping (single shared `snap_rational` for qsolve/limit/apart) or Medium-priority sign-scan all-roots numeric solve (both thin extensions).

### Codebase Scorecard — fifth pass (2026-10-13, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 1042 unit (+22 lean tests) + 268 tests/ (+4 pipeline REPL tests for lean, +4 matrix REPL tests) — all green; 14 doc-tests; no regressions |
| Numerical fidelity | 5 | Lean exporter f64 → ℝ is exact-rational (sound); matrix REPL commands sit over existing `Matrix` public API (no new numerics) |
| Maintainability | 4 | Matrix REPL parity landed in `repl.rs` without touching mature `matrix.rs` (symmetry check in REPL layer per Mature Module Policy); `repl.rs` gained a symmetry-tolerance convention |
| Docs alignment | 5 | MEMORY.md "Lean 4 Export Patterns" + new "pass #2" cross-port entry; SPEC.md command table +5; README updated; TODO.md done entries + counts; AGENTS.md tree (lean.rs); ARCHITECTURE.md (lean, MCP server gap noted) |
| Wiring & ergonomics | 5 | 5 new matrix REPL commands (`inv`/`transpose`/`trace`/`eigen`/`linsolve`); `lean <expr>` REPL command (28 unit + 5 integration tests); `examples/lean_demo.rs` |
| Footprint | 4 | No new deps; the Lean exporter is pure string formatting |

**Total: 28/30**. Remaining 2 points: clippy style debt in mature modules frozen by policy. Next lever: the **MCP server** gap (surfaced in pass #2; highest-impact remaining external-surface work — mathr's `dispatch_str` is the near-natural back-end for ≈14 tools).

**Total: 28/30**. Next lever: Medium backlog (dot/cross, number-base) — high user-feel, thin surface area.

### Codebase Scorecard — sixth pass (2026-10-30, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 1119 unit (+42 MCP) + 287 tests/ (+11 MCP protocol) + 4 docs_sync + 14 doc-tests — all green; 30 examples compile; no regressions in the 4 behaviours ported from unsolve |
| Numerical fidelity | 5 | MCP tools route through `repl::dispatch_str`, so results are identical to the CLI by construction; exact rationals preserved (`1/2 + 3/4` → `5/4`); no new numerics |
| Maintainability | 4 | `mcp.rs` is a leaf module: 16 tools are thin wrappers over the REPL, so the MCP surface cannot drift from the CLI. One new public symbol (`repl::format_value`, reused rather than a second formatter) |
| Docs alignment | 5 | MEMORY.md gained "MCP Server Patterns"; SPEC.md gained a full MCP section (transport, handshake, tool table, error codes, client config); README gained an MCP section + module row; AGENTS.md tree (`mcp.rs`); ARCHITECTURE.md module list + design decision 18e; TODO.md counts |
| Wiring & ergonomics | 5 | `mathr mcp` CLI subcommand; 16 tools; `examples/mcp_demo.rs`; prelude exports (`call_tool`, `serve_stdio`, `McpJson`) |
| Footprint | 4 | **No new dependencies** — hand-rolled JSON layer inside `mcp.rs`, matching the `serialize.rs` convention; 9 runtime deps unchanged |

**Total: 28/30**. The **MCP server gap** — flagged in passes #2 and #3 as the highest-impact remaining external-surface work — is now **closed**. Remaining 2 points: clippy style debt in mature modules (frozen by the Mature Module Policy) and the dependency footprint at its practical ceiling. Next lever: the remaining Medium backlog (selective exact function folding, sign-scan all-roots solve, Unicode input normalization, missing LaTeX forms) — all seeded from the unsolve parity passes.

### Codebase Scorecard — seventh pass (2026-11-01, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 1300 unit + 355 integration + 4 docs_sync tests green; 14 doc-tests pass (1 ignored); all 30 examples compile |
| Numerical fidelity | 5 | Durand–Kerner cleanup validated on exact factorizations; leading trimming is scale-relative; downstream `apart` multiplicity suite remains green |
| Maintainability | 4 | Three small solver helpers, no public API breakage; `cargo clippy --all-targets` exits 0 with the established frozen-module style debt |
| Docs alignment | 5 | README/SPEC/ARCHITECTURE/TODO/MEMORY updated; all four docs_sync checks pass |
| Wiring & ergonomics | 5 | Existing `poly-roots`, qsolve, and apart paths share the polished implementation; no dispatch changes needed |
| Footprint | 4 | No new dependencies or modules |

**Total: 28/30**. The last general root-polish item is closed. The remaining parser-precedence change is blocked on explicit approval under the Mature Module Policy; GPU/3D/animated plotting remain speculative rather than evidence-based improvements.

### Codebase Scorecard — eighth pass (2026-11-01, /30)
| Dimension | Score | Evidence |
|---|---|---|
| Correctness & coverage | 5 | 1302 unit + 356 integration + 4 docs_sync tests green; 14 doc-tests pass (1 ignored); all 30 examples compile |
| Numerical fidelity | 5 | Parser now follows conventional sign/power semantics; signed exponents and explicit negative bases are regression-tested |
| Maintainability | 4 | One localized `parse_factor` change; stale assumptions removed from poly/autodiff; mature-module policy approval recorded |
| Docs alignment | 5 | Grammar, examples, architecture note, TODO item, and MEMORY entries updated for the compatibility break |
| Wiring & ergonomics | 5 | All REPL/MCP/notebook expression paths inherit the corrected parser behavior without dispatch changes |
| Footprint | 4 | No dependencies or module growth |

**Total: 28/30**. All pre-existing actionable backlog items are closed. The 2026-11-01 competitive-intelligence pass replaced speculative GPU/3D work with four evidence-backed candidates: arbitrary-size FFT, BFGS, Brent root finding, and hyperdual second derivatives.

### High Priority
(seeded from the unsolve cross-port pass, 2026-10-03 — see MEMORY.md "Cross-Port Learnings from unsolve")
- [x] **Linear-substitution integration** — `∫f(a·v+b) dx = F(a·v+b)/a` for `exp`, `sin`/`cos`, `sinh`/`cosh`, `tan` (Func arm via `linear_in`/`shift_var`); `(a·v+b)^n` and `(a·v+b)^(-1)` (Pow arm); `c/(a·v+b) → (c/a)·ln|a·v+b|` (constant-over-denominator arm); `sqrt(a·v+b)` and `cbrt(a·v+b)` (Func arm via power-of-(ax+b) form). All logarithmic antiderivatives wrapped in `abs` so `ln|x|` (not `ln(x)`) — correct on both branches. 12 regression tests in `symbolic.rs::tests::integrate_linear_subst_*` covering sin/cos/exp/sinh/cosh/power/reciprocal/power-negative/reciprocal-coeff/reciprocal-neg-var/ln-abs/quarter/sqrt. (2026-10-03)
- [x] **Continued-fraction rational snapping** — `Rational::snap_rational` / `snap_coeff_expr` / `SNAP_MAX_DEN` (+ compact `10_000` for qsolve); `limit`, `qsolve`, `apart` share CF snap + 12-sig fallback. (2026-10-19; unsolve `limit.ts:40-66`)
- [x] **Expression-level cubic/quartic solving** — `qsolve` degree ≤ 4: RRT + exact `Rational` synthetic deflation for integer coefficients, quadratic finish, else Durand–Kerner + real-root certification. (2026-10-19; unsolve `solver.ts:216-265`)
- [x] **Taylor non-analytic error** — `taylor_series` returns `MathError::Domain` when any `f^(k)(a)` is non-finite (singularity at the expansion point) instead of `unwrap_or(0.0)` (`taylor.rs:30`). Regression tests `taylor_pole_errors_instead_of_silent_zero` (1/x at 0) and `taylor_log_at_zero_errors` (log x at 0); `taylor_analytic_far_from_pole_still_works` guards the valid case (1/x at 1). Error message names `f^k(a)` and the offending value. (2026-10-03)
- [x] **Lean 4 export** — `lean` module: `to_lean`/`to_lean_doc` render `Expr` as Lean 4 terms over `ℝ` (requires Mathlib). Exact ascribed literals (`(2 : ℝ)`), integral exponents via monoid `^ (ℕ)` (correct for negative bases — rpow has junk there), `Real.*` function mapping with `log(x,b)` → `Real.logb b x` arg swap, `floor`/`ceil` coerced `↑⌊·⌋`/`↑⌈·⌉`, unmapped functions (`gamma`, `erf`, …) as plain applications, NaN/±inf rejected as `Domain`. Parser-folded `pi`/`e` exported as exact decimal literals (silently upgrading to `Real.pi` would be unsound). `lean <expr>` REPL command prints a self-contained file (`import Mathlib` + `noncomputable def` with var binders). No import direction (parsing Lean syntax is out of scope). Validation: 28 unit tests (21-case golden corpus, wellformedness validator as compile-proxy — balanced parens/no NaN/inf/every numeral ascribed, π/e/τ exact-decimal soundness goldens, numeric edge cases, full mapped-function sweep) + 5 integration tests incl. real-pipeline exports (differentiate/apart/taylor/qsolve → Lean). `examples/lean_demo.rs` worked-example walkthrough. (2026-10-13)
- [x] **Like-term collection + coefficient collapse in simplify** — `2*x+3*x → 5*x`, `2*(2*x) → 4*x`, `x*2 → 2*x`, `x - x → 0`, `--x → x`; also fixes the simplify.rs doc-comment drift (claims like-term collection that didn't exist). New `preflatten` + `collect` pipeline; `split_coefficient(c, rest)` + `rebuild_term` helpers; tests for Monomial/multi-var/neg/distributing edge cases. Side-effect: canonicalize sort puts numeric lex-first so `1/(x+1)` now displays as `1/(1 + x)` — downstream tests in `tests/integration.rs::apart_repl_dispatch` and `src/apart.rs` doctest updated to accept both orderings. (2026-10-13)

### High Priority — cross-port pass #3 (2026-10-30)
(seeded from the unsolve cross-port pass #3 — see MEMORY.md "Cross-Port Learnings from unsolve (pass #3)")

- [x] **`for <var>` suffix on calculus/solve commands** — `diff`/`d/dx`/`integrate`/`solve`/`taylor`/`limit` accept a trailing `for <name>` naming the free variable of the operation. Fixes a **silent wrong answer**: `diff y*x^2 for y` used to return `x^2*for` (the words `for` and `y` were parsed as identifiers and multiplied in). New `strip_for_var` helper (`repl.rs:287-329`) strips the suffix *before* the existing trailing-token heuristics (single-letter var, numeric order), so `for y` is not misread as a variable or an order argument. Keyword matched case-insensitively; the name must be letters-only and the suffix must be final (`x^2 for x extra` is not a suffix). An explicit `for <var>` wins over the trailing single-letter heuristic in all five commands. 16 unit tests in `repl.rs::tests::for_var_*` / `strip_for_var_*`. (2026-10-30; unsolve `solve-problem.ts:119-124`, `solver.ts:405-440`)
- [x] **Integrate pulls out free (non-numeric) factors** — `∫ a·x dx = a·x²/2` for symbolic `a`, not just numeric constants. The `Mul` arm of `int_step` now tests *syntactic occurrence* (`mentions_var`) instead of `Expr::is_constant`, so any factor free of the integration variable is treated as a constant of integration. New `mentions_var` helper in `symbolic.rs`. 5 unit tests in `symbolic.rs::tests::integrate_pulls_free_factor_*` (both product positions, higher powers, function factors, and the guard that a variable-mentioning factor is *not* pulled out). (2026-10-30; unsolve `symbolic.ts:307-317`)
- [x] **`-(a - b)` flips to `b - a`** — negating a subtraction swaps its operands instead of leaving a `Neg` wrapping a `Sub`. Keeps solver answers readable: a linear isolate builds `neg(constant)`, so `solve 2*y + 3*x = 6 for y` reads `y = (6 - 3*x)/2` rather than `y = -(3*x - 6)/2`. Implemented in `preflatten`'s `Neg` arm (`simplify.rs:37-52`). 6 unit tests in `simplify.rs::tests::neg_of_sub_*` (operand flip, numeric case, no lingering `Neg` node, value preservation, double negation, self-cancel). (2026-10-30; unsolve `simplify.ts:218-219`)
- [x] **Exact `n*(e/n)` cancellation — already correct, verified not ported** — mathr's `collect` pass cancels `3*(x/3) → x` *structurally* (multiplies the coefficients) rather than by float pattern-matching, so it also reduces `(3*x)/3 → x` and `4*(x/2) → 2*x` — both of which unsolve leaves unreduced. Unsolve's rule compares raw floats (`simplify.ts:321,324`), a latent false-cancel hazard; mathr's coefficient-multiplication approach is strictly better. No code change; recorded to prevent a future "port" that would regress it. (2026-10-30)

### Medium Priority
- [x] **MCP (Model Context Protocol) server** — `mcp` module: JSON-RPC 2.0 over **stdio** (newline-delimited, protocol version 2025-06-18). Handshake (`initialize` → `notifications/initialized`), `tools/list`, `tools/call`, `ping`. 16 tools (`evaluate`, `differentiate`, `integrate`, `solve`, `nsolve`, `simplify`, `taylor`, `limit`, `poly-roots`, `number-theory`, `stats`, `matrix`, `minimize`, `apart`, `expand`) — each a thin wrapper that builds a command string and runs it through `repl::dispatch_str`, so the MCP surface cannot drift from the CLI. Hand-rolled JSON layer (no `serde` dep), matching the `serialize.rs` convention; exact rational results preserved (`1/2 + 3/4` → `5/4`). Tool *execution* failures return a successful response with `isError: true` (the MCP convention); only protocol problems use JSON-RPC error codes (`-32601` method not found, `-32602` invalid params). Launched with `mathr mcp`. 42 unit + 11 integration tests. (2026-10-30; unsolve `api/mcp.ts`)
- [x] **`d/dx <expr>` alias for `diff`** — `d/dx <expr> [order]` fixed wrt `x`; shared `parse_diff_tail` with `diff`. (2026-10-19)
- [x] **Selective exact function folding** — `sqrt`/`cbrt`/`log2`/`log10`/`abs` now fold **only when the result is provably exact**: `sqrt(16)`→`4` but `sqrt(2)` stays `sqrt(2)` (an irrational value has a shorter exact symbolic form, so collapsing it to a float discarded exactness). New `Rational::exact_sqrt`/`exact_cbrt`/`exact_log2`/`exact_log10`/`exact_factorial` helpers (exact square/cube roots of rationals, with an exact cross-check that `r² == original`), plus `integer_sqrt`/`integer_cbrt` (Newton + a float-error guard). The folding runs in a **new first pass** (`fold_exact_functions`) *before* `collect`, because `collect` turns `1/4` into the float `0.25`, after which the exactness of `sqrt(1/4)` can no longer be established. `try_eval_const` no longer folds exact-only functions at all, so the final pass cannot re-flatten them. 20 unit tests in `simplify.rs` + 10 in `rational.rs`. (2026-10-30; unsolve `simplify.ts:221-270`)
- [x] **Sign-scan all-roots numeric solve** — `solver::find_all_roots` + `ScanOptions`: samples an interval (default `[-30, 30]`, 2400 points), bisects every sign change, and returns **all** roots. New `nsolve <expr> [var] [lo hi]` REPL command (plus step-by-step output) — `nsolve sin(x)` finds `0, ±π, ±2π, …` where `solve` returns only the root nearest the guess. **Pole rejection** is the critical property: a sign flip across a pole looks exactly like a sign flip across a root, so a bisected bracket is accepted only when `f` actually lands near zero there — `nsolve 1/(x-1)` reports *no real roots* (the pole at x=1 is not a solution), and `nsolve tan(x)` reports only the true multiples of π, not the poles at ±π/2. Undefined points (NaN) break the scan rather than producing a spurious bracket. 10 unit tests in `solver.rs` + 9 integration tests. Also fixed `format_value(-0.0)` printing `-0`. (2026-10-30; unsolve `solver.ts:286-329`)
- [x] **Exact-rational quadratic path in qsolve (integer coeffs)** — `try_exact_rational_quadratic_roots`: near-integer `a,b,c`, perfect-square discriminant → `Rational` roots + `fmt_root` on `Div`; `6x^2-5x-1` → `1, -1/6`. Complex real part still via `snap_to_string` on f64 (`x^2+x+1` → `-1/2`). (2026-10-19; unsolve `solver.ts:380-451`)
- [x] **Exact-rational quadratic (general coeffs)** — the quadratic path now handles **rational** coefficients, not just integers, and keeps irrational roots as **exact radicals** instead of decimal approximations. `(1/2)x²−(3/4)x+1/8` → `3/4 + sqrt(5)/4` and `3/4 − sqrt(5)/4` (previously `1.30901699437`); `x²−2` → `sqrt(2)`; `x²−8` → `2*sqrt(2)`; `(2/3)x²+(1/3)x−1` → `1`, `−3/2`. New `exact_coeff` reads a coefficient as an exact rational (round-trip verified); the discriminant `b²−4ac` is computed in exact `Rational` arithmetic; `split_radical` pulls out the largest square *factor* (so `√8 → 2√2`, `√(5/16) → √5/4`) and `build_radical_root` writes the result as `p ± q·√d` with every coefficient kept rational. A perfect-square discriminant still yields plain rationals (the pre-existing integer path is untouched). 8 unit tests in `qsolve.rs` + 4 integration tests. Two existing tests were updated to the improved output (`irrational_discriminant_gives_decimals` → `..._gives_exact_radicals`, and the Lean-export pipeline now sees `Real.sqrt`). (2026-10-30; unsolve `solver.ts:380-451`)
- [x] **Exact cover-up fast path in apart** — when the denominator splits into **distinct rational** linear factors, the **cover-up rule** `Aᵢ = N(rᵢ)/D′(rᵢ)` (with `D′(rᵢ) = lead·Π_{j≠i}(rᵢ−rⱼ)`) gives the coefficients directly — no linear system, and **exact rationals**: `1/(x²−1)` → `1/2/(x−1) − 1/2/(x+1)` (was `0.5/(-1+x) + -(0.5/(1+x))`). The roots are found with the **rational-root theorem** in exact `Rational` arithmetic (candidates `p/q` with `p` dividing the constant term and `q` the leading coefficient, **both signs** of `p` tried). The path is declined — falling back to the numeric linear-system path — for repeated roots, complex roots, or irrational roots. `apart`/`apart_steps` no longer re-simplify the result, because `simplify` folds rationals to floats. 12 unit tests in `apart.rs` + 4 integration tests. The Lean exporter now sees exact rationals (`(1 : ℝ) / (x + (1 : ℝ))`) instead of decimals — more correct for Lean; one existing test was updated. (2026-10-30; unsolve `apart.ts:179-197`)
- [x] **Wire missing matrix commands** — `inv`, `transpose`, `trace`, `linsolve` (augmented `[A|b]`), full symmetric `eigen` command over the existing `symmetric_eig`. Symmetry check lives in the REPL layer (Mature Module Policy: `symmetric_eig` documents itself as unchecked, and changing it would alter library behavior beyond this scope). Numeric f64 (rational exact path deferred — `rational.rs` has no linear solve yet). Confirmed by the unsolve parity pass 2026-10-13. (unsolve `matrix.ts`, `solve-problem.ts:189-195`)
- [x] **Unicode input normalization** — `unicode` module: a pre-tokenizer that rewrites Unicode math into the ASCII the parser understands. Character substitution: `×`/`·`/`⋅`→`*`, `÷`/`∕`→`/`, `−`(U+2212)/en-dash/em-dash→`-`, `π`/`Π`→`pi`, `τ`→`tau`, `∞`→`inf`, `±`→`+`. Superscript runs → explicit powers: `x²`→`x^(2)`, `x⁻¹`→`x^(-1)`, `x¹²`→`x^(12)`, `xⁿ`→`x^(n)`, `x½`→`x^(1/2)`. The root signs need *argument delimiting*, not plain substitution — `√2`→`sqrt(2)` (not `sqrt2`, which reads as one identifier), so the root branch captures its argument and normalizes it separately. Unmapped glyphs pass through unchanged so the parser reports the real error. 17 unit tests in `unicode.rs` + 13 parser tests + 4 integration tests. (2026-10-30; unsolve `parser.ts:25-70`)
- [x] **Missing LaTeX forms** — `\dfrac`/`\tfrac` are now aliases of `\frac`; `\sqrt[n]{x}` is an n-th root (`8^(1/3)`); `\sin^2(x)` is a function power (`sin(x)^2`), with the exponent read *before* the argument via a new `read_tex_superscript` helper. 9 parser tests + 3 integration tests. (2026-10-30; unsolve `parser.ts:113-168`)
- [x] **`from a to b` suffix parsing with expression bounds** — `int`, `romberg` and `minimize` now share one bound parser (`repl::split_bounds`) that accepts **either** the historical positional form (`<expr> <a> <b>`) **or** an explicit `from <a> to <b>` suffix. The suffix's bounds are **full expressions** and may contain spaces (`from pi/2 to 2*pi`). This fixes an inconsistency: `romberg sin(x) 0 pi` used to fail with "bad b" while `int sin(x) 0 pi` worked — now both accept symbolic bounds. The positional form still takes the **last two tokens** as the bounds, so the trailing-token ambiguity class stays retired (`romberg x^2 0 1 12` still integrates `x^2*0` over `[1, 12]`). A malformed suffix (`from 0` with no `to`) is a clear error instead of silently treating "from" as part of the expression. 13 unit tests in `repl.rs` + 5 integration tests. The MCP `integrate` tool gained a definite-integral mode (`from`/`to`) and `minimize` now accepts expression bounds. (2026-10-30; unsolve `solve-problem.ts:115-133`)
- [x] **Numeric domain probes with named failure points** — `integrate_adaptive` and `integrate_romberg` now check the integrand at the endpoints and the midpoint and **name the offending point** instead of returning `NaN`; `minimize` checks `f(a)`/`f(b)` instead of silently mapping an undefined value to `+∞`. This fixed a genuine **hang**: `int sqrt(x) -1 1` recursed to the depth limit without converging (the error estimate is `NaN`, so the convergence test never passes) — the branching recursion is `2^max_depth`, i.e. up to 2³⁰ calls. Now `int sqrt(x) -1 1` reports "the integrand is not defined at x = -1". 6 unit tests in `calculus.rs` + 4 in `repl.rs` + 5 integration tests. (2026-10-30; unsolve `optim.ts:33-35`, `numeric.ts:51-57`)
- [x] **BigInt auto-upgrade for `is-prime`/`factor`** — both now fall back to `BigInt` when the argument does not fit in `u64` (the same pattern `fib`/`binom`/`fact` already used), so `is-prime 100000000000000000000` → `false` and `factor 100000000000000000000` → `2^20 · 5^20` instead of "could not parse as integer". Verified against known values (10^20 = 2^20·5^20; 2^127−1 is a Mersenne prime). The `u64` path is kept as the fast path. 4 unit tests + 2 integration tests. (2026-10-30; unsolve `solve-problem.ts:458-464`)
- [x] **Number-base conversion** — new `numbase` module: `to_base`/`from_base` for bases 2–36, **exact** via `BigInt` (a 40-digit decimal converts to hex and back without loss). Base names (`bin`/`oct`/`dec`/`hex`) or a number. New `base <n> <target>` REPL command (`base 255 hex` → `FF_16`, `base 42 36` → `16_36`) plus a `base <n> from <base> to <base>` form for converting *between* bases (`base FF from 16 to 10` → `255_10`). Negative values and both letter cases are handled; invalid digits are rejected. 11 unit tests + 4 repl tests + 2 integration tests. (2026-10-30; qalculate pass)
- [x] **Implicit differentiation** — `symbolic::implicit_diff` / `idiff_steps_str`: `dy/dx = −F_x/F_y` for `F(x,y)=0`; `idiff <expr> [= rhs] [<indep> <dep>]` REPL with step-by-step (F, ∂F/∂x, ∂F/∂y, result); infers two free vars (prefer `x` independent). (2026-10-03)
- [x] **Dot/cross product REPL commands** — `vector` module (`dot`, `cross`/`cross3`); REPL `dot a... | b...`, `cross a1 a2 a3 | b1 b2 b3`. (2026-10-19)
- [x] **Higher-order `diff` order argument** — `diff <expr> [var] [n]` / `d/dx <expr> [n]` via `symbolic::differentiate_n` (cap 64). (2026-10-19)
- [x] **Complex-number expression evaluation** — `ceval` module: `eval_complex`/`eval_complex_str` over the Expr AST (`Var("i")` = imaginary unit, no parser changes), principal-branch ln/sqrt/pow/exp/trig/inverse-trig/hyperbolic, exact integer powers, `cval <expr> [with <var>=<val>,...]` REPL command, pretty `format_complex` (snaps IEEE −0.0 and noise so `exp(i*pi)` renders `-1`). Validated by Euler identity, principal branches, identity cross-checks (`sin(2i) = i·sinh 2`), and `asin(2)` closed form. (2026-09-15)

### Competitive intelligence notes (2026-10-03)
- Source: Qalculate! / libqalculate v5.12 manual + release notes (web2md MCP unavailable — used WebSearch/fetch fallback).
- Gaps worth seeding (beyond existing Medium items): higher-order `diff` order arg; `dsolve` for simple 1st-order ODEs; LaTeX *output* (`to latex`); entrywise logical/bitwise on vectors; richer `multisolve`.
- Already covered or deferred: Σ/∏, interval arithmetic, number bases (Medium), units/currency (out of scope for pure math crate).
- **2026-10-03 — unsolve.org** (TypeScript port of mathr's symbolic core, newer than some Rust code). Reviewed 19 engine modules + 10 test files. Seeded High Priority (linear-substitution integration, continued-fraction snap, expression-level cubic solve, Taylor non-analytic error, like-term collection), 9 Medium items (exact function folding, sign-scan all-roots, exact-rational quadratic, exact cover-up, missing matrix commands, Unicode/LaTeX forms, `from a to b` suffix, domain probes, BigInt auto-upgrade) and 7 Low items (Durand–Kerner polish, degenerate-equation semantics, separate L/R limit probes, rational overflow, matrix det pivot, bare-equation auto-solve, `-x^2` precedence, symbolic-coefficient univariate expand). Full patterns in MEMORY.md "Cross-Port Learnings from unsolve".
- Cross-port gaps unsolve has that mathr *doesn't* need to copy: TS structural equality is string-based; mathr already has `canonicalize`/`equals`. TS is univariate-only; mathr has multivariate `poly.rs` + general `apart.rs` (repeated + irreducible quadratics). mathr's VAS isolation is strictly more rigorous than TS's numeric divisor hunt — wire VAS into the new expression-level cubic path instead of replicating TS's.
- **2026-10-13 — unsolve parity pass #2** (reviewing `api/mcp.ts` MCP surface + `api/lib/engine/solve-problem.ts` dispatch). **Confirmed + shipped this pass**: matrix commands (`inv`/`transpose`/`trace`/`linsolve`/`eigen`) — the highest-value unsolve gap; `repl.rs` now mirrors unsolve's `matrix` tool. **New gap surfaced**: **MCP server**. unsolve exposes 14 MCP tools over its engine (`evaluate`, `differentiate`, `integrate`, `solve`, `simplify`, `taylor`, `limit`, `poly-roots`, `number-theory`, `stats`, `matrix`, `minimize`, `apart`) — mathr has the HTTP notebook server but no MCP endpoint. MCP is the standard interop layer for AI agents and is the highest-impact remaining "external surface" gap. mathr's `repl::dispatch_str(&str, &Context) -> Result<String>` is a natural one-call mapping for almost every tool. Trivial secondary gap: `d/dx <expr>` alias for `diff` (unsolve has it; mathr's `diff` already exists). See MEMORY.md "Cross-Port Learnings from unsolve (2026-10-13)".

### Competitive intelligence notes (2026-11-01)
- Sources: current docs.rs pages for **Symbolica 3.0.2**, **rustfft 6.4.1**, **argmin 0.11.0**, **nalgebra 0.35.0**, and **statrs 0.19.1**. The required `web2md` MCP server is not configured in this environment, so the pass used the built-in search/fetch fallback.
- **Seeded candidate — arbitrary-size FFT**: mathr only has radix-2 Cooley–Tukey; rustfft's planner covers arbitrary and prime sizes in O(n log n). A mixed-radix/Bluestein path is a clearer user gap than GPU FFT.
- **Seeded candidate — quasi-Newton optimization**: mathr has golden section and Nelder–Mead; argmin exposes BFGS/L-BFGS, Newton-CG, trust-region, and line searches. mathr's symbolic/autodiff gradients make BFGS a natural next optimizer without a new dependency.
- **Seeded candidate — Brent root finding**: argmin includes Brent's bracketed method; mathr has robust bisection and faster secant/Newton but no safeguarded interpolation method. Brent is a compact, high-value solver addition.
- **Seeded candidate — second-order autodiff**: Symbolica exposes hyperdual derivatives; mathr's `Dual` covers values/first derivatives, gradients, and Jacobians but not Hessians. A small hyperdual type could provide exact second derivatives and Hessian-vector products.
- **Deferred**: Symbolica's large-expression/JIT/code-generation surface is out of scope for mathr's small-footprint goal; nalgebra's geometry transforms and statrs's full distribution catalogue exceed the crate's current focus. GPU FFT, 3D plotting, and animated plots remain parked.

### Low Priority
- [ ] **Arbitrary-size FFT** — mixed-radix or Bluestein path for non-power-of-two lengths (seeded by rustfft 6.4.1 research).
- [ ] **BFGS quasi-Newton optimization** — gradient-based N-D minimization using symbolic/autodiff gradients and a line search (seeded by argmin 0.11.0 research).
- [ ] **Brent root finding** — safeguarded interpolation combining bisection robustness with secant/inverse-quadratic speed (seeded by argmin 0.11.0 research).
- [ ] **Hyperdual second derivatives** — Hessian and Hessian-vector products over the Expr AST (seeded by Symbolica 3.0.2 research).
- [x] **Durand–Kerner polish** — leading coefficients below `1e-14 × max|coeff|` are stripped (scale-relative, so uniformly tiny polynomials are untouched); near-real roots are snapped onto the axis, near-integer components are snapped to integers, and repeated roots in the real-root API are deduplicated. `poly-roots 1 -2 1` now returns one clean `x = 1 (f = 0)` instead of two approximate roots, and `poly-roots 0 1 -5 6` / `poly-roots 1e-20 1 -5 6` both return `2, 3` rather than erroring or adding a spurious huge root. The complex-root API deliberately **retains multiplicity** because `apart` reconstructs repeated factors by counting roots. 5 unit tests + 2 integration tests. (2026-11-01; unsolve `poly.ts:48-100`)
- [x] **Degenerate-equation semantics in qsolve** — the two degenerate cases are now distinguished: an **identity** (`x − x = 0`) reports "every value is a solution (the equation is an identity)", a **nonzero constant** reports "no solution (the equation is the constant 5 ≠ 0)". Previously both produced the same generic "no variable dependence" error. 3 unit tests + 2 integration tests. (2026-10-30; unsolve `solver.ts:527-541`)
- [x] **Duplicate-root dedup in qsolve** — a repeated root is listed **once**: `qsolve (x-1)^2` → `x = 1` (was `x = 1, x = 1`, which reads as two distinct solutions). New `dedupe_roots` merges roots whose numeric values agree to a relative 1e-9; roots that are not plain numbers (exact radicals, complex expressions) are never merged. Distinct roots and complex-conjugate pairs are unaffected. 1 existing test updated to the improved output. (2026-10-30; unsolve `poly.ts:48-100`)
- [x] **Separate left/right limit probes + growth-trend pole detection** — the probe now collects the two sides **separately** (the old code interleaved `point−h, point+h` and then split the array in half, so "left" was really the coarse samples and "right" the fine ones). Pole detection uses a **growth trend** over the last three samples (each >5× the previous) instead of a single `>1e9` magnitude threshold — the threshold missed small-magnitude poles like `1/(x−1)` (which only reaches ~100 at the closest sample) and misfired on merely-large finite values. 4 unit tests in `limit.rs` + 2 integration tests. (2026-10-31; unsolve `limit.ts:221-263`)
- [x] **The `at` keyword for limits** — `limit sin(x)/x at 0` used to read the literal string "at" as the *variable*, so it computed the limit with respect to a variable named "at" and wrongly reported "does not exist". The variable is now inferred from the expression when `at` is present; an explicit `for <var>` still wins. 4 unit tests in `repl.rs` + 1 integration test.
- [x] **Rational overflow soft-fail** — `reduce_i128` narrowed the `i128` intermediates with an unchecked `as i64`, so `i64::MAX × 2` **silently wrapped to −2** — a wrong answer, not an error. The narrowing is now **checked** and **saturates** to the largest representable magnitude (with the correct sign), preserving the `Rational` invariant (a reduced fraction that fits in `i64`). `parse_rational` now distinguishes "not a number" from "a valid number too large for 64-bit", and the latter points at the `big` command. 4 unit tests + 2 integration tests. (2026-10-30; unsolve `rational.ts:21,38`)
- [x] **Matrix det scale-relative pivot threshold** — the pivot test used an **absolute** threshold (`best < 1e-14 → Ok(0.0)`), which reported a tiny-scaled but nonsingular matrix as singular: `det [[1e-15,0],[0,1]]` returned `0` instead of `1e-15`. The threshold is now **relative to the matrix's scale** (largest absolute entry × EPSILON × n), so tiny-scaled matrices keep their determinant while genuinely singular ones still report `0`. A related fix in `repl::format_value`: a value below the 10-decimal format's resolution (`< 1e-10`) now prints in scientific notation instead of rounding away to `0`. 3 unit tests in `matrix.rs` + 3 integration tests. (2026-10-31)
- [x] **Bare-equation auto-solve + friendlier `=` errors** — a top-level `=` now routes to the solver instead of failing with "unexpected trailing characters at byte 8": `x^2 - 4 = 0` → `x = 2, x = -2`. **Polynomial** equations go through the symbolic solver (`qsolve`), so complex roots (`x^2+1=0` → `x = i, x = -i`) and exact radicals (`x^2-2=0` → `x = sqrt(2)`) come out exactly; **non-polynomial** equations fall back to the numeric solver (`exp(x) = 5` → `root ≈ 1.609`). The `=` must be at the *top level* — a nested `(a = b)` or a comparison `!=`/`==` is not an equation. A leftover `=` in another command now names the character and points at `solve`: "unexpected '=' — for equations use: solve <lhs> = <rhs>". 6 unit tests in `repl.rs` + 4 integration tests. (2026-10-31; unsolve `parser.ts:466-483`)
- [x] **LaTeX output** — new `latex <expr>` REPL command exposing the existing `Expr::to_tex`: `latex x^2/2` → `\frac{x^{2}}{2}`, `latex sin(x)` → `\sin\left(x\right)`. The counterpart of mathr's LaTeX *input* support; the emitted TeX **round-trips** (it parses back to an equal expression) for fractions, powers, functions and radicals. 2 unit tests + 2 integration tests. (2026-10-31; qalculate pass)
- [x] **`-x^2` precedence** — prefix unary minus now applies to the complete power/postfix expression: `-x^2` parses as `-(x^2)` (conventional mathematics) rather than `(-x)^2`; explicit `(-x)^2` retains the old positive-square meaning. `parse_factor` collects all prefix signs before the atom and applies their parity after factorial/power expansion, so signed exponents (`2^-3`, `2^-x^2`) and right-associative powers keep working. This was an **approved breaking parser change**; the stale polynomial-expansion and autodiff assumptions were updated. 2 parser unit tests + 1 REPL integration test. (2026-11-01)
- [x] **Univariate expansion with symbolic coefficients** — new `poly::expand_in` / `to_uni_poly`: expand in **one** variable, treating every other sub-expression as an **opaque exact coefficient** (an `Expr`, not a float). The multivariate `expand` can only distribute when both sides are polynomials, so `(x + sin(y))^2` stayed unexpanded; now it gives `x^2 + 2·sin(y)·x + sin(y)·sin(y)`, and `(x + 2·y)^3` gives `x^3 + 6·y·x^2 + 12·y^2·x + 8·y^3`. The "termsOf" style the CASs use. New `expand-in <expr> for <var>` REPL command (plus step-by-step output); the variable is required. The coefficients are simplified and expanded on the way out, so repeated factors fold into powers. 7 unit tests in `poly.rs` + 4 integration tests, including numeric agreement with the original at sample points. (2026-10-31; unsolve `poly.ts:132-210`)
- [x] **Exact symbolic linear system solve** (qalculate's `multisolve`) — new `linsolve` module: solves 2–4 linear equations **exactly** by rational Gaussian elimination (partial pivoting + back-substitution over `Rational`). The equations are given as text — `2x + y = 5; x + 3y = 10` — either with an explicit `=` or with the right-hand side implicit (`2x + y - 5`). New `msolve <eq>; <eq>; …` REPL command; `linsolve` keeps its existing numeric augmented-matrix meaning. The answer stays **exact** (`(1/2)x + y = 1; x + (1/3)y = 1` → `x = 4/5, y = 3/5`, not decimals). A singular system is reported as such, and a non-linear equation is rejected rather than mis-solved. 10 unit tests in `linsolve.rs` + 6 integration tests. (2026-10-31; qalculate pass)
- [x] **Simple symbolic ODE solve (`dsolve`)** — new `dsolve` module: solves first-order ODEs `y' = rhs` symbolically by pattern recognition, producing `y` with the constant `C`. Three patterns: **exact derivative** (`y' = 2x` → `y = C + x²`), **separable** (`y' = x·y` → `y = C·e^(x²/2)`, `y' = y·sin(x)` → `y = C·e^(−cos x)`), and **linear** via integrating factor (`y' = −y` → `y = C/e^x`). Every solution is **verified** by differentiating it and comparing numerically along the solution curve (the RHS's `y` is evaluated from the solution itself, not an arbitrary point). Pattern order matters: `x·y` is both separable *and* affine in `y`, so separable is tried first. A failed integration *declines* the pattern (`Ok(None)`) instead of erroring — `y' = x − 2xy` needs `∫e^(−x²)·x`, which has no closed form here, so it correctly reports "cannot solve". New `dsolve y' = <rhs>` / `dsolve dy/dx = <rhs>` REPL command with step-by-step output. 8 unit tests, all solutions verified. (2026-10-31; qalculate pass)
- [x] **Symbolic quadratic solving** — `qsolve` module: `solve_symbolic(&Expr, var)` extracts coefficients via `to_poly` after simplification, classifies by discriminant (perfect-square → exact integer roots, zero → double root, negative → complex-conjugate roots as Exprs with the `i` literal, evaluable by `ceval`), 12-sig snapping otherwise. `qsolve <expr> [= rhs]` REPL command with step-by-step output (equation → coefficients → discriminant → roots); exactly one inferred variable (avoids the trailing-token ambiguity class). Cubic+ via radicals out of scope. (2026-09-15)
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