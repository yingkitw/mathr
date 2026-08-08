# Architecture

## Module Relationships

```
parser.rs ──parse──▶ expr.rs (Expr AST)
                        │
            ┌───────────┼───────────┐
            ▼           ▼           ▼
        eval.rs    symbolic.rs   simplify.rs
            │           │           │
            ▼           ▼           ▼
        (f64)    diff() / integrate()
                        │           → Expr
                        ▼
                    calculus.rs (numerical)
                    solver.rs   (root finding, 1-D and systems, polynomial root isolation)
                    taylor.rs   (uses symbolic + eval)
                    laurent.rs  (uses symbolic + eval, expansion around poles)

complex.rs ──▶ fft.rs (Cooley–Tukey, convolution, cross-correlation, windows)

matrix.rs       (standalone, f64 linear algebra: LU, Cholesky, SVD, power iteration,
                 symmetric eigenvalue decomposition via QR algorithm,
                 Hessenberg + real Schur decomposition)
stats.rs        (standalone, descriptive statistics)
numtheory.rs    (standalone, integer number theory + Miller–Rabin, CRT,
                 Jacobi symbol, continued fractions, Diophantine, discrete log)
rational.rs     (standalone, exact rational arithmetic with GCD reduction)
notebook.rs     (standalone, .mnb notebook format, JSON parse/serialize, cell eval,
                 cell types Math/Text, cell reordering, shared context via &mut Context)
server.rs       (standalone, minimal HTTP server for web notebook UI,
                 shared context across cells, /api/reset, /api/context)
ode.rs          (standalone, numerical ODE integration)
interpolate.rs  (standalone, Lagrange, Newton, cubic spline, Chebyshev,
                 Legendre, Gauss–Legendre)
special.rs      (standalone, Gamma, Beta, erf, sinc, incomplete gamma,
                 Bessel J_0/J_1/J_n)
fastmath.rs     (standalone, Chebyshev-based fast approximations of sin, cos,
                 tan, exp, log, sqrt, pow with argument reduction)
bigint.rs       (standalone, arbitrary-precision integers via num-bigint:
                 Miller–Rabin primality, Pollard's rho factorization, GCD,
                 factorial, Fibonacci fast doubling, binomial, mod_pow, totient.
                 REPL fact/fib/binom auto-upgrade on u64 overflow; `big` command
                 for inputs > u64::MAX)
autodiff.rs     (standalone, automatic differentiation via dual numbers:
                 Dual type with arithmetic, elementary functions, Expr eval,
                 derivative, gradient, jacobian)
mathml.rs        (standalone, W3C Presentation MathML export/import:
                 Expr→MathML and MathML→Expr, supporting mn/mi/mo/mrow/
                 mfrac/msup/msub/msqrt/mroot/mfenced/mstyle/mtext)

expr.rs ──canonicalize/equals──▶ expr.rs (canonical form comparison)

plot.rs ──uses──▶ eval.rs, expr.rs
repl.rs ──uses──▶ all of the above
error.rs ──used by──▶ all modules
```

## Data Flow

1. **Parse**: `Parser::parse(str) → Expr` — recursive-descent, handles precedence, implicit multiplication, function calls
2. **Evaluate**: `eval(&Expr, &Context) → f64` — tree-walking with variable/function context
3. **Symbolic**: `differentiate(&Expr, var) → Expr` (partial derivative w.r.t. one variable), `gradient(&Expr) → Vec<(String, Expr)>` (all partials), and `integrate(&Expr, var) → Expr` — apply calculus rules, then simplify
4. **Numerical**: `calculus::derivative/integrate_trap/integrate_simpson/integrate_adaptive/integrate_romberg` — operate on closures `F: Fn(f64) → f64`. **Fourier series** via `calculus::fourier_series/fourier_eval` — computes coefficients using Simpson's rule integration. **Monte Carlo** via `calculus::monte_carlo_integrate_1d/monte_carlo_integrate_nd` — reproducible LCG-based sampling with standard error estimates.
5. **Solve**: `solver::bisect/newton/secant/polynomial_roots` (single-variable), `solver::newton_system` (`n`-dimensional with central-difference Jacobian), and `solver::isolate_real_roots` (VAS root isolation with i128 exact arithmetic)
6. **FFT**: `fft::fft/ifft/rfft` — operates on `Vec<Complex<f64>>` or `Vec<f64>`
7. **Matrix**: `Matrix` — row-major `Vec<f64>`. Gaussian elimination for det/solve/inverse; **LU** with partial pivoting, **Cholesky** `A = L·Lᵀ` for SPD, **SVD** `A = U·Σ·Vᵀ` via Jacobi rotations, **power iteration** for the dominant eigenpair, **symmetric eigenvalue decomposition** via Householder tridiagonalisation + Wilkinson-shift QR iteration, **Hessenberg decomposition** via Householder reflections, **real Schur decomposition** via shifted QR on Hessenberg form, **Hilbert matrix** construction, **Tikhonov regularised solve** `(AᵀA + λI)x = Aᵀb` for ill-conditioned and rectangular systems
8. **Stats**: Functions on `&[f64]` slices. **Stochastic primitives**: `Rng` (reproducible LCG), `normal`/`exponential` sampling (Box–Muller, inverse-CDF), `normal_pdf`/`normal_cdf`/`exp_pdf`/`exp_cdf`, `moments` (skewness, excess kurtosis), `cumulants`
9. **Number theory**: Functions on `u64` integers (gcd, sieve, totient, Miller–Rabin, CRT, modular inverse, **discrete logarithm**) and on `i64` (extended GCD, **Jacobi symbol**, **continued fractions**, **linear Diophantine solver**)
10. **ODE**: `ode::euler/rk4/rkf45` — operate on closures `F: Fn(f64, f64) → f64`
11. **Taylor**: `taylor::taylor_series` — uses `symbolic::differentiate` + `eval::eval` to compute coefficients
12. **Laurent**: `laurent::laurent_series` — expands `g(x) = (x-a)^k · f(x)` via Taylor, divides by `(x-a)^k` to get principal + analytic parts
13. **Rational**: `rational::Rational` — exact `i64/i64` with GCD reduction, `i128` intermediate arithmetic, parsing from `"n/d"` or decimal strings. `rational::eval_rational` walks `Expr` AST and returns exact `Rational` when all leaves are integers (returns `None` for functions, variables, or non-integer constants)
14. **Notebook**: `notebook::Notebook` — `.mnb` JSON format with cells of TeX/math expressions; `eval_cell`/`eval_all` dispatch through `repl::dispatch_with_ctx` (mutating context for `let`/`fn` persistence); cells have types (`Math`/`Text`) and support reordering/duplication
15. **Web server**: `server::NotebookServer` — minimal HTTP server on `std::net::TcpListener`; serves single-page web UI at `GET /`, REST API at `/api/eval` (with step-by-step `steps` array, shared context, inline plot images as base64), `/api/notebook`, `/api/reset`, `/api/context`
16. **Step-by-step solving**: `repl::dispatch_steps` — returns `Vec<String>` of intermediate steps for `diff`, `solve`, `taylor`, `integrate`, `simplify`, `rat`, `laurent`, and plain eval; tries `rational::eval_rational` first for exact fraction results; falls back to `simplify` when evaluation fails on unbound variables; delegates all other REPL commands (`int`, `romberg`, `fft`, `det`, etc.) to `dispatch_inner`
17. **Symbolic integration**: `symbolic::integrate` — pattern-matches common elementary rules (polynomial, exp, ln, sin/cos/tan/sec, atan, asin)
18. **Gradient**: `symbolic::gradient` — collects all free variables and computes partial derivatives for each
19. **Interpolate**: `interpolate::lagrange_interp/newton_interp/CubicSpline/chebyshev_*/legendre_*/gauss_legendre` — point-wise polynomial, smooth C²-continuous cubic, Chebyshev series, Legendre polynomials, Gauss–Legendre quadrature
20. **Special**: `special::gamma/erf/sinc/bessel_j0/j1/jn` — Lanczos, continued-fraction, Maclaurin-series, and asymptotic-form approximations
21. **Plot**: `plot::plot_function/multi/scatter` — evaluates `Expr` over a range, renders PNG via `plotters`

## Design Decisions

- **No external math crates**: FFT, complex numbers, matrix ops, stats, special functions, and number theory are all implemented from scratch. Only `plotters` is used for rendering.
- **Expr as the universal AST**: One enum serves parsing, evaluation, symbolic differentiation, symbolic integration, simplification, and Taylor series.
- **Closures for numerical methods**: Solvers, integrators, and ODE steppers accept `Fn` closures, making them composable with the evaluator.
- **Inline tests**: Each module has `#[cfg(test)] mod tests` — no separate test files. **341 tests** total (plus 112 integration tests).
- **Error handling**: `MathError` (via `thiserror`) in the library, `anyhow` in the binary.
- **Prelude**: `mathr::prelude` re-exports the most common types for ergonomic `use mathr::prelude::*;`.
- **Numerical recipes**: Where possible (Gamma, sinc, Bessel, incomplete gamma) we use well-tested A&S or NR polynomial approximations; otherwise we use direct series / closed-form recurrences.
- **Matrix decomposition family**: Gaussian elimination for inverse/determinant; LU with partial pivoting for fast solves; Cholesky for symmetric positive-definite matrices; SVD via Jacobi rotations for rank/reconditioning; power iteration for the dominant eigenpair; symmetric eigenvalue decomposition via Householder tridiagonalisation + Wilkinson-shift QR iteration; Hessenberg decomposition via Householder reflections; real Schur decomposition via shifted QR on Hessenberg form with 2×2 block handling for complex conjugate eigenvalue pairs.

## Deployment

The crate produces both a library (`mathr`) and a binary (`mathr`). The binary is a thin CLI wrapper around library functions. The `notebook` subcommand starts a minimal HTTP server serving a Jupyter-like web UI with KaTeX math rendering, step-by-step solving, exact fraction arithmetic, shared context across cells (variables/functions persist), cell types (math/text with Markdown rendering), inline plots (base64 PNG), execution counters, and cell management (reorder, duplicate, toggle type).