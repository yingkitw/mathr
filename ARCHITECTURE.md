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
                    solver.rs   (root finding, 1-D and systems)
                    taylor.rs   (uses symbolic + eval)

complex.rs ──▶ fft.rs (Cooley–Tukey, convolution, cross-correlation, windows)

matrix.rs       (standalone, f64 linear algebra: LU, Cholesky, SVD, power iteration)
stats.rs        (standalone, descriptive statistics)
numtheory.rs    (standalone, integer number theory + Miller–Rabin, CRT,
                 Jacobi symbol, continued fractions, Diophantine, discrete log)
ode.rs          (standalone, numerical ODE integration)
interpolate.rs  (standalone, Lagrange, Newton, cubic spline, Chebyshev,
                 Legendre, Gauss–Legendre)
special.rs      (standalone, Gamma, Beta, erf, sinc, incomplete gamma,
                 Bessel J_0/J_1/J_n)

expr.rs ──canonicalize/equals──▶ expr.rs (canonical form comparison)

plot.rs ──uses──▶ eval.rs, expr.rs
repl.rs ──uses──▶ all of the above
error.rs ──used by──▶ all modules
```

## Data Flow

1. **Parse**: `Parser::parse(str) → Expr` — recursive-descent, handles precedence, implicit multiplication, function calls
2. **Evaluate**: `eval(&Expr, &Context) → f64` — tree-walking with variable/function context
3. **Symbolic**: `differentiate(&Expr, var) → Expr` and `integrate(&Expr, var) → Expr` — apply calculus rules, then simplify
4. **Numerical**: `calculus::derivative/integrate_trap/integrate_simpson/integrate_adaptive/integrate_romberg` — operate on closures `F: Fn(f64) → f64`
5. **Solve**: `solver::bisect/newton/secant/polynomial_roots` (single-variable) and `solver::newton_system` (`n`-dimensional with central-difference Jacobian)
6. **FFT**: `fft::fft/ifft/rfft` — operates on `Vec<Complex<f64>>` or `Vec<f64>`
7. **Matrix**: `Matrix` — row-major `Vec<f64>`. Gaussian elimination for det/solve/inverse; **LU** with partial pivoting, **Cholesky** `A = L·Lᵀ` for SPD, **SVD** `A = U·Σ·Vᵀ` via Jacobi rotations, **power iteration** for the dominant eigenpair
8. **Stats**: Functions on `&[f64]` slices
9. **Number theory**: Functions on `u64` integers (gcd, sieve, totient, Miller–Rabin, CRT, modular inverse, **discrete logarithm**) and on `i64` (extended GCD, **Jacobi symbol**, **continued fractions**, **linear Diophantine solver**)
10. **ODE**: `ode::euler/rk4/rkf45` — operate on closures `F: Fn(f64, f64) → f64`
11. **Taylor**: `taylor::taylor_series` — uses `symbolic::differentiate` + `eval::eval` to compute coefficients
12. **Symbolic integration**: `symbolic::integrate` — pattern-matches common elementary rules (polynomial, exp, ln, sin/cos/tan/sec, atan, asin)
13. **Interpolate**: `interpolate::lagrange_interp/newton_interp/CubicSpline/chebyshev_*/legendre_*/gauss_legendre` — point-wise polynomial, smooth C²-continuous cubic, Chebyshev series, Legendre polynomials, Gauss–Legendre quadrature
14. **Special**: `special::gamma/erf/sinc/bessel_j0/j1/jn` — Lanczos, continued-fraction, Maclaurin-series, and asymptotic-form approximations
15. **Plot**: `plot::plot_function/multi/scatter` — evaluates `Expr` over a range, renders PNG via `plotters`

## Design Decisions

- **No external math crates**: FFT, complex numbers, matrix ops, stats, special functions, and number theory are all implemented from scratch. Only `plotters` is used for rendering.
- **Expr as the universal AST**: One enum serves parsing, evaluation, symbolic differentiation, symbolic integration, simplification, and Taylor series.
- **Closures for numerical methods**: Solvers, integrators, and ODE steppers accept `Fn` closures, making them composable with the evaluator.
- **Inline tests**: Each module has `#[cfg(test)] mod tests` — no separate test files. **209 tests** total (plus 53 integration tests).
- **Error handling**: `MathError` (via `thiserror`) in the library, `anyhow` in the binary.
- **Prelude**: `mathr::prelude` re-exports the most common types for ergonomic `use mathr::prelude::*;`.
- **Numerical recipes**: Where possible (Gamma, sinc, Bessel, incomplete gamma) we use well-tested A&S or NR polynomial approximations; otherwise we use direct series / closed-form recurrences.
- **Matrix decomposition family**: Gaussian elimination for inverse/determinant; LU with partial pivoting for fast solves; Cholesky for symmetric positive-definite matrices; SVD via Jacobi rotations for rank/reconditioning; power iteration for the dominant eigenpair.

## Deployment

The crate produces both a library (`mathr`) and a binary (`mathr`). The binary is a thin CLI wrapper around library functions. No network, no state, no config files — just math.