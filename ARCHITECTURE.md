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
        (f64)    diff() → Expr   simplify() → Expr
                        │
                        ▼
                    calculus.rs (numerical)
                    solver.rs   (root finding)
                    taylor.rs   (uses symbolic + eval)

complex.rs ──▶ fft.rs (Cooley–Tukey, convolution, cross-correlation, windows)

matrix.rs       (standalone, f64 linear algebra)
stats.rs        (standalone, descriptive statistics)
numtheory.rs    (standalone, integer number theory + Miller–Rabin, CRT)
ode.rs          (standalone, numerical ODE integration)
interpolate.rs  (standalone, Lagrange + Newton interpolation)
special.rs      (standalone, Gamma, Beta, erf, sinc, incomplete gamma)

expr.rs ──canonicalize/equals──▶ expr.rs (canonical form comparison)

plot.rs ──uses──▶ eval.rs, expr.rs
repl.rs ──uses──▶ all of the above
error.rs ──used by──▶ all modules
```

## Data Flow

1. **Parse**: `Parser::parse(str) → Expr` — recursive-descent, handles precedence, implicit multiplication, function calls
2. **Evaluate**: `eval(&Expr, &Context) → f64` — tree-walking with variable/function context
3. **Symbolic**: `differentiate(&Expr, var) → Expr` — applies calculus rules, chain rule, then simplifies
4. **Numerical**: `calculus::derivative/integrate_*` — operates on closures `F: Fn(f64) → f64`
5. **Solve**: `solver::bisect/newton/secant/polynomial_roots` — operates on closures
6. **FFT**: `fft::fft/ifft/rfft` — operates on `Vec<Complex<f64>>` or `Vec<f64>`
7. **Matrix**: `Matrix` — row-major `Vec<f64>`, Gaussian elimination for det/solve/inverse
8. **Stats**: Functions on `&[f64]` slices
9. **Number theory**: Functions on `u64` integers
10. **ODE**: `ode::euler/rk4/rkf45` — operate on closures `F: Fn(f64, f64) → f64`
11. **Taylor**: `taylor::taylor_series` — uses `symbolic::differentiate` + `eval::eval` to compute coefficients
12. **Plot**: `plot::plot_function/multi/scatter` — evaluates `Expr` over a range, renders PNG via `plotters`

## Design Decisions

- **No external math crates**: FFT, complex numbers, matrix ops, stats, and number theory are all implemented from scratch. Only `plotters` is used for rendering.
- **Expr as the universal AST**: One enum serves parsing, evaluation, symbolic differentiation, simplification, and Taylor series.
- **Closures for numerical methods**: Solvers, integrators, and ODE stepper accept `Fn` closures, making them composable with the evaluator.
- **Inline tests**: Each module has `#[cfg(test)] mod tests` — no separate test files. 129 tests total.
- **Error handling**: `MathError` (via `thiserror`) in the library, `anyhow` in the binary.
- **Prelude**: `mathr::prelude` re-exports the most common types for ergonomic `use mathr::prelude::*;`.

## Deployment

The crate produces both a library (`mathr`) and a binary (`mathr`). The binary is a thin CLI wrapper around library functions. No network, no state, no config files — just math.
