# Agent Development Loop

This document defines the continuous improvement cycle for the **mathr** crate — a pure-Rust mathematics library and CLI for symbolic and numerical computation.

## Project Structure

```
.
├── src/
│   ├── lib.rs          # crate root, module declarations, prelude re-exports
│   ├── main.rs         # CLI entry point (single string dispatch via clap)
│   ├── expr.rs         # Expr AST (num, var, neg, add, sub, mul, div, pow, func) + canonicalize/equals
│   ├── parser.rs       # recursive-descent parser → Expr (with implicit multiplication, LaTeX/TeX)
│   ├── eval.rs         # tree-walking evaluator with Context (vars, funcs, constants)
│   ├── simplify.rs     # constant folding + algebraic identity simplifier
│   ├── symbolic.rs     # symbolic differentiation + indefinite integration rules
│   ├── calculus.rs     # numerical derivatives, trapezoidal/Simpson/adaptive quadrature,
│   │                   # gradients, Romberg with Richardson extrapolation
│   ├── solver.rs       # bisection, Newton–Raphson, secant, Durand–Kerner polynomial roots,
│   │                   # Newton's method for nonlinear systems
│   ├── fft.rs          # Cooley–Tukey radix-2 FFT (forward, inverse, 2D, real-input, spectra),
│   │                   # convolution, cross-correlation, window functions
│   ├── complex.rs      # Complex<T> type with arithmetic, abs, arg, powers
│   ├── interpolate.rs  # Lagrange, Newton, linear, cubic spline, Chebyshev,
│   │                   # Legendre polynomials, Gauss–Legendre quadrature
│   ├── matrix.rs       # Matrix type: arithmetic, determinant, inverse, linear solve,
│   │                   # rank, trace, transpose, LU, Cholesky, SVD, power iteration
│   ├── stats.rs        # mean, median, variance, stddev, quartiles, correlation, regression
│   ├── numtheory.rs    # GCD, LCM, extended GCD, modular inverse, modular exponentiation,
│   │                   # primality, factorization, sieve, binomial, factorial, Fibonacci,
│   │                   # Euler's totient, Miller–Rabin, Jacobi symbol, continued fractions,
│   │                   # linear Diophantine solver, discrete logarithm, CRT
│   ├── ode.rs          # Euler, RK4, RK4 systems, adaptive RKF45
│   ├── taylor.rs       # symbolic Taylor series expansion
│   ├── special.rs      # Gamma, Beta, erf, erfc, sinc, incomplete gamma P,
│   │                   # Bessel functions J_0, J_1, J_n
│   ├── plot.rs         # PNG plotting via plotters (single, multi, scatter)
│   ├── repl.rs         # interactive REPL (rustyline) + REPL dispatch
│   └── error.rs        # MathError + Result alias
├── examples/
│   └── debug_solve.rs  # example: using the solver API programmatically
├── tests/
│   └── integration.rs  # end-to-end CLI smoke tests (parse, dispatch, output)
├── Cargo.toml          # package metadata, deps (clap, anyhow, thiserror, rustyline,
│                       # plotters, num-traits), dev-dep (approx)
└── Cargo.lock
```

## The Loop

### 1. Complete Remaining TODO Items
Pick the next highest-priority item from `TODO.md` (or `ARCHITECTURE.md` if the task is architectural). Implement it with minimal, focused changes. Do not add speculative features.

### 2. Create Tests and Examples
For every new capability:
- Add inline `#[cfg(test)] mod tests` in the relevant source file — exercise the feature end-to-end
- Add unit tests for core math logic (use `approx` for float comparisons)
- Provide a minimal usage example in `examples/` if the feature is library-facing
- Add a CLI smoke test to `tests/integration.rs` if there's a CLI dispatch path

### 3. Ensure `cargo test` Passes
Run the full test suite:
```bash
cargo test                  # all inline unit tests + integration tests
cargo test --examples       # examples compile and run
cargo clippy                # lint pass (warnings acceptable but noted)
```
Fix any failures before proceeding.

### 4. Loop Back to Step 1
Return to `TODO.md` and pick the next item. Repeat until the backlog is clear.

### 5. Audit and Optimize
After each batch of features, perform a quality pass:
- **Maintainability**: Are functions small and well-named? Is the module structure logical?
- **Leanness**: Remove dead code, unused imports, and speculative abstractions
- **Wiring**: Ensure all new features are properly integrated into `lib.rs`, the `prelude`, and `repl.rs` CLI dispatch
- **Small footprint**: Avoid unnecessary dependencies; prefer standard library or lightweight crates
- **Consistency**: Match existing code style and patterns (Rust 2021 edition, `thiserror` for errors, `num-traits` for numerics)

### 6. Competitive Intelligence
Research similar Rust math libraries (e.g. `nalgebra`, `ndarray`, `rustfft`, `meval`, `symengine` bindings). Identify capabilities they have that this project lacks. Add the most valuable ones to the `TODO.md` brainstorming section. Prioritize features that provide clear competitive advantage.

### 7. Update Documentation
Keep all project docs aligned with the current implementation. Root docs (required):

- **`README.md`**: Quick start, CLI usage, feature list, crate API summary
- **`ARCHITECTURE.md`**: Module relationships, data flow, design decisions
- **`TODO.md`**: Mark completed items, move them to Done, keep brainstorming current
- **`SPEC.md`**: CLI subcommands, expression grammar, supported functions/constants

Update **`AGENTS.md`** if the loop itself evolves.

## Principles

- **Simplicity over flexibility**: Solve the problem at hand, not every hypothetical future problem
- **Surgical changes**: Touch only what you must; clean up only your own mess
- **Goal-driven**: Every change should have a verifiable success criterion
- **Test before ship**: No feature is complete until it has passing tests
- **Docs are code**: Documentation drift is a bug