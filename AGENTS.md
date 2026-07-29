# Agent Development Loop

This document defines the continuous improvement cycle for the **mathr** crate — a pure-Rust mathematics library and CLI for symbolic and numerical computation (FFT, calculus, equation solving, plotting).

## Project Structure

```
.
├── src/
│   ├── lib.rs          # crate root, module declarations, prelude re-exports
│   ├── main.rs         # CLI entry point (clap subcommands: eval, diff, simplify, integrate, solve, poly-roots, plot, fft, taylor, stats, gcd, lcm, is-prime, factor, fib, binom, fact, mr-prime, gamma, erf, conv, repl)
│   ├── expr.rs         # Expr AST (num, var, binary, unary, call, pow)
│   ├── parser.rs       # recursive-descent parser → Expr
│   ├── eval.rs         # tree-walking evaluator with Context (vars, funcs, constants)
│   ├── simplify.rs     # constant folding + algebraic identity simplifier
│   ├── symbolic.rs     # symbolic differentiation rules
│   ├── calculus.rs     # numerical derivatives, Simpson/trapezoidal/adaptive quadrature, gradients
│   ├── solver.rs       # bisection, Newton–Raphson, secant, Durand–Kerner polynomial roots
│   ├── fft.rs          # Cooley–Tukey radix-2 FFT (forward, inverse, 2D, real-input, spectra)
│   ├── complex.rs      # Complex<T> type with arithmetic, abs, arg, powers
│   ├── interpolate.rs  # Lagrange, Newton divided-difference, linear interpolation
│   ├── matrix.rs       # Matrix type, arithmetic, determinant, inverse, linear solve, trace
│   ├── stats.rs        # mean, median, variance, stddev, quartiles, correlation, regression
│   ├── numtheory.rs    # GCD, LCM, primality, factorization, binomial, factorial, Fibonacci, sieve, totient, Miller–Rabin, CRT
│   ├── ode.rs          # Euler, RK4, RK4 systems, adaptive RKF45
│   ├── taylor.rs       # symbolic Taylor series expansion
│   ├── special.rs      # Gamma, Beta, erf, erfc, sinc, incomplete gamma P
│   ├── plot.rs         # PNG plotting via plotters (single, multi, scatter)
│   ├── repl.rs         # interactive REPL (rustyline)
│   └── error.rs        # MathError + Result alias
├── examples/
│   └── debug_solve.rs  # example: using the solver API programmatically
├── Cargo.toml          # package metadata, deps (clap, anyhow, thiserror, rustyline, plotters, num-traits), dev-dep (approx)
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

### 3. Ensure `cargo test` Passes
Run the full test suite:
```bash
cargo test                  # all inline unit tests
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
- **Wiring**: Ensure all new features are properly integrated into `lib.rs`, `main.rs` CLI subcommands, and the `prelude`
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
