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
│   ├── laurent.rs      # Laurent series expansion around a point
│   ├── special.rs      # Gamma, Beta, erf, erfc, sinc, incomplete gamma P,
│   │                   # Bessel functions J_0, J_1, J_n
│   ├── fastmath.rs     # Chebyshev-based fast approximations (sin, cos, tan, exp,
│   │                   # log, sqrt, pow) with argument reduction + ChebyshevApprox struct
│   ├── bigint.rs       # Arbitrary-precision integers (num-bigint): primality,
│   │                   # factorization (Pollard's rho), GCD, factorial, Fibonacci,
│   │                   # binomial, mod_pow, totient
│   ├── autodiff.rs     # Automatic differentiation via dual numbers: Dual type,
│   │                   # derivative, gradient, jacobian, Expr evaluation
│   ├── plot.rs         # PNG plotting via plotters (single, multi, scatter),
│   │                   # plot_*_to_bytes for inline notebook rendering
│   ├── repl.rs         # interactive REPL (rustyline) + REPL dispatch + dispatch_steps
│   ├── notebook.rs     # .mnb notebook format, JSON parse/serialize, cell eval,
│   │                   # cell types (Math/Text), cell reordering, shared context
│   ├── server.rs       # minimal HTTP server for web notebook UI, shared context,
│   │                   # /api/reset, /api/context endpoints
│   ├── webui.html      # single-page web UI (KaTeX rendering, step-by-step solving)
│   ├── rational.rs     # exact rational arithmetic + eval_rational for AST evaluation
│   └── error.rs        # MathError + Result alias
├── examples/
│   ├── *.rs            # Rust API examples (calculus, fft, matrix, numtheory, etc.)
│   └── notebooks/
│       └── *.mnb       # math notebook files (demo, calculus, fractions, solving, etc.)
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

### 4. Harvest to MEMORY.md
After each completed feature, extract patterns and best practices:
- **Success patterns**: What worked well and should be repeated
- **Anti-patterns**: What to avoid in future implementations
- **Math domain knowledge**: Numerical stability pitfalls, convergence criteria, edge cases (branch cuts, singularities, overflow/underflow)
- **Rust patterns**: mathr-specific conventions for the Expr AST, parser, evaluator, and CLI dispatch
- **Testing patterns**: How to assert on float results, fixture expressions, regression cases

Add these to `MEMORY.md` with clear categories and references to specific files/lines.

### 5. Loop Back to Step 1
Return to `TODO.md` and pick the next item. Repeat until the backlog is clear.

### 6. Audit and Optimize
After each batch of features, perform a quality pass:
- **Maintainability**: Are functions small and well-named? Is the module structure logical?
- **Leanness**: Remove dead code, unused imports, and speculative abstractions
- **Wiring**: Ensure all new features are properly integrated into `lib.rs`, the `prelude`, and `repl.rs` CLI dispatch
- **Small footprint**: Avoid unnecessary dependencies; prefer standard library or lightweight crates
- **Consistency**: Match existing code style and patterns (Rust 2021 edition, `thiserror` for errors, `num-traits` for numerics)
- **Numerical accuracy**: Verify algorithms against known references (e.g. special functions vs. tables, quadrature vs. analytic integrals)

### 7. Competitive Intelligence
Research similar Rust math libraries (e.g. `nalgebra`, `ndarray`, `rustfft`, `meval`, `symengine` bindings). Identify capabilities they have that this project lacks. Add the most valuable ones to the `TODO.md` brainstorming section. Prioritize features that provide clear competitive advantage.

### 8. Update Documentation
Keep all project docs aligned with the current implementation. Root docs (required):

- **`README.md`**: Quick start, CLI usage, feature list, crate API summary
- **`ARCHITECTURE.md`**: Module relationships, data flow, design decisions
- **`TODO.md`**: Mark completed items, move them to Done, keep brainstorming current
- **`SPEC.md`**: CLI subcommands, expression grammar, supported functions/constants
- **`MEMORY.md`**: Harvested patterns, domain knowledge, technical conventions (enhanced)

Update **`AGENTS.md`** if the loop itself evolves.

## Memory System (MEMORY.md)

### Purpose
`MEMORY.md` is the institutional knowledge repository that accelerates development by:
- **Preventing wheel reinvention**: Reuse proven patterns instead of guessing
- **Domain knowledge preservation**: Capture numerical and symbolic math rules that may be counter-intuitive
- **Onboarding acceleration**: New contributors (human or AI) can understand patterns quickly
- **Quality consistency**: Ensure all features follow established conventions

### Structure
Organize `MEMORY.md` into these sections:

#### 1. Expr AST & Parser Patterns
- AST node design (`num`, `var`, `neg`, `add`, `sub`, `mul`, `div`, `pow`, `func`)
- Canonicalization and equality rules
- Recursive-descent parser conventions (implicit multiplication, LaTeX/TeX input)
- Operator precedence and associativity handling

#### 2. Evaluator & Simplifier Patterns
- `Context` structure (vars, funcs, constants)
- Tree-walking evaluation conventions
- Constant folding and algebraic identity simplification rules
- Handling of undefined/singular expressions

#### 3. Symbolic & Calculus Patterns
- Differentiation rules and chain-rule application
- Indefinite integration heuristics and their limits
- Taylor and Laurent series expansion conventions
- Numerical quadrature convergence criteria (trapezoidal, Simpson, adaptive, Romberg)

#### 4. Solver & Matrix Patterns
- Root-finding convergence and bracketing (bisection, Newton–Raphson, secant, Durand–Kerner)
- Linear algebra decomposition conventions (LU, Cholesky, SVD, power iteration)
- Conditioning and numerical stability notes

#### 5. Number Theory & Special Functions
- Primality and factorization algorithm choices (Miller–Rabin, sieve)
- Modular arithmetic edge cases (inverse existence, CRT assembly)
- Special function implementations (Gamma, Beta, erf, Bessel) and reference tables for validation

#### 6. Testing Patterns
- Using `approx` for float comparisons and choosing tolerances
- Reference fixtures for symbolic/numeric verification
- CLI smoke test structure in `tests/integration.rs`
- Regression cases for parser and evaluator edge cases

#### 7. CLI & REPL Patterns
- `clap` single-string dispatch conventions
- REPL command routing and step-by-step solving (`dispatch_steps`)
- Notebook (`.mnb`) format and cell evaluation
- HTTP server and web UI integration

## Principles

- **Simplicity over flexibility**: Solve the problem at hand, not every hypothetical future problem
- **Surgical changes**: Touch only what you must; clean up only your own mess
- **Goal-driven**: Every change should have a verifiable success criterion
- **Test before ship**: No feature is complete until it has passing tests
- **Docs are code**: Documentation drift is a bug
- **Numerical fidelity**: Never compromise on mathematical accuracy for convenience
- **Memory first**: Always check `MEMORY.md` before starting a new feature
- **Pattern harvesting**: After success, update `MEMORY.md` to share the learning

## File Positioning and Value

### README.md
- **Value**: User-facing documentation and project overview
- **Audience**: Users, contributors, stakeholders
- **Position**: Entry point for anyone discovering the project
- **Focus**: Features, quick start, CLI usage, crate API summary, architecture summary

### TODO.md
- **Value**: Feature roadmap and backlog management
- **Audience**: Development team (human and AI agents)
- **Position**: Development planning and prioritization
- **Focus**: What to build next, what's done, competitive intelligence

### ARCHITECTURE.md
- **Value**: Module relationships, data flow, and design decisions
- **Audience**: Contributors maintaining or extending the crate
- **Position**: Structural reference for the codebase
- **Focus**: Module boundaries, data flow, deployment topology

### SPEC.md
- **Value**: Interface specification for the CLI and expression language
- **Audience**: Users and contributors integrating with mathr
- **Position**: Contract definition for inputs and outputs
- **Focus**: CLI subcommands, expression grammar, supported functions/constants

### MEMORY.md
- **Value**: Institutional knowledge and pattern library
- **Audience**: Development team (accelerates onboarding and consistency)
- **Position**: Development acceleration and quality consistency
- **Focus**: Proven patterns, domain knowledge, technical conventions

### AGENTS.md (this file)
- **Value**: Development process and workflow definition
- **Audience**: AI agents and human developers following the development loop
- **Position**: Process automation and continuous improvement
- **Focus**: How we work, the loop, memory system, principles
- **Update**: This file should be updated when the development loop itself evolves or when new process patterns emerge

### MEMORY.md
- **Value**: Institutional knowledge and pattern library
- **Audience**: Development team (accelerates onboarding and consistency)
- **Position**: Development acceleration and quality consistency
- **Focus**: Proven patterns, domain knowledge, technical conventions
- **Update**: Must be updated after each completed feature to capture patterns and lessons learned

## How These Files Work Together

1. **README.md** tells stakeholders what the project is and how to use it
2. **SPEC.md** defines the CLI and expression-language contract
3. **ARCHITECTURE.md** describes how the modules fit together
4. **TODO.md** tells developers what to build next (driven by competitive intelligence)
5. **AGENTS.md** tells agents how to work through the TODO items with quality and memory
6. **MEMORY.md** captures what we learned so we don't repeat mistakes

The loop reinforces these files:
- Complete TODO → Test → Harvest to MEMORY → Optimize → Research → Update TODO

This creates a flywheel of continuous improvement with institutional knowledge preservation.
