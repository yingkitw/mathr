# Memory

This file captures institutional knowledge, patterns, and best practices for the mathr development team. Always check this file before implementing new features.

## Expr AST & Parser Patterns

### AST Node Design
- **Enum-based AST**: `Expr` enum with variants: `num`, `var`, `neg`, `add`, `sub`, `mul`, `div`, `pow`, `func`
- **Canonicalization**: `expr.rs:canonicalize()` normalizes expressions by sorting commutative operands (add, mul) and simplifying identities (x+0=x, x*1=x)
- **Equality checking**: `expr.rs:equals()` uses canonicalized forms for structural equality, not just syntax
- **Function calls**: `func` variant stores name and `Vec<Expr>` arguments, supporting variable-length args

### Parser Conventions
- **Recursive-descent**: `parser.rs:Parser` implements precedence climbing with lookahead tokens
- **Implicit multiplication**: `2x`, `3(x+1)`, `(x)(y)` → `2*x`, `3*(x+1)`, `x*y`
- **LaTeX/TeX support**: Parser recognizes `\frac`, `\sqrt`, `\sin`, `\pi`, `\left(\right)`, `^{...}`, `\Gamma`, `\log_2`, etc.
- **Operator precedence**: Standard mathematical precedence with right-associative exponentiation
- **Number parsing**: Supports integers, decimals, and scientific notation (`1.5e3`, `2E-2`)

## Evaluator & Simplifier Patterns

### Context Structure
- **Unified evaluation**: `eval.rs:Context` holds `HashMap<String, f64>` for variables, `HashMap<String, fn(&[f64]) -> f64>` for functions
- **Built-in constants**: pi, e, tau, inf, nan registered in default context
- **User bindings**: REPL `let x = expr` and `fn f(x) = expr` persist in mutable context

### Tree-Walking Evaluation
- **Recursive traversal**: `eval.rs:eval()` walks AST depth-first, substituting variables and evaluating functions
- **Error propagation**: Returns `Result<f64, MathError>` with domain validation (sqrt negative, log zero)
- **Exact rational mode**: `rational.rs:eval_rational()` walks same AST but returns `Rational` when all leaves are integers

### Constant Folding & Simplification
- **Identity application**: `simplify.rs:simplify()` applies x+0=x, x*1=x, x-x=0, x/x=1 (when x ≠ 0)
- **Constant propagation**: Reduces sub-expressions with only numeric leaves to single constants
- **Algebraic identities**: Applies sin^2+cos^2=1, log(e^x)=x, e^(ln x)=x where safe
- **Fallback on unbound variables**: Notebook evaluation falls back to simplification when eval fails due to undefined variables

## Symbolic & Calculus Patterns

### Differentiation Rules
- **Product/quotient/chain**: `symbolic.rs:differentiate()` implements all standard rules
- **Partial derivatives**: `symbolic.rs:gradient()` collects all free variables and computes partials for each
- **Variable binding**: Differentiation is with respect to a specified variable; all others treated as constants
- **Simplification chain**: After differentiation, always run `simplify()` to reduce complexity

### Integration Heuristics
- **Pattern matching**: `symbolic.rs:integrate()` recognizes polynomial, exponential, trigonometric, and inverse-trigonometric primitives
- **Indefinite only**: Does not handle definite integrals or boundary conditions
- **Limited scope**: Cannot integrate arbitrary expressions; returns original if no pattern matches

### Series Expansions
- **Taylor series**: `taylor.rs:taylor_series()` uses symbolic differentiation + evaluation at expansion point
- **Laurent series**: `laurent.rs:laurent_series()` expands `g(x) = (x-a)^k·f(x)` via Taylor, then divides by `(x-a)^k`
- **Truncation**: Both series truncate at specified order; no radius of convergence analysis

### Numerical Quadrature
- **Trapezoidal**: Simple composite rule, O(h²) accuracy
- **Simpson's**: Higher-order O(h⁴) accuracy, used for Fourier coefficient computation
- **Adaptive**: `calculus.rs:integrate_adaptive()` recursively subdivides until tolerance met
- **Romberg**: Richardson extrapolation on trapezoidal rule, exponential convergence for smooth functions
- **Monte Carlo**: `calculus.rs:monte_carlo_integrate_1d/nd()` uses reproducible LCG-based sampling with standard error estimates
- **Fourier**: `calculus.rs:fourier_series()` computes coefficients via Simpson integration

## Solver & Matrix Patterns

### Root-Finding Convergence
- **Bracketing required**: Bisection requires f(a)·f(b) < 0
- **Derivative-based**: Newton–Raphson needs initial guess and analytical derivative
- **Secant method**: Derivative-free alternative using finite differences
- **Polynomial roots**: `solver.rs:polynomial_roots()` (Durand–Kerner) finds all complex roots simultaneously
- **VAS isolation**: `solver.rs:isolate_real_roots()` uses i128 exact arithmetic for real root isolation with integer coefficients
- **Nonlinear systems**: `solver.rs:newton_system()` uses central-difference Jacobian, requires good initial guess

### Linear Algebra Decompositions
- **Matrix storage**: Row-major `Vec<f64>` with `rows` and `cols` fields
- **LU decomposition**: Partial pivoting for stability, `matrix.rs:lu()` returns L and U matrices
- **Cholesky**: Requires symmetric positive-definite; `matrix.rs:cholesky()` returns lower triangular L where A = L·Lᵀ
- **SVD**: Jacobi rotations diagonalize; `matrix.rs:svd()` returns U, Σ, Vᵀ
- **Power iteration**: Finds dominant eigenpair; `matrix.rs:power_iteration()` for largest eigenvalue
- **QR algorithm**: `matrix.rs:sym_eig()` for full symmetric eigenvalue decomposition via Householder tridiagonalisation + Wilkinson-shift QR
- **Hessenberg**: `matrix.rs:hessenberg()` for general matrices A = Q·H·Qᵀ via Householder reflections
- **Real Schur**: `matrix.rs:schur()` for real Schur decomposition A = Q·T·Qᵀ via shifted QR on Hessenberg form
- **Tikhonov regularisation**: `matrix.rs:tikhonov_solve()` solves ill-conditioned systems via (AᵀA + λI)x = Aᵀb

## Number Theory & Special Functions

### Primality & Factorization
- **Miller–Rabin**: `numtheory.rs:is_prime()` uses deterministic bases for n < 3.3e24
- **Sieve of Eratosthenes**: `numtheory.rs:sieve()` generates primes up to n
- **Factorization**: Trial division up to sqrt(n) for small numbers
- **GCD/LCM**: Euclidean algorithm for GCD, derived LCM
- **Extended GCD**: Returns coefficients for Bezout identity

### Modular Arithmetic
- **Inverse existence**: `numtheory.rs:mod_inverse()` checks gcd(a, n) = 1 first
- **CRT assembly**: `numtheory.rs::crt()` solves simultaneous congruences when moduli are coprime
- **Discrete logarithm**: Baby-step giant-step algorithm for finding x in g^x ≡ h (mod p)
- **Jacobi symbol**: `numtheory.rs::jacobi()` for quadratic residuosity testing
- **Linear Diophantine**: `numtheory.rs::diophantine()` solves a·x + b·y = c
- **Continued fractions**: `numtheory.rs::continued_fraction()` for rational and real-valued approximants

### Special Function Implementations
- **Gamma/Lanczos**: `special.rs:gamma()` uses Lanczos approximation for x > 0
- **Incomplete gamma**: `special.rs::incomplete_gamma_p()` for regularized lower incomplete gamma
- **erf/erfc**: Via incomplete gamma relationship
- **Bessel functions**: `special.rs::bessel_j0/j1/jn()` use Maclaurin series + asymptotic + forward recurrence
- **sinc**: Normalized sinc function sin(πx)/(πx)

## Testing Patterns

### Float Comparisons
- **Use approx crate**: `dev-dependencies: approx = "0.5"` for `assert_abs_diff_eq!`, `assert_relative_eq!`
- **Tolerance selection**: Use absolute tolerance for near-zero values, relative for non-zero
- **Edge cases**: Test at boundaries (0, ±inf, NaN, very large/small numbers)

### Reference Fixtures
- **Symbolic verification**: Compare derivative of known function with expected result
- **Numeric benchmarks**: Test quadrature against analytically integrable functions
- **Regression cases**: Capture parser/evaluator edge cases as unit tests

### Test Structure
- **Inline tests**: Each module has `#[cfg(test)] mod tests` section
- **Integration tests**: `tests/integration.rs` for CLI smoke tests (parse, dispatch, output)
- **Example tests**: `cargo test --examples` verifies examples compile and run

## CLI & REPL Patterns

### Dispatch Architecture
- **Single-string dispatch**: `clap` with `Arg<Action::SetTrue>` for flags
- **Command routing**: `repl.rs:dispatch()` routes to library functions based on first token
- **Step-by-step**: `repl.rs:dispatch_steps()` returns `Vec<String>` of intermediate steps for educational display

### Notebook Integration
- **Context persistence**: `notebook.rs:eval_all()` uses `&mut Context` across cells
- **Cell types**: `CellType::Math` (evaluated) and `CellType::Text` (Markdown-rendered)
- **JSON serialization**: `.mnb` format with `cells` array, each cell has `id`, `input`, `output`, `cell_type`
- **Server API**: `server.rs:NotebookServer` serves `/api/eval`, `/api/notebook`, `/api/reset`, `/api/context`

### REPL Commands
- **Variable binding**: `let x = expr` stores in context
- **Function definition**: `fn f(x) = expr` stores closure in context
- **Context inspection**: `vars` and `funcs` commands show bindings
- **Reset**: `clear` command resets context to defaults

## Fast Approximation Patterns

### Chebyshev Approximations
- **Argument reduction**: `fastmath.rs:fast_sin/cos/tan` reduce input to principal range before approximation
- **ChebyshevApprox struct**: Stores coefficients and approximation degree
- **Clenshaw evaluation**: Efficient evaluation of Chebyshev series
- **Trade-off**: Accuracy ~1e-6 to 1e-10 vs. standard library speed

## Code Conventions

### Error Handling
- **Library errors**: `error.rs:MathError` enum with `thiserror` derives
- **Binary errors**: `anyhow` for CLI error context
- **Result propagation**: Use `?` operator consistently

### Type Design
- **Generic where beneficial**: `Complex<T>`, `Rational` (wraps i64/i64 with i128 intermediates)
- **Closure-based APIs**: Numerical methods accept `Fn(f64) -> f64` for composability
- **Prelude re-exports**: `lib.rs:pub use` common types for ergonomic imports

### Documentation
- **Doc comments**: `///` for public APIs with examples
- **Inline comments**: Minimal; prefer self-documenting code
- **README features**: Keep feature list aligned with TODO.md Done section

## Numerical Stability Notes

### Conditioning
- **Ill-conditioned systems**: Use Tikhonov regularisation for Hilbert matrices
- **Pivot strategies**: LU uses partial pivoting; Cholesky requires positive-definite check

### Convergence Criteria
- **Default tolerances**: 1e-10 for most iterative methods
- **Max iterations**: 100 for solvers, 1e6 for adaptive quadrature
- **Divergence detection**: Return `Err(MathError::NotConvergent)` when limits exceeded

### Edge Cases
- **Singularities**: Guard against division by zero, sqrt of negative, log of zero
- **Overflow/underflow**: Check for f64 limits in number theory (factorials, binomials)
- **NaN propagation**: Allow NaN in intermediate results but validate at boundaries

## Performance Considerations

### Algorithm Selection
- **FFT**: Cooley–Tukey radix-2 for powers of 2; `fft.rs::rfft()` exploits real symmetry
- **Matrix operations**: Use decomposition for repeated solves (LU, Cholesky)
- **Root finding**: Bisection for robustness, Newton–Raphson for speed when derivative available

### Memory Management
- **Allocation minimization**: Reuse buffers where possible in FFT
- **Clone avoidance**: Pass references (`&[f64]`) instead of owned vectors
- **Stack allocation**: Prefer fixed-size arrays for small matrices

## Competitive Intelligence

### Similar Rust Crates
- **nalgebra**: More comprehensive linear algebra, heavier dependency footprint
- **ndarray**: n-dimensional arrays, different API design
- **rustfft**: Faster FFT with SIMD, but mathr has simpler pure-Rust implementation
- **meval**: Expression evaluation only, no symbolic or numerical capabilities
- **symengine**: Symbolic computation via C++ bindings, not pure Rust

### mathr Advantages
- **Pure Rust**: No C/C++ dependencies
- **Symbolic + numerical**: Rare combination in single crate
- **Small footprint**: Minimal dependencies (clap, anyhow, thiserror, rustyline, plotters, num-traits)
- **Education-friendly**: Step-by-step solving, web notebook, REPL

### Feature Gaps (Brainstorming)
- **Arbitrary precision**: BigDecimal for exact arithmetic beyond i64
- **Interval arithmetic**: Rigorous bounds for numerical methods
- **Automatic differentiation**: Dual numbers for gradient computation
- **GPU acceleration**: wgpu for FFT, matrix operations
- **More special functions**: Elliptic integrals, hypergeometric, polygamma
- **Symbolic limits**: L'Hôpital's rule, asymptotic analysis
- **ODE solvers**: Implicit methods (BDF), stiff solvers
- **Sparse matrices**: Compressed storage formats
- **Interpolation**: More spline types (B-spline, Hermite)
- **Plotting**: Surface plots, 3D visualization, animation

## Integration Patterns

### Adding New Features
1. **Implement in module**: Add function to appropriate `.rs` file
2. **Add tests**: Inline `#[cfg(test)]` + integration test if CLI exposed
3. **Update SPEC.md**: Document syntax and behavior
4. **Wire in repl.rs**: Add dispatch case and optional step-by-step
5. **Update TODO.md**: Move to Done section
6. **Add to prelude**: Export from `lib.rs` if library-facing
7. **Harvest patterns**: Update this MEMORY.md with lessons learned

### Module Dependencies
- **expr.rs** is foundational: parser → expr → {eval, symbolic, simplify}
- **calculus.rs** and **solver.rs** use eval for closure generation
- **taylor.rs** and **laurent.rs** depend on symbolic + eval
- **fft.rs** depends on complex.rs
- **plot.rs** depends on eval + expr
- **repl.rs** uses all modules
- **server.rs** uses repl dispatch for notebook evaluation
- **error.rs** used by all modules

### Testing Workflow
```bash
# Run all tests
cargo test

# Test specific module
cargo test parser::tests

# Run integration tests
cargo test --test integration

# Run examples
cargo test --examples

# Lint pass
cargo clippy
```

## Common Pitfalls

### Parser
- **Implicit multiplication ambiguity**: `2 3` parsed as `2*3`, but avoid in practice
- **Function argument parsing**: Distinguish `f x y` (parse error) from `f(x, y)` (function call)

### Evaluation
- **Undefined variables**: Return `Err(MathError::UnknownVariable)` rather than panicking
- **Domain errors**: Validate before computation (sqrt negative, log zero)

### Symbolic
- **Variable scope**: Ensure differentiation variable exists in expression
- **Integration limits**: Return original expression when no rule matches

### Numerical
- **Convergence failures**: Always provide error path with max iteration check
- **Stability issues**: Use decomposition for ill-conditioned linear systems

### Matrix
- **Dimension mismatches**: Validate row/col counts before operations
- **Singular matrices**: Check for near-zero determinants before inversion

## Development Workflow

1. **Check MEMORY.md**: Review patterns before implementing
2. **Pick TODO item**: Select next highest-priority feature
3. **Implement minimally**: Focused changes without speculative features
4. **Add tests**: Inline unit tests + integration tests if CLI-exposed
5. **Run cargo test**: Ensure all tests pass
6. **Harvest to MEMORY.md**: Extract patterns and domain knowledge
7. **Update docs**: Align README, SPEC, TODO, ARCHITECTURE
8. **Loop**: Return to TODO for next item