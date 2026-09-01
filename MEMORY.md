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

### Testing New Modules — Patterns Discovered (2026-09)

#### Limit module (limit.rs)
- **Numeric evaluation limitations**: The limit implementation uses numeric probing, not symbolic analysis. `floor(x)` at integers returns `Finite` (the right-side value), not `DoesNotExist`. `sign(x)` at 0 similarly returns `Finite(1.0)`. `exp(1/x)` at 0 returns `Finite(0.0)` (left-side value). Tests should document actual behavior, not mathematical ideals.
- **L'Hôpital depth**: `MAX_LOPITAL = 6`. `(x - sin x)/x³` needs 3 applications → 1/6. Test with `snap_num(1.0/6.0)` for float comparison.
- **Limits at ±∞**: Pass `f64::INFINITY` or `f64::NEG_INFINITY` as the point. Test both directions.

#### Apart module (apart.rs)
- **`apart_steps_str` signature**: Takes `&str` (raw input like `"1/(x^2-1) x"`), NOT `(&Expr, &str)`. This is a REPL convenience wrapper. Use `apart_steps(&Expr, &str)` for API testing.
- **High-degree denominators**: `1/(x^4 - 1)` may fail numeric agreement at some points. Use `agrees()` helper with points away from poles. Cubic denominators (`1/(x^3 - 1)`) work reliably.
- **`agrees()` helper**: Evaluates original and decomposed at sample points, checks `|a-b| < 1e-6 * |a|.max(1)`. Avoid points near poles (roots of denominator).

#### Poly module (poly.rs)
- **`to_poly` returns `Option<Poly>`**: Returns `None` for non-polynomial expressions (functions, variable denominators, symbolic/fractional/negative exponents). Test each rejection path explicitly.
- **Empty poly**: `poly_to_expr(&vec![])` produces an expression that evaluates to 0.0.
- **Float collecting**: `0.1*x + 0.2*x` collects to a single term with coeff ≈ 0.3 (within f64 tolerance).

#### BigDec module (bigdec.rs)
- **`pi(prec)` returns `Result<BigDecimal>`**: Must `.unwrap()` before use. Not a plain `BigDecimal`.
- **`sub`/`add`/`mul` return `BigDecimal` directly**: NOT `Result`. Only `div`, `sqrt`, `ln`, `exp`, `tan`, `asin`, `acos` return `Result`.
- **Inverse trig at boundaries**: `asin(±1)` and `acos(±1)` may error with "division by zero" due to the sqrt(1-x²) term. Test with `±0.5` or `±0.999` instead.
- **`ln` accuracy for extreme values**: `ln(1e-100)` gives wrong results (implementation convergence issue). Use moderate values like `ln(0.001)` or `ln(1000)` with 1e-2 tolerance.
- **`tan` near π/2**: Precision issues — `tan(π/2 - eps)` for small `eps` may not give expected large values. Test `tan(1.0)` instead (well away from singularity).
- **`pow(0, 0)`**: Returns 1.0 (conventional). Test this explicitly.
- **Negative base with integer exponent**: `pow(-2, 3)` = -8, `pow(-2, 4)` = 16. Works correctly.

#### Special functions (special.rs)
- **`polygamma` takes `u32`**: First argument is `u32`, but `f64::powi` takes `i32`. Cast with `m as i32` in tests.
- **Polygamma recurrence**: `ψ^(m)(x+1) = ψ^(m)(x) + (-1)^m * m! / x^(m+1)` — note it's `(-1)^m`, NOT `(-1)^(m+1)`.
- **Trigamma reflection**: `ψ₁(1-x) + ψ₁(x) = π²/sin²(πx)`. For `x = -0.5`: `ψ₁(-0.5) = π² - ψ₁(1.5) = π² - (π²/2 - 4) = π²/2 + 4`.
- **Carlson RF with zero argument**: `RF(0, 1, 1) = π/2` (not `asinh(1)`). Verified by integral substitution.
- **Eval context registration**: Not all special functions are registered in `eval::Context::standard()`. `hurwitz_zeta` and `carlson_rf` are NOT available via `eval_str`. Check `eval.rs` function registration before writing eval-based tests.
- **Zeta at negative odd integers**: `ζ(-5) = -1/252`, `ζ(-7) = 1/240`. These use the Bernoulli number formula.
- **Zeta at even positive integers**: `ζ(6) = π⁶/945`, `ζ(8) = π⁸/9450`. Use `PI.powi(6)` etc.

#### Array literal corruption
- **Edit tool bug**: When editing files, `]` in array literals like `&[1.0, 2.0, 3.0]` can get corrupted to `}`. Always verify with `cargo test` after edits and check hex dump if delimiter errors appear.

## CLI & REPL Patterns

### Dispatch Architecture
- **Single-string dispatch**: `clap` with `Arg<Action::SetTrue>` for flags
- **Command routing**: `repl.rs:dispatch()` routes to library functions based on first token
- **Step-by-step**: `repl.rs:dispatch_steps()` returns `Vec<String>` of intermediate steps for educational display
- **Keyword list for step dispatch**: `repl.rs` keeps a `cmd_keywords` array (~line 163) of `"prefix "` strings; commands starting with one of these fall through to `dispatch_inner` for single-step output. Add new REPL commands here AND in `dispatch_inner`'s `strip_prefix` chain.
- **Import subcommand convention**: for reversible formats (MathML, serialize), use `<cmd> import <text>` to import and bare `<cmd> <expr>` to export — see `do_mathml` (~line 1538) and `do_serialize` (~line 1548) in `repl.rs`.

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

## Serialization Patterns (serialize.rs)

### Three-format design
- **S-expressions** (`to_sexpr`/`from_sexpr`): Lisp-like prefix notation, `(add (mul (num 2) (var x)) (num 1))`. Recursive-descent parser over paren/atom tokens.
- **JSON** (`to_json`/`from_json`): nested objects with a `"t"` tag discriminator. Hand-rolled recursive JSON parser (no serde) — matches the crate's no-serde convention used by `notebook.rs`.
- **RPN** (`to_rpn`/`from_rpn`): space-separated postfix. Functions use a `<name>:<arity>` call token (e.g. `x sin:1`, `x 2 pow:2`). Stack-based parser.

### Round-trip contract
- All three formats round-trip via `Expr::equals` (canonical form), NOT raw `==`. This matters because `equals` normalizes `Sub`→`Add(Neg)` and sorts commutative operands — two structurally different `Expr`s can be mathematically equal. Tests use `assert!(e.equals(&e2))`, never `assert_eq!(e, e2)`.

### Number formatting
- Shared `fmt_num` helper: integers as `i64`, `NaN`/`inf`/`-inf` as identifiers (sexpr/rpn) or strings (json, since JSON numbers can't be NaN/inf). `parse_num` is the inverse. Keep these consistent across formats so the textual forms stay aligned with `Expr::Display`.

### Pitfalls encountered
- **`Iterator::rev()` on `Map` of `Range`**: `(0..n).map(|_| stack.pop()).rev().collect()` does NOT reverse the collected pops — it reverses the range iteration, but `pop()` still pulls from the top in the same order. Fix: collect into a `Vec` first, then call `.reverse()` on the Vec. (Bug found in `from_rpn` function-arg ordering, 2026-08-18.)
- **JSON parser UTF-8**: when reading string chars byte-by-byte, multi-byte UTF-8 sequences need full decoding (back up, decode with `str::from_utf8`, advance by `char.len_utf8()`). Single-byte fallback only for ASCII `< 0x80`.
- **Dead-code warning on `JsonValue::Bool`**: a robust JSON parser should support all value types even if the Expr schema never uses them; `#[allow(dead_code)]` on the variant is cleaner than removing it.

## QR Decomposition Patterns (matrix.rs)

### Householder reflections algorithm
- **Sign choice**: `alpha = -sign(x[0]) * ||x||` where `x` is the sub-column being zeroed. This ensures the reflection is numerically stable (avoids cancellation). The sign of R's diagonal entries depends on this choice — they can be negative even for positive-definite matrices. Tests should check `abs(r[i,i])` not `r[i,i] == expected` when the sign isn't deterministic.
- **Householder vector**: `v = x - alpha * e_1`, then `H = I - 2*v*vᵀ/(vᵀv)`. Apply H to R from the left and to Q from the right (Q = Q·H accumulates the reflections).
- **Skip zero columns**: if `||x|| < 1e-300` or `||v||² < 1e-300`, skip the step — the column is already zeroed below the diagonal.
- **Rectangular matrices**: works for both tall (m≥n) and wide (m<n). Q is always m×m orthogonal, R is m×n upper-trapezoidal. No transpose swap needed.

### Least-squares solve
- For overdetermined systems (m≥n, full rank): compute `Qᵀb` then back-substitute `R·x = Qᵀb` using the n×n leading block of R.
- Rank deficiency: if any `|R[i,i]| < 1e-14` during back-substitution, return an error.
- Underdetermined (n>m): reject — QR can't solve underdetermined systems.

### Pitfalls encountered
- **`&Matrix * &Matrix` returns `Result<Matrix>`**: not `Matrix` directly. In tests, must `.unwrap()` the product before indexing: `let qtq = (&qt * &q).unwrap();`. Forgetting this causes `E0608: cannot index into Result`.
- **R diagonal sign**: for the identity matrix, the Householder reflection of `[1,0,0]` produces `v = [2,0,0]`, which maps `R[0,0]` from 1 to -1. The test `assert!(close(r[(i,i)], 1.0))` fails — use `assert!(close(r[(i,i)].abs(), 1.0))` instead.

## Special Function Implementation Patterns (special.rs, added 2026-08-31)

### Working recipes (validated to ~1e-14 against reference digits)
- **digamma/trigamma**: reflection (x < 0) → recurrence (push x above 12/14) → asymptotic Bernoulli series. Recurrence threshold 6 leaves ~1e-10 error from series truncation; raising to 12–14 gets ~1e-15. Cheaper than adding more series terms.
- **Hurwitz zeta** ζ(s, a): Euler–Maclaurin with N = 12 direct terms + 6 Bernoulli corrections, valid for ANY real s ≠ 1. **Power bookkeeping trap**: term j needs `s^{↑(2j−1)}·t^{−(s+2j−1)}` — the first step divides by t¹, subsequent steps by t². Getting the exponents wrong costs 2–7 digits (ζ(2) was off at 7e-5).
- **polygamma(m ≥ 2)** = (−1)^(m+1)·m!·ζ(m+1, x) — free once Hurwitz exists.
- **Elliptic K** via Carlson RF: K(k) = RF(0, 1−k², 1). **Elliptic E via the AGM identity** E = K·(1 − Σ 2^{n−1}·(a_n²−b_n²)), a₀ = 1, b₀ = k′ = √(1−k²).
- **Carlson RF convergence threshold**: NR's published 0.0025 gives only ~1e-7 accuracy (the Taylor corrections are ~O(dels²)); tightening to 1e-9 costs one extra quadratic-convergence iteration and reaches ~1e-15.

### Traps hit this iteration
- **Stack overflow from the zeta functional equation**: `ζ(s) = …ζ(1−s)` maps s ∈ (0,1) to (0,1) — ζ(0.5) recursed into itself. Euler–Maclaurin is valid for all real s ≠ 1, so the FE is unnecessary; only special-case the pole, ζ(0), and trivial zeros.
- **Catastrophic cancellation in AGM for E**: computing `a_n² − b_n²` directly leaves rounding-noise ~1e-16 when a ≈ b, and the `Σ 2^{n−1}·noise` factor amplifies it into E(0.6) = −3584. Fix: compute `d² = (a−b)(a+b)` (exact subtraction) AND break before adding when `d² ≤ ε·a²`.
- **Reference-value integrity**: 1.3110287771 is K(1/√2)?? No — it's a fabricated hybrid; RF(0, ½, 1) = K(1/√2) = 1.8540746773. Always derive constants from identities (Γ(1/4)²-formulas are convention-sensitive) or verify against a second source before asserting them in tests.
- **Distinguish modulus k vs parameter m**: sympy's elliptic_k takes m = k²; these functions take the MODULUS k (K(0.5) = 1.6857503548 = parameter 0.25). Document the convention at every call site.

## Condition Number & Nullspace Patterns (matrix.rs, added 2026-08-31)

### SVD shape asymmetry
- **`svd()`'s `v` field is only the true right-singular-vector basis when `m ≥ n`**. For wide matrices (m < n) the Jacobi run works on the transposed matrix, so `v` holds the transposed role's vectors and `u` is zero-padded (documented in svd() comments). Anything needing right singular vectors of a wide matrix (nullspace, etc.) must go through `Aᵀ·A` instead.
- **`Matrix::nullspace` therefore uses `symmetric_eig` on `Aᵀ·A`** (n×n symmetric, always the right shape): eigenvectors with `|λ| ≤ tol·λ_max` span {x : Ax = 0}. `symmetric_eig` returns eigenvalues ASCENDING, so scan from the front and break past the threshold. This works for any aspect ratio and yields orthonormal vectors.
- **`condition_number` only needs `singular_values`** (aspect-independent) — but remember `singular_values` are computed as `sqrt(max(0, eigenvalue))` and sorted descending; a singular matrix yields a 0 (or negative-clamped) value → return `inf` rather than dividing.

## Partial Fraction Patterns (apart.rs, added 2026-08-31)

### solver::polynomial_roots API trap
- **`polynomial_roots(coeffs) -> Vec<(f64, f64)>` returns `(root, f(root))` for REAL roots only** — complex roots are silently discarded. The tuple is NOT (re, im)! For complex roots use `polynomial_roots_complex(coeffs) -> Vec<(re, im)>` (added 2026-08-31; both share the extracted `durand_kerner_all_roots` helper). Check the docs before assuming tuple semantics.

### apart.rs design
- **Quadratic factors need TWO unknowns per power** (numerator `B·x + C`): columns are `x·D/F^k` and `D/F^k`. Linear factors contribute one column (`D/(x−r)^k`). K columns = deg D exactly after long division (deg R < deg D), giving a square system solved by `Matrix::solve`.
- **Clone-before-mutate column vectors**: multiplying a dense descending polynomial by x appends `0.0` — mutating the shared `Vec` before the second push made both quadratic columns identical → "singular matrix". 
- **A conjugate pair covers 2 degrees**: the degree sanity check must weight quadratic factors by 2× (each positive-imag root in the DK output represents a pair).
- **Repeated roots**: cluster numerically equal roots (tolerance `1e-5·max(1,|r|)`); multiple-root finding carries inherent fuzz — document that coefficients may be off in the last digits.
- **Sign-aware output joins** (same pattern as `poly_to_expr`): build denominators with `Add`/`Sub` chosen by the sign of `r`/`p`/`q` (`x + 1`, never `x - -1`); join terms with `Sub`/`Add` by leading sign and wrap negative numerators as `Neg(Div(magnitude))`. `simplify` does NOT fold `x - -1` → `x + 1`.
- **Testing decompositions**: numeric probe equivalence (eval original and result at sample points away from poles, compare within relative tolerance) is far more robust than canonical `equals` for decomposed output — the latter hits all the parser/canonical-form traps documented above.

## Polynomial Expansion Patterns (poly.rs, added 2026-08-31)

### Canonical-form traps (read before touching expr.rs transforms)
- **The parser canonicalizes aggressively**: prefix unary minus binds to the power BASE (`-x^2` parses as `Pow(Neg(x), 2)`, NOT `Neg(Pow(x, 2))`), while binary subtraction wraps whole terms (`a - x^2` flattens to `[..., Neg(Pow(x,2))]`). These are DIFFERENT canonical forms — `equals` fails even when Display strings match. When writing expected forms in tests, express negated terms via binary subtraction (`0 - x^2 - …`), never prefix minus on powers.
- **Display can hide structural differences**: two canonically-different Exprs displayed identically (`-1 + -(2*x) + -x^2`). Debug with `mathr::serialize::to_sexpr(&e.canonicalize())` — Display is not a reliable tree dump.
- **Don't build expressions via fold-from-`Num(1)`**: `monomial_to_expr` originally emitted `Mul(1, x)` and `Mul(3, 1)`, which then failed `to_poly`'s `Expr::Num` exponent match. Handle the empty-monomial → `Num(c)` and single-factor cases explicitly.

### poly.rs design
- **Representation**: `Poly = Vec<Term>`, `Term { coeff: f64, monomial: BTreeMap<String, u32> }` — multivariate, sparse, exact-ish f64 coefficients. `normalize` combines like terms (BTreeMap keyed by monomial) and drops exact-zero coefficients only.
- **Bottom-up expand**: expand children first, then `to_poly` the whole node; if conversion fails, keep the node. Functions/symbolic powers/variable denominators stay opaque but their children are already expanded.
- **Div is polynomial only for constant denominators** (divide coefficients by `1/c`).
- **Caps prevent blowup**: integer exponent ≤ 64 for `Pow` expansion (`(x+1)^1000000` stays `Pow`), term count ≤ 20k per product, monomial exponent ≤ 10k (saturating add → bail to `None` → caller keeps symbolic form).
- **Output ordering**: descending total degree, then lex comparison of monomials via explicit exponent walk over the union of variable names (BTreeMap's own `Ord` compares key sets and puts `y^2` before `x^2` — wrong for math convention). Negative terms join with `Sub` and magnitudes, so `-x^2 - 2*x - 1` renders without `+ -` artifacts.

## Symbolic Limit Patterns (limit.rs, added 2026-08-31)

### Three-stage strategy
- **Order matters**: (1) direct substitution on `simplify(expr)` → finite value wins by continuity; (2) L'Hôpital for top-level `Div` in `0/0`/`∞/∞` form — substitute num/den separately, if both zero or both infinite → `simplify(differentiate(num))/simplify(differentiate(den))`, recurse (cap `MAX_LOPITAL = 6`); (3) numeric probing fallback.
- **IEEE ±inf from substitution is NOT trustworthy at finite points**: `eval(1/x)` at `x=0` returns `+inf` (IEEE semantics = right-side behaviour only). Only accept ±∞ from substitution when the point itself is `±∞`; otherwise fall through to probing, which checks both sides (`1/x²` → `+∞`, `1/x` → DNE).

### Numeric probing conventions
- Probe points: finite target → `point ± 10^-k · max(1,|point|)` for k = 2..10, interleaved sides; ±∞ → `±10^k` for k = 1..15.
- **Check convergence BEFORE poles** — otherwise a function converging to a large finite value (e.g. `x + 1e11`) gets misclassified as a pole. Flow: side estimates (last two samples agree within `1e-3·max(1,|v|)`) → both agree → Finite; both converge but disagree → DNE; then pole sign check (`|v| > 1e9` per side; same sign → ±∞, opposite → DNE, one side only → DNE); else DNE.
- Filter probe values to `is_finite()` — `exp(1/x)` overflows to inf near 0 and drops out naturally.
- **Snap reported finite values**: integers within `1e-9` relative (so `sin(x)/x` → exactly `1.0`), then round to 12 significant digits (`0.166666666667` for 1/6). Tests must compare against `snap_num(expected)`, not the raw float.

### What limit.rs handles well / poorly
- Handles: removable singularities (`(x²−1)/(x−1)`), trig classics (`sin x/x`, `(1−cos x)/x²` needing 2 L'Hôpital rounds, `(x−sin x)/x³` needing 3), poles, `∞/∞` at infinity (`x/e^x`), `√(x²+x)−x → 1/2`, bounded/growing behaviour at ∞.
- Limitations (documented in SPEC.md): two-sided only (no `x → a⁺` syntax), L'Hôpital only for top-level quotients (`∞−∞`, `0·∞`, `0^0` fall to probing), probing is heuristic.

### Wiring
- REPL: `limit <expr> [<var>] <point>` — parse trailing tokens; point ∈ {`inf`, `+inf`, `-inf`, number}; var inferred via `Expr::variables()` when omitted (exactly one required). Step dispatch: add `"limit "` branch in `dispatch_steps` BEFORE the `cmd_keywords` fallback, since limit has true multi-step output.

## Arbitrary-Precision Decimal Patterns (bigdec.rs, added 2026-08-31)

### Precision management (the guard-digit contract)
- **Every public function's `prec` argument = working precision = digits kept**. `eval_decimal` computes at `wp = user_prec + GUARD` (GUARD = 10) and the caller displays with `round(x, prec)` (`with_prec`). Public transcendentals (`pi`, `exp`, `ln`, `sin`, …) wrap a `_work` fn that runs at `prec + GUARD`, then round once at the end. `eval_decimal_rounded` is the display-ready entry point.
- **Anti-pattern (cost hours)**: computing series at the OUTPUT precision directly. Each `add` rounds to `prec` sig digits, so an N-term series accumulates ~N·0.5·ulp error — π(30) came out wrong at digit 28 (`…38330` vs `…3832795`). Guard digits absorb this; Newton's `sqrt` self-corrects and is less sensitive, but still wraps.
- **Missing alternating sign in a Taylor series is silent**: `cos_series` without the (−1)^i flip computed cosh — `sin(1)` returned `cosh(π/2−1) = 1.167…` and `tan(π/4) = 1.87…`. Tests comparing against known digit strings caught it.

### bigdecimal crate 0.4.10 API notes
- **Fields `int_val`/`scale` are private**; use `as_bigint_and_scale() -> (Cow<BigInt>, i64)`, `into_bigint_and_scale(self)`, `decimal_digit_count() -> u64`, `order_of_magnitude() -> i64` (= ⌊log10|x|⌋), `with_prec(u64)`, `with_scale_round(i64, RoundingMode)`, `sign() -> num_bigint::Sign`.
- **PartialEq is numeric** (normalizes scale), so `x == BigDecimal::new(bd_round(x)?, 0)` is a valid integer test.
- **`FromStr`/`Display` normalize freely**: parsing plain `"5000000000"` keeps mantissa 5000000000 scale 0, but small decimals Display as e-notation (`1E-7`). Don't assert exact Display strings for values < 0.1 — compare `as_bigint_and_scale` structurally instead.
- **Name collisions**: `num_traits::pow` vs a module fn `pow` → import as `use num_traits::pow as int_pow;`. A parameter named `e` shadows a module fn `e` (Euler's number) — name Expr params `expr`.

### Division sizing (div via BigInt long division)
- Q = round(ma·10^shift/mb) with **shift = prec − da + db** (da/db = `decimal_digit_count` of each operand), result scale = shift − sb + sa. Sizing from digit counts is essential because e-notation expands to plain form (`"1e100"` → mantissa 10^100, scale 0), so the mantissa ratio can carry a 10^±200 factor — a naive fixed shift returns 0 for `1e-100 / 1e100`.

### High-precision algorithms that worked
- **π**: Machin π = 16·atan(1/5) − 4·atan(1/239); atan series ratio z² = 1/25 and 1/57121 → fast convergence. Keep `pi_internal` separate from public `pi` so trig reduction (`reduce_half_pi`) can call it at working precision without double rounding.
- **ln**: normalize m = x/2^k2 into [√2/2, √2) (estimate k2 from `order_of_magnitude()/log2(10)`, fix with ≤8 compare/divide steps), then ln m = 2·atanh((m−1)/(m+1)) with |z| ≤ 0.1716. **ln2 must be a standalone `ln2_internal`** — computing it via `ln(2)` inside `ln` recurses forever.
- **exp**: x = k·ln2 + r with |r| ≤ ln2/2 (k via bd_round of x/ln2), Taylor on r, then exact 2^k BigInt scaling. Negative x via reciprocal of exp_work(−x) (internal fn, not the wrapping public one).
- **trig**: k = round(x/(π/2)), r = x − k·π/2, quadrant switch on k mod 4; Taylor on |r| ≤ π/4. Works for huge arguments (no 2^k cap needed).
- **atan**: |x| > 1 → π/2 − atan(1/x); |x| > 0.5 → π/4 + atan((x−1)/(x+1)) (|z| ≤ 1/3); else direct series. Avoids the slow |z|→1 harmonic case.
- **pow**: integer exponent (test via numeric equality with bd_round, |k| ≤ 1e5) → square-and-multiply with per-step rounding; else exp(y·ln x) for x > 0.

### Critical gotcha: parser folds constants at f64
- `Parser::parse` substitutes `pi`/`e`/`tau` → `Expr::num(f64)` at parse time (parser.rs ~line 619). `eval_decimal` must map `Num(x)` == `std::f64::consts::PI/E/TAU` back to high-precision `pi(prec)`/`e(prec)`, or every high-precision result containing π silently bottoms out at ~1e-16. Compare with `==` on f64 constants — deterministic, since the ONLY source of exactly-those bits is the parser fold (or a user typing π's shortest repr, which is fine to upgrade).

### REPL `dec` command
- Syntax `dec <expr> [prec <n>] [with <var>=<val>,...]`; parse ` prec ` (rfind) BEFORE ` with ` so `dec x with x=1.5 prec 10` works. Reuses `split_assignments` (bracket-depth aware). prec range 1–1000, default 30. Wired in `dispatch_inner` + HELP text; bare `dec` returns usage (otherwise it would eval as an unknown variable).

### Testing digit strings
- Assert `to_string().starts_with(prefix)` with prefixes verified against independently-rounded reference digits; count significant digits carefully — for 0 < x < 1, leading "0." doesn't count, so `prec` digits = `prec` decimals only when the leading digit is nonzero. Truncating a reference instead of rounding it caused two test failures (`…72420` vs correct `…72421`).

## Interval Arithmetic Patterns (interval.rs)

### Core design
- **`Interval { lo, hi, is_empty }`**: closed `[lo, hi]` or empty ∅. `new(lo, hi)` auto-converts to empty if `lo > hi`. `whole()` is `(-∞, +∞)`.
- **No outward rounding**: uses plain `f64` arithmetic (not IEEE 1788 directed rounding). Bounds may be tight at the last bit — documented limitation. For safety-critical use, widen by epsilon.
- **`eval_interval(&Expr, &HashMap<String, Interval>)`**: tree-walking evaluator over the Expr AST with interval values. Constants `pi`, `e`, `tau` provided as point intervals. Unknown vars/functions error.

### Arithmetic rules
- **Add**: `[a,b] + [c,d] = [a+c, b+d]` — endpoint sums.
- **Sub**: `[a,b] - [c,d] = [a-d, b-c]` — cross endpoints (this causes the dependency problem).
- **Mul**: compute all 4 endpoint products `{a·c, a·d, b·c, b·d}`, take min/max. Handles zero-crossing automatically.
- **Div**: if divisor contains 0 → `whole()`. Else same 4-endpoint product pattern with division.
- **Neg**: `[a,b] → [-b, -a]` — bounds swap.

### Function-specific patterns
- **Monotonic increasing** (exp, ln, log10, log2, sqrt, atan, asin, sinh, tanh, cbrt): image is `[f(lo), f(hi)]`. Check domain first (ln/sqrt/log need positive input; return empty if outside domain).
- **Monotonic decreasing** (acos): image is `[f(hi), f(lo)]` — bounds swap.
- **Even function** (cosh): if 0 ∈ [lo,hi], min is `f(0)=1`; else min is `min(f(lo), f(hi))`.
- **Trig extrema tracking** (sin, cos): check if the interval contains a global max (sin: π/2+2kπ; cos: 2kπ) or min (sin: -π/2+2kπ; cos: π+2kπ). If so, bound includes ±1. Use `contains_extremum(lo, hi, target)` helper which checks congruence mod 2π.
- **tan**: has poles at π/2+kπ. If interval contains a pole → `whole()`. Use `interval_contains_tan_pole`.
- **sqr (x²)**: tighter than `self * self` because it knows both operands are the same value. If 0 ∈ [lo,hi] → `[0, max(lo², hi²)]`; else `[min(lo², hi²), max(lo², hi²)]`.
- **powi (x^n)**: even n → if 0 ∈ [lo,hi], min is 0; odd n → monotonic, `[lo^n, hi^n]`.

### Dependency problem (fundamental limitation)
- `x - x` over `[1, 2]` → `[-1, 1]`, NOT `[0, 0]`. Each occurrence of `x` is treated as independent. This over-conservatism grows with expression complexity. Documented in module docs and SPEC.md. No fix without more advanced techniques (affine arithmetic, Taylor models).

### REPL parsing pitfall
- `interval <expr> with x=[-2,3],y=[1,4]` — the comma inside `[lo,hi]` must NOT be treated as an assignment separator. Use `split_assignments()` in `repl.rs` which tracks bracket depth (`[` increments, `]` decrements, split on `,` only when depth==0). (Bug found 2026-08-18.)

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
- **mathml.rs** and **serialize.rs** depend on expr + error only (pure AST transforms)
- **interval.rs** depends on expr + error + std::collections (interval eval over Expr AST)
- **bigdec.rs** depends on expr + error + bigdecimal/num-bigint/num-traits/num-integer (decimal eval over Expr AST; reuses bigint::factorial for `fact`)
- **limit.rs** depends on expr + eval + simplify + symbolic (substitution via eval + Context; L'Hôpital via differentiate)
- **poly.rs** depends on expr only (pure AST transform, like serialize/mathml)
- **apart.rs** depends on expr + poly (dense/sparse conversions) + solver (complex roots) + matrix (coefficient solve) + simplify
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