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

matrix.rs       (standalone, f64 linear algebra: LU, Cholesky, SVD, QR, power iteration,
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
mcp.rs          (standalone, Model Context Protocol server over stdio:
                 JSON-RPC 2.0, 14 tools backed by repl::dispatch_str)
unicode.rs      (standalone, Unicode math input normalization: √ × ÷ π τ −,
                 superscript exponents x² → x^(2); runs before tokenizing)
ode.rs          (standalone, numerical ODE integration)
interpolate.rs  (standalone, Lagrange, Newton, cubic spline, Chebyshev,
                 Legendre, Gauss–Legendre)
special.rs      (standalone, Gamma, Beta, erf, sinc, incomplete gamma,
                 Bessel J_0/J_1/J_n, digamma/trigamma/polygamma,
                 harmonic numbers, Hurwitz/Riemann zeta (Euler–Maclaurin),
                 elliptic integrals K/E/F/E_inc (Carlson RF + AGM))
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
serialize.rs     (standalone, Expr serialization in three interchangeable
                 formats: S-expressions, JSON, RPN. All round-trip via
                 Expr::equals. Hand-rolled JSON parser — no serde dep)
lean.rs          (standalone, Lean 4 export over ℝ — requires Mathlib:
                 literals as exact-ascribed ℝ, integral powers via monoid
                 ^ (ℕ), rpow otherwise, elementary fns mapped to Real.*,
                 log(x,b) → Real.logb b x (arg swap), floor/ceil coerced
                 ↑⌊·⌋/↑⌈·⌉, unmapped fns emitted as plain applications.
                 to_lean (bare term) + to_lean_doc (runnable file with
                 import Mathlib + noncomputable def with var binders))
interval.rs      (standalone, interval arithmetic for rigorous bounds:
                 Interval type with arithmetic, sqr/powi with zero-crossing
                 awareness, elementary functions with monotonicity/extrema
                 tracking, set ops, eval_interval over Expr AST)
bigdec.rs        (standalone, arbitrary-precision decimal arithmetic on the
                 bigdecimal crate: precision-controlled add/sub/mul/div via
                 BigInt long division, sqrt/exp/ln/trig with argument
                 reduction + Taylor/atanh series + 10 guard digits, pi via
                 Machin, eval_decimal over Expr AST. `dec` REPL command)
limit.rs         (standalone, symbolic limits: direct substitution,
                 L'Hôpital's rule for 0/0 and ∞/∞ quotients (recursive,
                 max 6), numeric probing fallback with pole/divergence
                 detection, limits at ±∞. `limit` REPL command + steps)
poly.rs          (standalone, polynomial expansion: multivariate to_poly/
                 poly_to_expr conversions with like-term collection and
                 degree/lex ordering; expand distributes products and
                 integer powers, keeps functions/symbolic powers opaque.
                 `expand` REPL command + steps)
apart.rs         (partial fraction decomposition of rational functions:
                 dense coefficient conversion from poly, long division,
                 Durand–Kerner factorization via polynomial_roots_complex,
                 conjugate pairs → quadratics, root clustering for
                 multiplicities, Gaussian-elimination coefficient solve.
                 `apart` REPL command + steps)

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
7. **Matrix**: `Matrix` — row-major `Vec<f64>`. Gaussian elimination for det/solve/inverse; **LU** with partial pivoting, **Cholesky** `A = L·Lᵀ` for SPD, **SVD** `A = U·Σ·Vᵀ` via Jacobi rotations, **QR** `A = Q·R` via Householder reflections (orthogonal Q, upper-triangular R, least-squares solve), **power iteration** for the dominant eigenpair, **symmetric eigenvalue decomposition** via Householder tridiagonalisation + Wilkinson-shift QR iteration, **Hessenberg decomposition** via Householder reflections, **real Schur decomposition** via shifted QR on Hessenberg form, **Hilbert matrix** construction, **Tikhonov regularised solve** `(AᵀA + λI)x = Aᵀb` for ill-conditioned and rectangular systems
8. **Stats**: Functions on `&[f64]` slices. **Stochastic primitives**: `Rng` (reproducible LCG), `normal`/`exponential` sampling (Box–Muller, inverse-CDF), `normal_pdf`/`normal_cdf`/`exp_pdf`/`exp_cdf`, `moments` (skewness, excess kurtosis), `cumulants`
9. **Number theory**: Functions on `u64` integers (gcd, sieve, totient, Miller–Rabin, CRT, modular inverse, **discrete logarithm**) and on `i64` (extended GCD, **Jacobi symbol**, **continued fractions**, **linear Diophantine solver**)
10. **ODE**: `ode::euler/rk4/rkf45` — operate on closures `F: Fn(f64, f64) → f64`
11. **Taylor**: `taylor::taylor_series` — uses `symbolic::differentiate` + `eval::eval` to compute coefficients
12. **Laurent**: `laurent::laurent_series` — expands `g(x) = (x-a)^k · f(x)` via Taylor, divides by `(x-a)^k` to get principal + analytic parts
13. **Rational**: `rational::Rational` — exact `i64/i64` with GCD reduction, `i128` intermediate arithmetic, parsing from `"n/d"` or decimal strings. `rational::eval_rational` walks `Expr` AST and returns exact `Rational` when all leaves are integers (returns `None` for functions, variables, or non-integer constants)
14. **Notebook**: `notebook::Notebook` — `.mnb` JSON format with cells of TeX/math expressions; `eval_cell`/`eval_all` dispatch through `repl::dispatch_with_ctx` (mutating context for `let`/`fn` persistence); cells have types (`Math`/`Text`) and support reordering/duplication
15. **Web server**: `server::NotebookServer` — minimal HTTP server on `std::net::TcpListener`; serves single-page web UI at `GET /`, REST API at `/api/eval` (with step-by-step `steps` array, shared context, inline plot images as base64), `/api/notebook`, `/api/reset`, `/api/context`
16. **Step-by-step solving**: `repl::dispatch_steps` — returns `Vec<String>` of intermediate steps for `diff`, `solve`, `taylor`, `integrate`, `simplify`, `rat`, `laurent`, and plain eval; tries `rational::eval_rational` first for exact fraction results; falls back to `simplify` when evaluation fails on unbound variables; delegates all other REPL commands (`int`, `romberg`, `fft`, `det`, etc.) to `dispatch_inner`
17. **Symbolic integration**: `symbolic::integrate` — pattern-matches common elementary rules (polynomial, exp, ln, sin/cos/tan/sec, atan, asin). A product factor that does not *mention* the integration variable is treated as a constant and pulled out (`mentions_var` occurrence test, not `is_constant`) — `∫ a·x dx = a·x²/2` for symbolic `a`, not just numeric constants
17b. **`for <var>` suffix**: `repl::strip_for_var` — a trailing `for <name>` on `diff`/`d/dx`/`integrate`/`solve`/`taylor`/`limit` names the free variable of the operation (the unknown to isolate for `solve`). Stripped *before* the trailing-token heuristics (single-letter var, numeric order/point/guess) so `for` never leaks into the expression; keyword matched case-insensitively, name must be letters-only, suffix must be final. Wins over the trailing single-letter heuristic. Applied on both the `dispatch_steps` and `dispatch_inner` paths
17c. **Negation of a subtraction**: `simplify::preflatten` flips `-(a - b)` to `b - a` so a negated constant never lingers as a `Neg` wrapping a `Sub` (keeps linear-isolate answers readable: `y = (6 - 3*x)/2`, not `y = -(3*x - 6)/2`)
18. **Gradient**: `symbolic::gradient` — collects all free variables and computes partial derivatives for each
18b. **Implicit differentiation**: `symbolic::implicit_diff` / `idiff_steps_str` — `d(dep)/d(indep) = −F_indep/F_dep` for `F=0`; REPL `idiff`
18c. **Higher-order differentiation**: `symbolic::differentiate_n` — repeated `differentiate` (order capped at 64); REPL `diff <expr> [var] [order]` and `d/dx <expr> [order]`
18d. **Vector ops**: `vector::dot` / `vector::cross` / `vector::cross3` — numeric dot and cross products from comma-separated REPL lists; `dot` / `cross` commands
18e. **MCP server**: `mcp::serve_stdio` — Model Context Protocol over **stdio** (newline-delimited JSON-RPC 2.0, protocol version 2025-06-18). Handshake (`initialize` → `notifications/initialized`), `tools/list`, `tools/call`, `ping`. 14 tools (`evaluate`, `differentiate`, `integrate`, `solve`, `nsolve`, `simplify`, `taylor`, `limit`, `poly-roots`, `number-theory`, `stats`, `matrix`, `minimize`, `apart`, `expand`), each a thin wrapper that builds a command string and runs it through `repl::dispatch_str` — so the MCP surface cannot drift from the CLI. Hand-rolled JSON layer (no `serde` dep), matching the `serialize.rs` convention. Tool *execution* failures return a successful response with `isError: true` (the MCP convention); only protocol problems (unknown tool/method, bad JSON) use JSON-RPC error codes. Launched with `mathr mcp`
18f. **Unicode input normalization**: `unicode::normalize` runs in `Parser::parse` **before** tokenizing, rewriting Unicode math into the ASCII the parser understands: `×`→`*`, `÷`→`/`, `−`(U+2212)/en-dash/em-dash→`-`, `π`→`pi`, `τ`→`tau`, `∞`→`inf`, and superscript runs → explicit powers (`x²`→`x^(2)`, `x⁻¹`→`x^(-1)`). The root signs need argument delimiting, not plain substitution: `√2`→`sqrt(2)` (not `sqrt2`, which would read as one identifier), so the root branch captures its argument and normalizes it separately. Unmapped glyphs are passed through unchanged so the parser reports the *real* error rather than a mangled one
18g. **Selective exact function folding**: `simplify` folds `sqrt`/`cbrt`/`log2`/`log10`/`abs` **only when the result is provably exact** — `sqrt(16)`→`4` but `sqrt(2)` stays `sqrt(2)`, because an irrational value has a shorter exact symbolic form and collapsing it to a float discards exactness. Exactness is decided by `Rational::exact_sqrt`/`exact_cbrt`/`exact_log2`/`exact_log10` (exact roots of rationals, each cross-checked exactly). The folding is a **separate first pass** (`fold_exact_functions`) that runs *before* `collect`, because `collect` turns `1/4` into the float `0.25` — after which the exactness of `sqrt(1/4)` can no longer be established. `try_eval_const` deliberately does **not** fold exact-only functions, so the final constant-folding pass cannot re-flatten them
18h. **All-roots sign scan**: `solver::find_all_roots` + `ScanOptions` — samples an interval (default `[-30, 30]`, 2400 points), bisects every sign change, and returns all roots sorted and de-duplicated. `nsolve <expr> [var] [lo hi]` REPL command (plus step-by-step output). **Pole rejection is the critical property**: a sign flip across a pole is indistinguishable from a sign flip across a root, so a bisected bracket is accepted only when `f` lands near zero at the bisection point — `1/(x-1)` reports *no real roots*, and `tan(x)` reports only the multiples of π, not the poles at ±π/2. Undefined points (NaN) break the scan instead of producing a spurious bracket. Bounds are stripped *before* the `for <var>` suffix, so `nsolve sin(t) for t -7 7` parses correctly
18i. **Exact radical quadratic roots**: `qsolve` keeps irrational quadratic roots **exact** rather than collapsing them to decimals — `x²−2` → `sqrt(2)`, `x²−8` → `2*sqrt(2)`, `(1/2)x²−(3/4)x+1/8` → `3/4 + sqrt(5)/4`. The coefficients are read as exact rationals (`exact_coeff`, round-trip verified), the discriminant `b²−4ac` is computed in exact `Rational` arithmetic, and `split_radical` pulls out the largest square *factor* of the discriminant (`√8 → 2√2`, `√(5/16) → √5/4`) so the inner radicand is square-free. `build_radical_root` writes the result as `p ± q·√d` with every coefficient kept rational. A perfect-square discriminant still yields plain rationals; a negative one still yields complex `i`-literal roots
18j. **Exact cover-up partial fractions**: `apart` uses the **cover-up rule** `Aᵢ = N(rᵢ)/D′(rᵢ)` (with `D′(rᵢ) = lead·Π_{j≠i}(rᵢ−rⱼ)`) when the denominator splits into *distinct rational* linear factors — no linear system is needed, and the coefficients come out as exact rationals (`1/(x²−1)` → `1/2/(x−1) − 1/2/(x+1)`). The roots are found with the rational-root theorem in exact `Rational` arithmetic (candidates `p/q` with `p` dividing the constant term and `q` the leading coefficient, **both signs** of `p`). The path is declined for repeated roots, complex roots, or irrational roots, which fall back to the numeric factorization + linear-system path. `apart`/`apart_steps` no longer re-simplify the result — `simplify` folds rationals to floats, which would undo the exactness
18k. **Shared bound parsing**: `repl::split_bounds` gives `int`, `romberg` and `minimize` one bound parser that accepts either the positional form (`<expr> <a> <b>`) or an explicit `from <a> to <b>` suffix whose bounds are full *expressions* (so `pi`, `pi/2`, `2*pi` work and may contain spaces). The suffix form is preferred when present because it is unambiguous; the positional form still takes the **last two tokens** as the bounds, which keeps the trailing-token ambiguity class retired. A malformed suffix is a clear error rather than being read as part of the expression. This also fixed an inconsistency: `romberg sin(x) 0 pi` used to fail with "bad b" while `int sin(x) 0 pi` worked
19. **Interpolate**: `interpolate::lagrange_interp/newton_interp/CubicSpline/chebyshev_*/legendre_*/gauss_legendre` — point-wise polynomial, smooth C²-continuous cubic, Chebyshev series, Legendre polynomials, Gauss–Legendre quadrature
20. **Special**: `special::gamma/erf/sinc/bessel_j0/j1/jn` — Lanczos, continued-fraction, Maclaurin-series, and asymptotic-form approximations
21. **Plot**: `plot::plot_function/multi/scatter` — evaluates `Expr` over a range, renders PNG via `plotters`
22. **Serialize**: `serialize::to_sexpr/from_sexpr`, `to_json/from_json`, `to_rpn/from_rpn` — three interchangeable textual encodings of `Expr`. All round-trip via `Expr::equals`. JSON uses a hand-rolled recursive-descent parser (no serde dependency, matching the crate convention). `serialize <fmt> <expr>` and `serialize <fmt> import <text>` REPL commands.
23. **Interval**: `interval::Interval` — closed `[lo, hi]` bounds with outward-conservative arithmetic (+, -, *, /, neg), `sqr`/`powi` with zero-crossing awareness, elementary functions (sin, cos, tan, exp, ln, sqrt, abs, asin, acos, atan, sinh, cosh, tanh, cbrt) with monotonicity and extrema tracking for trig, set operations (intersect, hull, contains, overlaps), and `eval_interval` for rigorous bounds over the `Expr` AST. `interval <expr> with <var>=[lo,hi],...` REPL command.
24. **Arbitrary-precision decimals**: `bigdec` — precision-controlled arithmetic over `bigdecimal::BigDecimal` (`add`/`sub`/`mul` round to N significant digits; `div` uses BigInt long division sized from operand digit counts), transcendentals (`pi`, `e`, `sqrt`, `exp`, `ln`, `log2`, `log10`, `sin`, `cos`, `tan`, `atan`, `asin`, `acos`, `sinh`, `cosh`, `tanh`, `pow`) computed at `prec + 10` guard digits and correctly rounded, `eval_decimal`/`eval_decimal_rounded` over the `Expr` AST with variable bindings. `dec <expr> [prec <n>] [with <var>=<val>,...]` REPL command.
25. **Limits**: `limit::limit`/`limit::limit_steps` — three-stage strategy: direct substitution on the simplified expression, L'Hôpital's rule for `0/0`/`∞/∞` quotients (simplify → differentiate → recurse, max 6), numeric probing near the target (both sides, geometric steps) with convergence/pole/oscillation classification. Results are `LimitValue::{Finite, PosInfinity, NegInfinity, DoesNotExist}`; finite values are snapped to integers/12 significant digits. `limit <expr> [<var>] <point>` REPL command; step-by-step via `dispatch_steps`.
26. **Polynomial expansion**: `poly::expand` — bottom-up: children expand first, then a node converts to a dense multivariate polynomial (`Vec<Term>`, `coeff × var^exp` map) when possible and back to a collected sum. Safety caps: exponent ≤ 64, term count ≤ 20k, monomial exponent ≤ 10k. `to_poly`/`poly_to_expr` are public building blocks for future partial fractions. `expand <expr>` REPL command; step-by-step via `dispatch_steps`.
27. **Partial fractions**: `apart::apart` — long division (`N = Q·D + R`), numeric factorization of `D` via `solver::polynomial_roots_complex` (conjugate pairs → quadratics `x² + px + q`; root clustering → multiplicities), then a square K×K linear system (columns = `D/F^k`, and `x·D/F^k` for quadratics) solved with `Matrix::solve`. Sign-aware output joins and coefficient snapping. `apart <expr> [<var>]` REPL command; step-by-step via `dispatch_steps`.
28. **Condition number & nullspace**: `Matrix::condition_number` (2-norm κ = σ_max/σ_min from `svd()`; `inf` when singular) and `Matrix::nullspace` (orthonormal basis from the eigendecomposition of `Aᵀ·A` via `symmetric_eig`, relative tolerance — works for wide matrices where `svd()`'s V is not the full right-vector basis). `cond`/`null` REPL commands.
29. **Sparse matrices & iterative solvers**: `sparse::Csr`/`sparse::Csc` — compressed-row/column storage built from coordinate triplets (duplicates summed, canonical ascending-column order) or dense matrices; `matvec`, `transpose` (scatter into the opposite layout), row-wise sparse × sparse `multiply`, `conjugate_gradient` (structurally verified symmetry, non-positive-curvature breakdown detection, relative-residual stopping), `conjugate_gradient_jacobi` (diagonal preconditioning through a shared preconditioned core), and `bicgstab` (general nonsymmetric systems, four explicit breakdown guards). `spcg [jacobi]`/`spbicg` REPL commands; `sparse_demo` example.
30. **B-spline & Hermite interpolation**: `bspline` — `basis_function` (Cox–de Boor recursion; right-endpoint evaluated as left-hand limit), `BSpline` (de Boor evaluation with an affine `t0/t_scale` parameter map so data-x evaluation matches interpolants built on normalized parameters), `cubic_bspline_interp` (de Boor knot averaging → Schoenberg–Whitney-guaranteed banded collocation, dense solve, degree reduction for n < 4); `interpolate::CubicHermite` — piecewise cubic Hermite with explicit slopes. `bspline`/`hermite` REPL commands; `bspline_demo` example.
31. **Complex evaluation**: `ceval::eval_complex` — tree-walks the `Expr` AST over `Complex<f64>` with zero parser changes (`Var("i")` is the imaginary unit; implicit multiplication makes `1 + 2i` natural; other variables bind real values from `Context`). Principal-branch `ln`/`sqrt`/powers/inverse trig; exact repeated multiplication for integer powers, `exp(w·ln z)` otherwise. IEEE −0.0 canonicalized at the root (`−0.0` breaks `atan2` sign selection, flipping `sqrt(-1)` to `−i`). `cval <expr> [with <var>=<val>,...]` REPL command; `format_complex` snaps noise for display.
32. **Symbolic polynomial solve (deg ≤ 4)**: `qsolve::solve_symbolic` — `simplify` → `to_poly` ascending coefficients; degree ≤ 2 uses the quadratic formula (exact rationals when the discriminant is a perfect square; complex `i`-literal roots via `ceval`; else 12-sig snap). Degree 3–4 with integer coefficients: rational-root candidates + exact `Rational` synthetic deflation, then finish at degree ≤ 2; otherwise Durand–Kerner numeric roots with real-root certification. One inferred variable (no trailing-var ambiguity). `qsolve <expr> [= rhs]` REPL command; step-by-step via `dispatch_steps`.
33. **Distributions & hypothesis tests**: `dists` — Student-t/chi²/F/binomial/Poisson/uniform pdf/cdf/pmf composed from `special::incomplete_gamma_p` and `beta_inc` (continued fraction), quantiles via generic bracketed bisection (`invert_cdf`), hypothesis tests returning statistic/df/p-value (one-sample/Welch/paired t, chi-squared, one-way ANOVA, Mann–Whitney U, Wilcoxon signed-rank, Kruskal–Wallis — all tie-corrected), Spearman rank correlation, seeded percentile bootstrap confidence intervals. `dist`/`qtile`/`ttest`/`chitest`/`anova`/`mwu`/`wilcoxon`/`kw`/`boot` REPL commands.
34. **Nonlinear least squares**: `curvefit::curve_fit` — Levenberg–Marquardt with numeric central-difference Jacobian, Marquardt diagonal scaling, adaptive damping sub-loop, parameter standard errors from σ²(JᵀJ)⁻¹; takes `FnMut` models (Expr-evaluating closures). `fit` REPL command; `curvefit_demo` example.
35. **Logistic regression**: `logit::logistic_regression` — IRLS/Newton on the log-likelihood with a 1e-10 ridge (separation stays well-defined: coefficients diverge honestly with huge Wald standard errors), overflow-safe sigmoid and log-likelihood, `predict_proba`. `logit` REPL command; `logit_demo` example.
36. **Stiff ODE solvers**: `stiff` — backward Euler (L-stable), implicit trapezoidal (A-stable, 2nd order), and BDF2 (A-stable, trapezoidal startup step), all sharing one Newton core (`G(y) = y − base − hb·f(t_new, y) = 0`) with central-difference numeric Jacobians solved via `Matrix::solve`. `stiff_demo` example.
37. **Monotone interpolation (PCHIP)**: `pchip::Pchip` — Fritsch–Carlson piecewise cubic Hermite with non-uniform-grid weighted-harmonic-mean slopes, zero slope at local extrema, endpoint clipping; never overshoots the knot range (the THE property), n = 2 degenerates to the exact line. `pchip` REPL command.
38. **Derivative-free optimization**: `optim` — `golden_section` (1-D unimodal bracket, φ-shrink, reuses the leftover interior point) and `nelder_mead` (N-D simplex with reflection/expansion/contraction/shrink, convergence on simplex diameter). `minimize <expr> [var] a b` REPL command.
39. **Summation & product**: `sumprod::summation/product` — Σ/∏ over inclusive integer bounds with three evaluation tiers: symbolic closed forms (var-independent → `expr·count`; `var` → arithmetic series; `var²`/`var³` → Faulhaber; `c^var` → geometric — all O(1) i128-exact where the value fits f64 integers), an exact rational term loop (`eval_rational` accumulation, capped and cross-checked against the f64 pass because i64-backed Rationals silently overflow on denominator-heavy sums), and a plain f64 loop. Empty ranges: Σ = 0, ∏ = 1. `sum <expr> [<var>] <a> <b>` / `prod <expr> [<var>] <a> <b>` REPL commands (minimize-style trailing bounds).

## Design Decisions

- **No external math crates**: FFT, complex numbers, matrix ops, stats, special functions, and number theory are all implemented from scratch. Only `plotters` is used for rendering.
- **Expr as the universal AST**: One enum serves parsing, evaluation, symbolic differentiation, symbolic integration, simplification, and Taylor series.
- **Closures for numerical methods**: Solvers, integrators, and ODE steppers accept `Fn` closures, making them composable with the evaluator.
- **Inline tests**: Each module has `#[cfg(test)] mod tests` — no separate test files. **689 tests** total (plus 213 integration tests).
- **Error handling**: `MathError` (via `thiserror`) in the library, `anyhow` in the binary.
- **Prelude**: `mathr::prelude` re-exports the most common types for ergonomic `use mathr::prelude::*;`.
- **Numerical recipes**: Where possible (Gamma, sinc, Bessel, incomplete gamma) we use well-tested A&S or NR polynomial approximations; otherwise we use direct series / closed-form recurrences.
- **Matrix decomposition family**: Gaussian elimination for inverse/determinant; LU with partial pivoting for fast solves; Cholesky for symmetric positive-definite matrices; SVD via Jacobi rotations for rank/reconditioning; QR via Householder reflections for least-squares; power iteration for the dominant eigenpair; symmetric eigenvalue decomposition via Householder tridiagonalisation + Wilkinson-shift QR iteration; Hessenberg decomposition via Householder reflections; real Schur decomposition via shifted QR on Hessenberg form with 2×2 block handling for complex conjugate eigenvalue pairs.

## Deployment

The crate produces both a library (`mathr`) and a binary (`mathr`). The binary is a thin CLI wrapper around library functions. The `notebook` subcommand starts a minimal HTTP server serving a Jupyter-like web UI with KaTeX math rendering, step-by-step solving, exact fraction arithmetic, shared context across cells (variables/functions persist), cell types (math/text with Markdown rendering), inline plots (base64 PNG), execution counters, cell management (reorder, duplicate, toggle type), and a CSS-responsive layout for tablet/phone viewports.