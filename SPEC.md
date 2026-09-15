# Specification

## CLI Subcommands

| Command | Syntax | Description |
|---------|--------|-------------|
| `eval` | `mathr eval <expr> [--set name=value]` | Evaluate expression |
| `diff` | `mathr diff <expr> [--var x] [--simplify]` | Symbolic derivative (partial if multivariate) |
| `simplify` | `mathr simplify <expr>` | Constant-fold & simplify |
| `integrate` | `mathr integrate <expr> [--var x]` | Symbolic indefinite integral |
| `integrate-num` | `mathr integrate-num <expr> a b [--var x] [--n N] [--adaptive] [--romberg LEVELS]` | Numerical integral |
| `solve` | `mathr solve <expr> [--var x] [--guess 0] [--bisect A B] [--max-iter 100] [--tol 1e-10]` | Root finding |
| `solve-system` | `mathr solve-system <sys> [--guess x0,y0,...]` | Newton's method for nonlinear systems |
| `poly-roots` | `mathr poly-roots <coeffs...>` | Polynomial roots (Durand–Kerner) |
| `isolate-roots` | `mathr isolate-roots <ints...>` | Real root isolation (VAS, integer coefficients) |
| `plot` | `mathr plot <expr> [-o path] [--var x] [--a -τ] [--b τ] [--samples 800]` | PNG plot |
| `fft` | `mathr fft <samples...> [--complex] [--inverse] [--magnitude] [--power]` | FFT |
| `conv` | `mathr conv <a...> <b...>` | FFT convolution of two signals |
| `taylor` | `mathr taylor <expr> [--var x] [--around 0] [--order 5]` | Taylor series |
| `laurent` | `mathr laurent <expr> [a] [pole_order] [n_positive]` | Laurent series around a pole |
| `rat` | `mathr rat <a> <op> <b>` | Exact rational arithmetic |
| `notebook` | `mathr notebook [file.mnb] [port]` | Web notebook UI (Jupyter-like) |
| `fourier` | `mathr fourier <expr> <L> <N> [x]` | Fourier series on [-L, L] with N terms |
| `mc` | `mathr mc <expr> <a> <b> <N> [seed]` | Monte Carlo integral over [a, b] |
| `sample` | `mathr sample <dist> <params...> <N> [seed]` | Random sampling (uniform/normal/exponential) |
| `dist` | `mathr dist <dist> <x> <params...>` | PDF and CDF (`normal`/`exp`/`uniform`/`t`/`chi2`/`f`/`binom`/`poisson`) |
| `qtile` | `mathr qtile <dist> <p> <params...>` | Quantile / critical value (`normal [mu sigma]`/`uniform`/`t`/`chi2`/`f`) |
| `logit` | `mathr logit <y0 1...> with <x1...> [\| <x2...>]` | Logistic regression (IRLS) with Wald standard errors and fitted probabilities |
| `fit` | `mathr fit <expr> [in <var>] [with p=v,...] x1 y1 ...` | Levenberg–Marquardt nonlinear least-squares curve fit |
| `ttest` | `mathr ttest <mu> <values...>` | One-sample t-test (two-sided) |
| `ttest` | `mathr ttest <values...> \| <values...>` | Two-sample Welch t-test |
| `ttest` | `mathr ttest paired <values...> \| <values...>` | Paired t-test |
| `chitest` | `mathr chitest <observed...> [\| <expected...>]` | Chi-squared test (uniform expectation by default) |
| `anova` | `mathr anova <group1...> \| <group2...> [\| ...]` | One-way ANOVA F-test |
| `mwu` | `mathr mwu <group1...> \| <group2...>` | Mann–Whitney U test (nonparametric two-sample) |
| `wilcoxon` | `mathr wilcoxon <before...> \| <after...>` | Wilcoxon signed-rank test (paired) |
| `kw` | `mathr kw <group1...> \| <group2...> [\| ...]` | Kruskal–Wallis H test (nonparametric ANOVA) |
| `boot` | `mathr boot <mean\|median> <values...> [seed <n>]` | Bootstrap 95% confidence interval (10000 resamples, default seed 42) |
| `stats` | `mathr stats <data...>` | Descriptive statistics |
| `matrix` | `mathr matrix <op> <rows...>` | `lu`/`qr`/`cholesky`/`svd`/`eig`/`symlig`/`hessenberg`/`schur`/`rank`/`cond`/`null`/`det`/`solve` |
| `tikhonov` | `mathr tikhonov <rows...> \| <b...> <lambda>` | Tikhonov-regularised solve |
| `spcg` | `mathr spcg [jacobi] <rows...> \| <b...>` | Sparse conjugate-gradient solve of `Ax = b` (A must be symmetric positive-definite; `jacobi` = diagonal preconditioning) |
| `spbicg` | `mathr spbicg [ilu] <rows...> \| <b...>` | Sparse BiCGStab solve (nonsymmetric OK; `ilu` = ILU(0) preconditioning, default Jacobi) |
| `interp` | `mathr interp <op> ...` | `lagrange`/`newton`/`spline`/`chebyshev`/`legendre` |
| `pchip` | `mathr pchip x1 y1 x2 y2 ... x_at` | Monotone piecewise cubic Hermite interpolation (no overshoot) |
| `bspline` | `mathr bspline x1 y1 x2 y2 ... x_at` | Clamped B-spline interpolant (cubic for 4+ points) |
| `hermite` | `mathr hermite x1 y1 d1 x2 y2 d2 ... x_at` | Piecewise cubic Hermite with given slopes |
| `minimize` | `mathr minimize <expr> [var] a b` | Golden-section minimum of expr on `[a, b]` |
| `gcd` | `mathr gcd <n1> <n2> [...]` | GCD of integers |
| `lcm` | `mathr lcm <n1> <n2> [...]` | LCM of integers |
| `is-prime` | `mathr is-prime <n>` | Primality test |
| `factor` | `mathr factor <n>` | Prime factorization |
| `fib` | `mathr fib <n>` | nth Fibonacci number |
| `binom` | `mathr binom <n> <k>` | Binomial coefficient C(n,k) |
| `fact` | `mathr fact <n>` | Factorial n! |
| `mr-prime` | `mathr mr-prime <n> [--rounds 20]` | Miller–Rabin primality test |
| `jacobi` | `mathr jacobi <a> <n>` | Jacobi symbol (a/n) |
| `cf` | `mathr cf <p> <q>` | Continued fraction of p/q |
| `diophantine` | `mathr diophantine <a> <b> <c>` | Solve a·x + b·y = c |
| `dlog` | `mathr dlog <g> <h> <p>` | Discrete logarithm `g^x ≡ h (mod p)` |
| `special` | `mathr special <op> <x>` | `gamma`/`erf`/`erfc`/`sinc`/`bessel_j0`/`bessel_j1`/`bessel_j` |
| `fast` | `mathr fast <func> <x> [y]` | Chebyshev fast approx (`sin`/`cos`/`tan`/`exp`/`log`/`sqrt`/`pow`) |
| `big` | `mathr big <op> <args>` | Big integer ops for inputs > u64::MAX (`prime`/`factor`/`gcd`/`lcm`/`modpow`/`totient`). `fact`/`fib`/`binom` auto-upgrade on overflow. |
| `ad` | `mathr ad <expr> at <var>=<val>` | Automatic differentiation (dual numbers) — returns `f(x)` and `f'(x)` |
| `ad` | `mathr ad grad <expr> with <var>=<val>,...` | Gradient of a multivariate expression |
| `ad` | `mathr ad jacobian <f1>, <f2>, ... with <var>=<val>,...` | Jacobian matrix of a system |
| `repl` | `mathr repl` | Interactive REPL |

## Expression Grammar

```
expr   := term (('+' | '-') term)*
term   := factor (('*' | '/' | 'mod') factor)*  -- includes implicit multiplication
factor := unary ('!')* ('^' factor)?            -- postfix factorial, then right-assoc power
unary  := ('+' | '-')? atom
atom   := number | ident | ident '(' args ')' | '(' expr ')' | '|' expr '|'
args   := expr (',' expr)*
```

### Implicit Multiplication

`2x` → `2*x`, `3(x+1)` → `3*(x+1)`, `(x)(y)` → `x*y`

### Postfix Factorial

`n!` → `factorial(n)`, e.g. `5!` = 120, `(2+3)!` = 120, `3!^2` = 36

Factorial binds tighter than `^`: `n!^2` = `(n!)^2`

### Absolute Value

`|x|` → `abs(x)`, e.g. `|-5|` = 5, `|sin(pi)|` = 0

### Infix Modulo

`a mod b` → `mod(a, b)`, e.g. `7 mod 3` = 1

### Function-Call Notation

These number-theory functions are available as function calls in expressions:

| Notation | Example | Result |
|----------|---------|--------|
| `gcd(a, b)` | `gcd(12, 8)` | 4 |
| `lcm(a, b)` | `lcm(4, 6)` | 12 |
| `C(n, k)` | `C(5, 2)` | 10 |
| `factorial(n)` | `factorial(5)` | 120 |

### TeX Notation

| TeX | Equivalent | Example |
|-----|------------|---------|
| `\binom{n}{k}` | `C(n, k)` | `\binom{5}{2}` = 10 |
| `\gcd(a, b)` | `gcd(a, b)` | `\gcd(12, 8)` = 4 |
| `\lcm(a, b)` | `lcm(a, b)` | `\lcm(4, 6)` = 12 |
| `\frac{a}{b}` | `a / b` | `\frac{1}{2}` = 0.5 |
| `\sqrt{x}` | `sqrt(x)` | `\sqrt{4}` = 2 |

### MathML

W3C Presentation MathML is supported for both export and import.

**Export** (`mathml <expr>`):
```
mathr> mathml x^2 + 1
<math xmlns="http://www.w3.org/1998/Math/MathML"><mrow><msup><mi>x</mi><mn>2</mn></msup><mo>+</mo><mn>1</mn></mrow></math>
```

**Import** (`mathml import <MathML>`):
```
mathr> mathml import <mfrac><mn>1</mn><mn>2</mn></mfrac>
1/2
```

Supported MathML elements: `<mn>`, `<mi>`, `<mo>`, `<mrow>`, `<mfrac>`, `<msup>`, `<msub>`, `<msqrt>`, `<mroot>`, `<mtext>`, `<mstyle>`, `<mfenced>`, `<math>`

### Serialization

`Expr` can be serialized to and from three interchangeable textual formats via the `serialize` module. All three round-trip: `from_*(to_*(e)).equals(e)`.

**S-expressions** (Lisp-like prefix notation):
```
mathr> serialize sexpr 2*x + 1
(add (mul (num 2) (var x)) (num 1))
mathr> serialize sexpr import (add (mul (num 2) (var x)) (num 1))
2*x + 1
```
Node forms: `(num <n>)`, `(var <name>)`, `(neg <e>)`, `(add|sub|mul|div|pow <a> <b>)`, `(func <name> <arg>...)`. Non-finite numbers use `NaN`, `inf`, `-inf`.

**JSON** (nested objects):
```
mathr> serialize json x^2
{"t":"pow","a":{"t":"var","v":"x"},"b":{"t":"num","v":2}}
mathr> serialize json import {"t":"pow","a":{"t":"var","v":"x"},"b":{"t":"num","v":2}}
x^2
```
Schema: `{"t":"num","v":<n>}` (or `"v":"NaN"|"inf"|"-inf"` for non-finite), `{"t":"var","v":"<name>"}`, `{"t":"neg","e":{...}}`, `{"t":"add|sub|mul|div|pow","a":{...},"b":{...}}`, `{"t":"func","n":"<name>","a":[{...},...]}`.

**RPN** (postfix, space-separated):
```
mathr> serialize rpn 2*x + 1
2 x * 1 +
mathr> serialize rpn import 2 x * 1 +
2*x + 1
```
Operators: `+ - * / ^` (binary), `neg` (unary). Functions use a `<name>:<arity>` call token, e.g. `x sin:1` for `sin(x)`, `x 2 pow:2` for `pow(x, 2)`.

### Interval Arithmetic

Interval arithmetic computes guaranteed bounds on a function's output over a range of inputs. Instead of a single value, each variable is assigned an interval `[lo, hi]`, and operations propagate worst-case bounds through the expression.

**REPL** (`interval <expr> with <var>=[lo,hi],...`):
```
mathr> interval x^2 + 1 with x=[-2,3]
[1, 10]
mathr> interval sin(x) with x=[0,6.283185307179586]
[-1, 1]
mathr> interval x*y with x=[1,2],y=[3,4]
[3, 8]
```

**Limitations**:
- No IEEE 1788 outward rounding — bounds use plain `f64` and may be tight at the last bit. Widen by a small epsilon for safety-critical use.
- **Dependency problem**: `x - x` over `[1, 2]` yields `[-1, 1]`, not `[0, 0]`, because each occurrence of `x` is treated independently. This is fundamental to interval arithmetic.
- Division by an interval containing zero returns the whole real line `[-∞, ∞]`.

**Supported functions**: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `exp`, `ln`, `log`, `log10`, `log2`, `sqrt`, `abs`, `sqr`, `cbrt`, `floor`, `ceil`, `round`, `sign`, `fract`, `min`, `max`, `pow`, `mod`. Trig functions track global extrema (±1) when the input interval spans a peak or trough.

### Arbitrary-Precision Decimals

The `bigdec` module evaluates expressions with exact decimal arithmetic at a chosen number of significant digits (default 30, max 1000). Constants `pi`, `e`, and `tau` are computed at full target precision (not the f64 approximations used by normal evaluation).

**REPL** (`dec <expr> [prec <n>] [with <var>=<val>,...]`):
```
mathr> dec pi prec 50
3.1415926535897932384626433832795028841971693993751
mathr> dec sqrt(2) prec 30
1.414213562373095048801688724210
mathr> dec 1/3 prec 10
0.3333333333
mathr> dec x*2 + 1 with x=1.5 prec 10
4.000000000
```

**Supported functions**: `sqrt`, `cbrt`, `exp`, `ln`, `log`, `log2`, `log10`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `abs`, `floor`, `ceil`, `round`, `sign`, `min`, `max`, `pow(a, b)`, `fact(n)`.

**How it works**: all intermediate steps keep 10 guard digits beyond the requested precision (MPFR-style); π comes from Machin's formula (16·atan(1/5) − 4·atan(1/239)), exp/ln use argument reduction plus Taylor/atanh series, and trig uses π/2 quadrant reduction with Taylor series.

**Limitations**:
- Decimal exponents beyond ±1,000,000 (i.e. 10^±1000000) are rejected.
- Results are correctly rounded to `prec` significant digits but are not interval-certified (no outward rounding).

### Complex Evaluation

`cval` evaluates expressions over the complex numbers. The parser needs no changes: `i` parses as a variable and is interpreted as the imaginary unit (implicit multiplication makes `1 + 2i` natural); other variables bind real values via `with`. Integer powers are exact; `ln`, `sqrt`, powers, and inverse trigonometric functions use the principal branch (`arg` in `(-pi, pi]`).

```
mathr> cval exp(i*pi)
-1
mathr> cval (1 + 2i) * (3 - i)
5 + 5i
mathr> cval sqrt(-1)
i
mathr> cval ln(-1)
3.1415926536i
mathr> cval x*i - y with x=2, y=-1
1 + 2i
```

Supported functions: `sqrt`, `exp`, `ln`, `log`, `log10`, `log2`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `sinh`, `cosh`, `tanh`, `abs`. Tiny components (below `1e-12`) are snapped away in output.

### Symbolic Quadratic Solve

`qsolve` solves `expr = 0` (or `lhs = rhs`) symbolically for a single inferred variable when the equation is linear or quadratic after expansion. Discriminant classification: perfect-square discriminant → exact integer roots; zero → double root; negative → complex-conjugate roots written with the `i` literal (evaluate them with `cval`); otherwise roots are snapped to 12 significant digits.

```
mathr> qsolve x^2 - 4
x = 2, x = -2
mathr> qsolve x^2 = 2x + 3
x = 3, x = -1
mathr> qsolve x^2 - 2x + 5
x = 1 + 2*i, x = 1 - 2*i
mathr> qsolve 2x - 6
x = 3
```

Limitations: degree > 2 is rejected (use `poly-roots` for numeric roots), non-polynomial expressions are rejected (use the numeric `solve`), and irrational roots are decimal approximations.

### Summation & Product

`sum` and `prod` evaluate Σ/∏ over inclusive integer bounds. Bounds come from the last two tokens; the variable is inferred (or given as a single-letter token before the bounds). Evaluation tiers: O(1) closed forms (var-independent → `expr·count`, `var` → arithmetic series, `var^2`/`var^3` → Faulhaber, `c^var` → geometric), an exact rational term loop (result shown as a fraction, e.g. `Σ 1/k, k=1..10 = 7381/2520`), else f64. Empty range: Σ = 0, ∏ = 1. Max 1,000,000 terms.

```
mathr> sum k 1 100
5050
mathr> sum 1/k 1 10
7381/2520
mathr> sum 2^x 0 10
2047
mathr> prod x 1 5
120
```

### Limits

The `limit` module computes two-sided limits `lim x→a f(x)` for finite points and for `x → ±∞`, using a three-stage strategy:

1. **Direct substitution** on the simplified expression (continuity).
2. **L'Hôpital's rule** for quotients in `0/0` or `∞/∞` form — differentiates numerator and denominator and recurses (up to 6 times).
3. **Numeric probing** — evaluates both sides approaching the target; detects finite convergence, poles (`±∞`), opposite-side divergence, and oscillation (`does not exist`).

**REPL** (`limit <expr> [<var>] <point>` — the variable may be omitted when the expression has exactly one):
```
mathr> limit sin(x)/x x 0
lim sin(x)/x as x → 0 = 1
mathr> limit (x^2 - 1)/(x - 1) 1
lim (x^2 - 1)/(x - 1) as x → 1 = 2
mathr> limit 1/x^2 0
lim 1/x^2 as x → 0 = +∞
mathr> limit 1/x 0
lim 1/x as x → 0 = does not exist
mathr> limit (2*x + 1)/(x + 5) inf
lim (2*x + 1)/(x + 5) as x → +∞ = 2
mathr> limit x/exp(x) x inf
lim x/exp(x) as x → +∞ = 0
```

Points: any number, or `inf`/`-inf`. The notebook step-by-step view shows substitution results, L'Hôpital applications, and the final verdict.

**Limitations**:
- Two-sided limits only (no one-sided `x → a⁺` syntax; sides are analysed internally).
- L'Hôpital applies to top-level quotients only; other indeterminate forms (`∞−∞`, `0·∞`, `0^0`) fall through to numeric probing.
- Numeric probing classifies behaviour by sampling; exotic functions may be misclassified near pathological points.

### Polynomial Expansion

The `poly` module distributes products and non-negative integer powers into a collected sum of monomials. Multivariate expressions are supported; like terms are collected; output is ordered by descending degree. Non-polynomial parts (function calls, symbolic powers, variable denominators) are left intact while their polynomial children still distribute.

**REPL** (`expand <expr>`):
```
mathr> expand (x+1)^3
x^3 + 3*x^2 + 3*x + 1
mathr> expand (x+y)*(x-y)
x^2 - y^2
mathr> expand (x+2)*(x+3)
x^2 + 5*x + 6
```

**Limitations**: integer exponents up to 64 are expanded; larger (or symbolic) exponents stay as `Pow`. Division only distributes when the denominator is a non-zero constant. Term count is capped at 20,000.

### Partial Fraction Decomposition

The `apart` module decomposes a rational function `N(x)/D(x)` into a polynomial quotient plus a sum of fractions over the linear and irreducible quadratic factors of `D`:

**REPL** (`apart <expr> [<var>]` — the variable may be omitted when the expression has exactly one):
```
mathr> apart 1/(x*(x+1))
-(1/(x + 1)) + 1/x
mathr> apart (x^2+1)/(x-1)
x + 1 + 2/(x - 1)
mathr> apart 1/(x^3+x^2)
1/(x + 1) - 1/x + 1/x^2
```

**Method**: polynomial long division first (`N = Q·D + R`), then numeric factorization of `D` (Durand–Kerner complex roots; conjugate pairs become irreducible quadratics; repeated roots are clustered for multiplicity), then a square linear system for the unknown coefficients solved by Gaussian elimination. Coefficients near integers are snapped for display.

**Limitations**:
- Univariate with a polynomial denominator only (functions in the denominator, or multiple variables, are rejected).
- Denominator degree is capped at 32; repeated roots carry the numerical fuzz inherent to multiple-root finding (coefficients may be off in the last digits).
- Factorization is numeric, so exact rational coefficients are approximated (typically ~12 significant digits).

### Curve Fitting (Levenberg–Marquardt)

Fits `y = f(x, params)` to data by minimizing squared residuals with the LM damped Gauss–Newton algorithm (numeric central-difference Jacobian, Marquardt diagonal scaling).

```
mathr> fit a*x + b with a=1, b=1 0 1 1 3 2 5
a = 2 ± 0
b = 1 ± 0
sse = 0
iterations = 4 (converged)
mathr> fit a*exp(-b*x) in x with a=3, b=0.5 0 3 1 1.4957 2 0.7458 3 0.3717 4 0.1852
a = 3.0000170211 ± 0.0001123634
b = 0.6960330223 ± 0.000051804
sse = 0.0000000401
iterations = 5 (converged)
```

Syntax: `fit <expr> [in <var>] [with p=v, ...] x1 y1 x2 y2 ...`

- Data is the trailing run of `x y` pairs (independent variable defaults to `x`; `in <var>` selects another).
- Parameters come from `with p=v, ...` initial guesses; if omitted, they are inferred as every free variable except the data variable (initial value 1).
- Output: fitted parameters with approximate 1-sigma standard errors (`σ²·(JᵀJ)⁻¹`, only when there are residual degrees of freedom), SSE, iteration count, convergence status.

**Limitations**: local minima (e.g. frequency fitting `sin(w·x)` has a spurious optimum near `w = 0` — start near the expected value); models are evaluated over the expression evaluator, so convergence tolerances are ~1e-12.

### Logistic Regression

Binary classification via IRLS (Newton on the log-likelihood), with Wald standard errors and overflow-safe log-likelihood. The `logit` REPL command takes the 0/1 response followed by `with` and one or more predictor columns separated by `|`:

```
mathr> logit 0 0 1 0 1 1 with 1 1 1 2 2 2
intercept = -2.0794415395 ± 2.7386127872
b1 = 1.3862943598 ± 1.7320508074
log-lik = -3.8190850098
iterations = 5 (converged)
p = 0.3333333333, 0.3333333333, 0.3333333334, 0.6666666666, 0.6666666666, 0.6666666667
```

Balanced grouped designs have a closed form (empirical logits): the example gives `b1 = 2 ln 2`, `b0 = −3 ln 2`. Complete separation reports huge coefficients/standard errors (the MLE does not exist); a `1e-10` ridge keeps the solves well-defined.

**Limitations**: no regularization options (beyond the numerical ridge), no factor predictors (encode numerically first), normal-approximation Wald intervals.

### Stiff ODE Solvers

Library-only (`stiff` module; the explicit solvers in `ode` are likewise library-only): backward (implicit) Euler — first order, L-stable — implicit trapezoidal / Crank–Nicolson — second order, A-stable — and BDF2 — second order, A-stable, fixed step with one trapezoidal startup step — for systems `dy/dt = f(t, y)`. Each step solves the implicit equation by Newton's method with a central-difference numeric Jacobian; linear solves go through `Matrix::solve`.

- `implicit_euler_system(f, t0, t1, y0, n) -> Vec<f64>` — state at `t1`
- `implicit_euler_trajectory(...) -> Vec<(f64, Vec<f64>)>` — full trajectory
- `trapezoid_system` / `trapezoid_trajectory` — same shapes for Crank–Nicolson
- `bdf2_system` / `bdf2_trajectory` — same shapes for fixed-step BDF2

All are unconditionally stable on stiff decays (e.g. `y' = −1000(y−1)` with `h = 0.01`, where explicit Euler amplifies by 9× per step and blows up). Validated by mass conservation through the Robertson problem, exact amplification factors, and observed convergence orders (2.01 / 4.00 / 3.96 error ratios for Euler / trapezoidal / BDF2).

### Optimization

Derivative-free minimization (`optim` module):

- `golden_section(f, a, b, tol) -> (x, fx)` — 1-D unimodal minimization on a bracket; the bracket shrinks by φ = (√5−1)/2 per iteration.
- `nelder_mead(f, start, opts) -> OptResult` — N-dimensional simplex method (reflection/expansion/contraction/shrink), returning `x`, `fx`, `iterations`, `converged`. Tuning via `OptOptions { max_iter, tol, init_step }`.

The `minimize <expr> [var] a b` REPL command runs golden-section search on an expression over `[a, b]` (the single-letter variable may be given just before the bracket; omit both it and brackets ending in a bare variable to default to `x`):

```
mathr> minimize x^2 - 3*x + 2 x 0 5
x* = 1.5000000137
f(1.5000000137) = -0.25
mathr> minimize sin(x) x 3 5
x* = 4.7123889909
f(4.7123889909) = -1
```

**Limitations**: golden-section assumes unimodality on the bracket; Nelder–Mead finds local minima (validated on Rosenbrock to 1e-4); near-flat minima are noise-limited to ~1e-8 accuracy.

### Sparse Matrices & Conjugate Gradient

Sparse storage for large systems (`sparse` module):

- `Csr::from_triplets(rows, cols, &[(r, c, v)...])` — CSR from coordinate triplets; duplicates summed, rows in canonical ascending-column order.
- `Csr::from_dense(&Matrix)` / `to_dense()` — dense conversion (exact zeros dropped).
- `matvec(&[f64])` — matrix–vector product; `transpose() -> Csc` — CSC layout (same matrix, column-compressed).
- `Csc::to_csr()`, `Csc::matvec`, `Csc::to_dense` — CSC counterparts (a CSC's buffers reinterpreted as CSR give the transpose matrix).
- `multiply(&Csr)` — sparse × sparse product (row-wise SpGEMM, canonical output).
- `conjugate_gradient(&Csr, b, tol, max_iter) -> CgResult { x, iterations, residual }` — iterative solve of `A·x = b` for symmetric positive-definite `A`. Symmetry is verified structurally (1e-10 relative); indefiniteness surfaces as a non-positive-curvature breakdown. Converged when `‖r‖₂ ≤ tol · max(1, ‖b‖₂`.
- `conjugate_gradient_jacobi(&Csr, b, tol, max_iter)` — Jacobi (diagonal) preconditioning; requires a nonzero diagonal. Provably performs the same iteration count as plain CG under diagonal scaling `(S·A·S, S·b) ↔ (A, b)`.
- `bicgstab(&Csr, b, tol, max_iter) -> BicgstabResult { x, iterations, residual }` — Jacobi-preconditioned BiCGStab for general (including nonsymmetric) square matrices; no symmetry requirement. Explicit breakdown guards (`rho = 0`, `rhat·v = 0`, `A·shat = 0`, `omega = 0`) return `NotConvergent`.
- `Ilu0::factorize(&Csr) -> Ilu0` — incomplete LU with zero fill-in over A's sparsity pattern (row-wise IKJ); errors on missing diagonals and zero pivots. `ilu.solve(&r)` applies `M⁻¹` via forward/back substitution. Exact for no-fill matrices (tridiagonal, diagonal).
- `bicgstab_ilu(&Csr, b, tol, max_iter, &Ilu0)` — ILU(0)-preconditioned BiCGStab; converges in 1 iteration when the factorization is exact.

```
mathr> spcg 2 1 | 1 3 | 3 4
x = [1, 1]  (2 iterations, residual 0.00e0)
mathr> spcg jacobi 2 1 | 1 3 | 3 4
x = [1, 1]  (2 iterations, residual 5.53e-16)
mathr> spbicg 2 1 | 0 3 | 3 4
x = [0.8333333333, 1.3333333333]  (2 iterations, residual 2.78e-17)
mathr> spbicg ilu 2 1 | 0 3 | 3 4
x = [0.8333333333, 1.3333333333]  (1 iterations, residual 1.07e-14)
```

**Limitations**: CG requires symmetric positive-definite matrices (nonsymmetric input is rejected; indefinite input breaks down); BiCGStab has no monotone residual guarantee and can stagnate on strongly nonsymmetric problems; ILU(0) is fill-limited (strongly coupled grids may need drop-tolerance ILU) and breaks down on zero pivots; the REPL accepts dense row input and sparsifies internally (triplets are the library-level entry point for genuinely large systems).

### Monotone Interpolation (PCHIP)

`pchip x1 y1 x2 y2 ... x_at` evaluates the Fritsch–Carlson monotone piecewise cubic Hermite interpolant at `x_at`. Unlike natural cubic splines it never overshoots the local knot range on monotone data, has zero slope at local extrema, and reproduces linear data exactly. Evaluation outside the knots clamps to the endpoint values.

```
mathr> pchip 0 0 1 0 2 0.3 3 4 4 4.2 2.5
pchip(2.5) = 2.1719391026
```

Accuracy is ~O(h³) on smooth data (the slopes approximate f′). The library type is `pchip::Pchip` (`new`, `eval`, `slopes`).

### B-Spline & Hermite Interpolation

More spline types (`bspline` module + `interpolate::CubicHermite`):

- `bspline x1 y1 x2 y2 ... x_at` — clamped B-spline interpolant evaluated at `x_at`: cubic for 4+ points (de Boor knot averaging guarantees a nonsingular banded collocation system), degree reduced for 2-3 points (linear/quadratic). Passes through every data point; clamps outside `[x_0, x_{n-1}]`.
- `hermite x1 y1 d1 x2 y2 d2 ... x_at` — piecewise cubic Hermite with user-supplied slopes at every knot; knot- and slope-exact, clamps outside the knots.
- Library: `basis_function(i, p, knots, t)` (Cox–de Boor, left-limit at the final knot), `BSpline::new(degree, knots, coeffs)` / `eval` (de Boor, affine parameter map), `cubic_bspline_interp(xs, ys)`, `interpolate::CubicHermite::new(xs, ys, ds)` / `eval` / `derivative`.

```
mathr> bspline 0 0 1 1 2 4 3 9 1.5
bspline(1.5) = 2.25 (degree 3)
mathr> hermite 0 0 0 1 1 2 0.5
hermite(0.5) = 0.25
```

**Limitations**: interpolation only (no least-squares/approximation B-splines, no 2-D surfaces); uniform-in-x parameterization (no chord-length option); a 4-point clamped cubic reproduces any cubic exactly, but general accuracy is ~O(h⁴); Hermite requires slopes as input (use `pchip` for estimated monotone slopes).

### Distributions & Hypothesis Testing

Continuous distributions with PDF and CDF: normal (also via `stats`), exponential, uniform, Student's t, chi-squared, F. Discrete distributions with PMF and CDF: binomial, Poisson.

```
mathr> dist t 2.228 10
pdf = 0.0423946247
cdf = 0.9749941141
mathr> dist binom 5 10 0.5
pmf = 0.24609375
cdf = 0.623046875
```

The library also exposes `normal_ppf` (inverse standard normal CDF), `student_t_ppf`, `chi2_ppf`, `f_ppf`, and `uniform_ppf` quantile functions (critical values), and `beta_inc` (regularized incomplete beta function). The `qtile <dist> <p> <params...>` REPL command gives critical values directly:

```
mathr> qtile t 0.975 10
q = 2.228138852
mathr> qtile normal 0.975 10 2
q = 13.91992797
```

Hypothesis tests return the statistic, degrees of freedom, and p-value:

- `ttest <mu> <values...>` — two-sided one-sample t-test
- `ttest <values...> | <values...>` — two-sided Welch (unequal-variance) t-test
- `ttest paired <values...> | <values...>` — two-sided paired t-test on differences
- `chitest <observed...>` — chi-squared test against uniform cell probabilities
- `chitest <observed...> | <expected...>` — chi-squared goodness-of-fit
- `anova <g1...> | <g2...> [| ...]` — one-way ANOVA F-test
- `mwu <g1...> | <g2...>` — Mann–Whitney U rank-sum test (normal approximation, continuity + tie corrected; reliable for group sizes >= 8)
- `wilcoxon <before...> | <after...>` — Wilcoxon signed-rank paired test (zeros dropped, normal approximation; needs >= 2 non-zero differences)
- `kw <g1...> | <g2...> [| ...]` — Kruskal–Wallis H test by ranks (chi-squared with `k−1` df; reliable for group sizes >= 5)

Also: `spearman_corr(x, y)` (library; Pearson correlation of average ranks) and
`boot <mean|median> <values...> [seed <n>]` — percentile bootstrap 95% confidence
interval from 10000 seeded resamples:

```
mathr> boot median 1 2 3 4 5 6 7 8 9
median = 5
95% ci = [2, 8]
iters = 10000 (seed 42)
```

```
mathr> ttest 0 1 2 3 4 5
one-sample t-test: stat=4.2426406871
df=4
p=0.0132355996
mathr> anova 1 2 3 | 4 5 6 | 7 8 9
one-way ANOVA: stat=27
df=(2, 6)
p=0.001
```

**Method**: chi-squared and Poisson CDFs reuse the regularized incomplete gamma; Student-t, F, and binomial CDFs use the regularized incomplete beta (continued fraction). ANOVA uses the classical between/within sum-of-squares decomposition.

### Numbers

- Integers: `42`
- Decimals: `3.14`
- Scientific: `1.5e3`, `2E-2`

## Constants

| Name | Value |
|------|-------|
| `pi`, `PI` | π ≈ 3.14159... |
| `e` | e ≈ 2.71828... |
| `tau` | τ = 2π |
| `inf`, `Inf`, `Infinity` | +∞ |
| `nan`, `NaN` | NaN |

## Built-in Functions

| Category | Functions |
|----------|-----------|
| Trig | `sin`, `cos`, `tan`, `asin`, `acos`, `atan` |
| Hyperbolic | `sinh`, `cosh`, `tanh` |
| Exp/Log | `exp`, `ln`, `log(x,b)`, `log2`, `log10` |
| Roots | `sqrt`, `cbrt` |
| Rounding | `floor`, `ceil`, `round`, `fract` |
| Other | `abs`, `sign`, `min(...)`, `max(...)`, `pow(x,y)`, `mod(x,y)` |
| Special | `gamma`, `erf`, `erfc`, `sinc`, `bessel_j0`, `bessel_j1`, `bessel_j(n,x)`, `digamma`, `trigamma`, `polygamma(m,x)`, `harmonic(n)`, `zeta(s)`, `hurwitz(s,a)`, `elliptic_k(k)`, `elliptic_e(k)`, `elliptic_f(phi,k)`, `elliptic_e_inc(phi,k)` |

## REPL Commands

| Command | Description |
|---------|-------------|
| `<expr>` | Evaluate |
| `let x = <expr>` | Bind variable |
| `fn f(x) = <expr>` | Define function |
| `diff <expr> [var]` | Symbolic derivative |
| `pdiff <expr> <var>` | Partial derivative |
| `gradient <expr>` | Gradient (all partials) |
| `integrate <expr> [var]` | Symbolic integration |
| `simplify <expr>` | Simplify |
| `int <expr> a b` | Numerical integral |
| `romberg <expr> a b` | Romberg-integrated |
| `solve <expr> [var] [guess]` | Root finding |
| `plot <expr> a b [out.png]` | PNG plot |
| `taylor <expr> [a] [order]` | Taylor series |
| `laurent <expr> [a] [k] [N]` | Laurent series around a pole |
| `rat <a> <op> <b>` | Exact rational arithmetic |
| `cval <expr> [with <var>=<val>,...]` | Complex evaluation (`i` = imaginary unit, principal branch) |
| `qsolve <expr> [= rhs]` | Symbolic linear/quadratic solve (one variable, roots exact when the discriminant is a perfect square) |
| `sum <expr> [<var>] <a> <b>` | Σ over inclusive integer bounds (exact rational results when possible; closed forms for arithmetic/geometric/Faulhaber patterns) |
| `prod <expr> [<var>] <a> <b>` | ∏ over inclusive integer bounds |
| `fourier <expr> L N [x]` | Fourier series on [-L, L] |
| `mc <expr> a b N [seed]` | Monte Carlo integral |
| `sample <dist> <params...> N [seed]` | Random sampling |
| `dist <dist> <x> <params...>` | PDF and CDF (normal/exponential/uniform/t/chi2/f/binom/poisson) |
| `ttest <mu> <values...>` | One-sample t-test |
| `ttest <values...> \| <values...>` | Two-sample Welch t-test |
| `ttest paired <values...> \| <values...>` | Paired t-test |
| `chitest <observed...> [\| <expected...>]` | Chi-squared test |
| `anova <group1...> \| <group2...> [\| ...]` | One-way ANOVA |
| `fft <numbers...>` | Magnitude spectrum |
| `conv <a...> x <b...>` | Convolution |
| `stats <numbers...>` | Descriptive statistics |
| `poly-roots <coeffs...>` | Polynomial roots |
| `isolate-roots <ints...>` | Real root isolation (VAS) |
| `lu <rows...>` | LU decomposition (rows separated by `\|`) |
| `qr <rows...>` | QR decomposition (Householder reflections; prints Q and R) |
| `tikhonov <rows...> \| <b...> <lambda>` | Tikhonov-regularised solve |
| `spcg [jacobi] <rows...> \| <b...>` | Sparse conjugate-gradient solve (SPD; `jacobi` = diagonal preconditioning; reports iterations and residual) |
| `spbicg [ilu] <rows...> \| <b...>` | Sparse BiCGStab solve (nonsymmetric OK; `ilu` = ILU(0) preconditioning, default Jacobi) |
| `cholesky <rows...>` | Cholesky decomposition |
| `svd <rows...>` | Singular value decomposition |
| `eig <rows...>` | Dominant eigenpair (power iteration) |
| `symlig <rows...>` | Full symmetric eigenvalue decomposition (QR algorithm) |
| `hessenberg <rows...>` | Hessenberg decomposition `A = Q·H·Qᵀ` |
| `schur <rows...>` | Real Schur decomposition `A = Q·T·Qᵀ` |
| `rank <rows...>` | Matrix rank |
| `cond <rows...>` | 2-norm condition number σ_max/σ_min (SVD); `inf` for singular |
| `null <rows...>` | Orthonormal nullspace basis of `{x : A·x = 0}` |
| `det <rows...>` | Matrix determinant |
| `spline x1 y1 x2 y2 ... x_at` | Cubic spline at `x_at` |
| `chebyshev n [x]` | Chebyshev `T_n(x)` (or `n` nodes) |
| `legendre n [x]` | Legendre `P_n(x)` (or `n`-point Gauss–Legendre) |
| `gcd / lcm / is-prime / factor / fib / binom / fact / mr-prime` | Number theory |
| `jacobi <a> <n>` | Jacobi symbol |
| `cf <p> <q>` | Continued fraction |
| `diophantine <a> <b> <c>` | Linear Diophantine solver |
| `dlog <g> <h> <p>` | Discrete logarithm |
| `mathml <expr>` / `mathml import <ml>` | Presentation MathML export/import |
| `serialize <fmt> <expr>` / `serialize <fmt> import <t>` | Expression serialization (`fmt`: `sexpr`/`json`/`rpn`) |
| `interval <expr> with <var>=[lo,hi],...` | Rigorous bounds via interval arithmetic |
| `vars` / `funcs` | Show bindings |
| `clear` | Reset context |
| `help` | Help text |
| `quit` | Exit |

## Notebook File Format (`.mnb`)

A math notebook is a JSON file with a `cells` array. Each cell has an `id`, `input` (math expression or TeX), `output` (evaluation result), and `cell_type` (`"math"` or `"text"`):

```json
{
  "cells": [
    { "id": 0, "input": "let x = 5", "output": "x = 5", "cell_type": "math" },
    { "id": 1, "input": "x * 3", "output": "= 15", "cell_type": "math" },
    { "id": 2, "input": "# Quadratic formula", "output": "# Quadratic formula", "cell_type": "text" }
  ]
}
```

The `cell_type` field is optional in old `.mnb` files (defaults to `"math"`).

Start the web UI with `mathr notebook [file.mnb] [port]` (default port 3000).

### Web UI Features

- **Shared context** — variables and functions defined via `let`/`fn` in one cell persist across subsequent cells (like Jupyter kernels)
- **Cell types** — Math cells (evaluated, KaTeX-rendered) and Text cells (Markdown-rendered documentation)
- **Inline plots** — `plot` commands render PNG images directly in the notebook via base64-encoded responses
- **Cell management** — add, delete, duplicate, move up/down, toggle type
- **Execution status** — each cell shows running/done/error status with `In [n]:` execution counters
- **Context panel** — collapsible panel showing bound variables and user functions
- **Reset & Run All** — resets the shared context and re-evaluates all cells in order
- **Markdown rendering** — text cells render Markdown (headings, lists, code, blockquotes) via marked.js
- **KaTeX rendering** — input expressions and output results are rendered as math notation
- **Step-by-step solving** — `POST /api/eval` returns a `steps` array with intermediate steps for `diff`, `solve`, `taylor`, `integrate`, `simplify`, `rat`, `laurent`
- **Exact fraction arithmetic** — expressions with integer fractions (e.g. `\frac{1}{2} + \frac{3}{4}`) are evaluated exactly as `Rational`, returning `5/4` instead of `1.25`
- **Live input preview** — each cell shows a rendered math/Markdown preview as you type
- **Keyboard shortcuts** — Shift/Cmd/Ctrl+Enter to run a cell, Alt+Enter to run and add a new cell

### API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| `GET` | `/` | Serve web UI HTML |
| `POST` | `/api/eval` | Evaluate expression (updates shared context); returns `{input, output, steps, image?}` where `image` is base64 PNG for `plot` commands |
| `GET` | `/api/notebook` | Get current notebook as JSON |
| `POST` | `/api/notebook` | Replace notebook state (auto-saves to file) |
| `POST` | `/api/save` | Save notebook to file |
| `POST` | `/api/reset` | Reset the shared evaluation context |
| `GET` | `/api/context` | Get current variables and user functions as `{vars, funcs}` |

## Error Handling

All library functions return `Result<T, MathError>` with variants:
- `Parse` — syntax errors
- `Eval` — evaluation errors (wrong arg count, etc.)
- `UnknownVariable` / `UnknownFunction`
- `Domain` — domain errors (sqrt of negative, log of 0)
- `NotConvergent` — solver/ODE failed to converge
- `InvalidArgument` — bad input dimensions/values
- `Io` — I/O errors
- `Plot` — rendering errors
- `Other` — catch-all