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
| `dist` | `mathr dist <dist> <x> <params...>` | PDF and CDF (normal/exponential) |
| `stats` | `mathr stats <data...>` | Descriptive statistics |
| `matrix` | `mathr matrix <op> <rows...>` | `lu`/`cholesky`/`svd`/`eig`/`symlig`/`hessenberg`/`schur`/`rank`/`det`/`solve` |
| `tikhonov` | `mathr tikhonov <rows...> \| <b...> <lambda>` | Tikhonov-regularised solve |
| `interp` | `mathr interp <op> ...` | `lagrange`/`newton`/`spline`/`chebyshev`/`legendre` |
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
| Special | `gamma`, `erf`, `erfc`, `sinc`, `bessel_j0`, `bessel_j1`, `bessel_j(n,x)` |

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
| `fourier <expr> L N [x]` | Fourier series on [-L, L] |
| `mc <expr> a b N [seed]` | Monte Carlo integral |
| `sample <dist> <params...> N [seed]` | Random sampling |
| `dist <dist> <x> <params...>` | PDF and CDF |
| `fft <numbers...>` | Magnitude spectrum |
| `conv <a...> x <b...>` | Convolution |
| `stats <numbers...>` | Descriptive statistics |
| `poly-roots <coeffs...>` | Polynomial roots |
| `isolate-roots <ints...>` | Real root isolation (VAS) |
| `lu <rows...>` | LU decomposition (rows separated by `\|`) |
| `tikhonov <rows...> \| <b...> <lambda>` | Tikhonov-regularised solve |
| `cholesky <rows...>` | Cholesky decomposition |
| `svd <rows...>` | Singular value decomposition |
| `eig <rows...>` | Dominant eigenpair (power iteration) |
| `symlig <rows...>` | Full symmetric eigenvalue decomposition (QR algorithm) |
| `hessenberg <rows...>` | Hessenberg decomposition `A = Q·H·Qᵀ` |
| `schur <rows...>` | Real Schur decomposition `A = Q·T·Qᵀ` |
| `rank <rows...>` | Matrix rank |
| `det <rows...>` | Matrix determinant |
| `spline x1 y1 x2 y2 ... x_at` | Cubic spline at `x_at` |
| `chebyshev n [x]` | Chebyshev `T_n(x)` (or `n` nodes) |
| `legendre n [x]` | Legendre `P_n(x)` (or `n`-point Gauss–Legendre) |
| `gcd / lcm / is-prime / factor / fib / binom / fact / mr-prime` | Number theory |
| `jacobi <a> <n>` | Jacobi symbol |
| `cf <p> <q>` | Continued fraction |
| `diophantine <a> <b> <c>` | Linear Diophantine solver |
| `dlog <g> <h> <p>` | Discrete logarithm |
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