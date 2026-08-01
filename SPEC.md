# Specification

## CLI Subcommands

| Command | Syntax | Description |
|---------|--------|-------------|
| `eval` | `mathr eval <expr> [--set name=value]` | Evaluate expression |
| `diff` | `mathr diff <expr> [--var x] [--simplify]` | Symbolic derivative |
| `simplify` | `mathr simplify <expr>` | Constant-fold & simplify |
| `integrate` | `mathr integrate <expr> [--var x]` | Symbolic indefinite integral |
| `integrate-num` | `mathr integrate-num <expr> a b [--var x] [--n N] [--adaptive] [--romberg LEVELS]` | Numerical integral |
| `solve` | `mathr solve <expr> [--var x] [--guess 0] [--bisect A B] [--max-iter 100] [--tol 1e-10]` | Root finding |
| `solve-system` | `mathr solve-system <sys> [--guess x0,y0,...]` | Newton's method for nonlinear systems |
| `poly-roots` | `mathr poly-roots <coeffs...>` | Polynomial roots (Durand–Kerner) |
| `plot` | `mathr plot <expr> [-o path] [--var x] [--a -τ] [--b τ] [--samples 800]` | PNG plot |
| `fft` | `mathr fft <samples...> [--complex] [--inverse] [--magnitude] [--power]` | FFT |
| `conv` | `mathr conv <a...> <b...>` | FFT convolution of two signals |
| `taylor` | `mathr taylor <expr> [--var x] [--around 0] [--order 5]` | Taylor series |
| `stats` | `mathr stats <data...>` | Descriptive statistics |
| `matrix` | `mathr matrix <op> <rows...>` | `lu`/`cholesky`/`svd`/`eig`/`rank`/`solve` |
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
| `repl` | `mathr repl` | Interactive REPL |

## Expression Grammar

```
expr   := term (('+' | '-') term)*
term   := factor (('*' | '/') factor)*         -- includes implicit multiplication
factor := unary ('^' factor)?                   -- right-associative
unary  := ('+' | '-')? atom
atom   := number | ident | ident '(' args ')' | '(' expr ')'
args   := expr (',' expr)*
```

### Implicit Multiplication

`2x` → `2*x`, `3(x+1)` → `3*(x+1)`, `(x)(y)` → `x*y`

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
| `integrate <expr> [var]` | Symbolic integration |
| `simplify <expr>` | Simplify |
| `int <expr> a b` | Numerical integral |
| `romberg <expr> a b` | Romberg-integrated |
| `solve <expr> [var] [guess]` | Root finding |
| `plot <expr> a b [out.png]` | PNG plot |
| `taylor <expr> [a] [order]` | Taylor series |
| `fft <numbers...>` | Magnitude spectrum |
| `conv <a...> x <b...>` | Convolution |
| `stats <numbers...>` | Descriptive statistics |
| `poly-roots <coeffs...>` | Polynomial roots |
| `lu <rows...>` | LU decomposition (rows separated by `\|`) |
| `cholesky <rows...>` | Cholesky decomposition |
| `svd <rows...>` | Singular value decomposition |
| `eig <rows...>` | Dominant eigenpair |
| `rank <rows...>` | Matrix rank |
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