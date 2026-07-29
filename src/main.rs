//! Command-line entry point. Each subcommand maps to a feature exposed by
//! the library so the same code drives both REPL and scripted use.

use anyhow::{Context as AnyContext, Result};
use clap::{Parser as ClapParser, Subcommand};
use maths::complex::Complex;
use maths::error::MathError;
use maths::eval::{eval, Context};
use maths::expr::Expr;
use maths::parser::Parser;
use maths::solver::SolveOptions;

#[derive(ClapParser)]
#[command(
    name = "maths",
    version,
    about = "Mathematical computation in your terminal: FFT, calculus, algebra and plotting.",
    long_about = None
)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Evaluate an expression.
    Eval {
        /// Expression to evaluate, e.g. `sin(pi/4) + 2^3`.
        expr: String,
        /// Bind variables as `name=value`, e.g. `x=3`.
        #[arg(long, value_name = "NAME=VALUE", num_args = 0..)]
        set: Vec<String>,
    },
    /// Symbolic derivative with respect to `var` (default `x`).
    Diff {
        expr: String,
        #[arg(long, default_value = "x")]
        var: String,
        #[arg(long, default_value_t = false)]
        simplify: bool,
    },
    /// Simplify an expression (constant folding, identities).
    Simplify { expr: String },
    /// Numerical integral over `[a, b]`.
    Integrate {
        expr: String,
        #[arg(long, default_value = "x")]
        var: String,
        a: f64,
        b: f64,
        /// Number of panels (Simpson's rule). Use `--adaptive` instead.
        #[arg(long)]
        n: Option<usize>,
        #[arg(long, default_value_t = false)]
        adaptive: bool,
    },
    /// Solve `expr = 0` for `var` starting from `guess`.
    Solve {
        expr: String,
        #[arg(long, default_value = "x")]
        var: String,
        #[arg(long, default_value_t = 0.0)]
        guess: f64,
        /// Bisection search on `[a, b]` instead of Newton.
        #[arg(long, num_args = 2, value_names = ["A", "B"])]
        bisect: Option<Vec<f64>>,
        #[arg(long, default_value_t = 100)]
        max_iter: usize,
        #[arg(long, default_value_t = 1e-10)]
        tol: f64,
    },
    /// Find roots of a polynomial: coefficients highest-degree first.
    PolyRoots { coefficients: Vec<f64> },
    /// Plot `f(x)` to a PNG file.
    Plot {
        expr: String,
        /// Output PNG path (default `plot.png`).
        #[arg(short, long, default_value = "plot.png")]
        output: String,
        #[arg(long, default_value = "x")]
        var: String,
        #[arg(long, default_value_t = -6.283185307179586)]
        a: f64,
        #[arg(long, default_value_t = 6.283185307179586)]
        b: f64,
        #[arg(long, default_value_t = 800)]
        samples: usize,
        #[arg(long, default_value_t = 60)]
        height: u32,
        #[arg(long, default_value_t = 40)]
        width: u32,
    },
    /// Compute the FFT of comma- or whitespace-separated samples.
    Fft {
        /// Real samples. If `--complex`, the data is interleaved re/im.
        samples: Vec<String>,
        #[arg(long, default_value_t = false)]
        complex: bool,
        #[arg(long, default_value_t = false)]
        inverse: bool,
        #[arg(long, default_value_t = false)]
        magnitude: bool,
        #[arg(long, default_value_t = false)]
        power: bool,
    },
    /// Start the interactive REPL.
    Repl,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::Eval { expr, set } => cmd_eval(&expr, &set),
        Cmd::Diff { expr, var, simplify } => cmd_diff(&expr, &var, simplify),
        Cmd::Simplify { expr } => cmd_simplify(&expr),
        Cmd::Integrate { expr, var, a, b, n, adaptive } => cmd_integrate(&expr, &var, a, b, n, adaptive),
        Cmd::Solve { expr, var, guess, bisect, max_iter, tol } => {
            cmd_solve(&expr, &var, guess, bisect, max_iter, tol)
        }
        Cmd::PolyRoots { coefficients } => cmd_poly_roots(&coefficients),
        Cmd::Plot { expr, output, var, a, b, samples, height, width } => {
            cmd_plot(&expr, &output, &var, a, b, samples, height, width)
        }
        Cmd::Fft { samples, complex, inverse, magnitude, power } => cmd_fft(&samples, complex, inverse, magnitude, power),
        Cmd::Repl => maths::repl::run().map_err(Into::into),
    }
}

fn cmd_eval(expr: &str, set: &[String]) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    let mut ctx = Context::standard();
    for s in set {
        let (k, v) = s.split_once('=').context("`--set` requires name=value")?;
        let v: f64 = v.parse().context("could not parse value as number")?;
        ctx.set(k, v);
    }
    let v = eval(&e, &ctx).map_err(map_err)?;
    println!("{}", v);
    Ok(())
}

fn cmd_diff(expr: &str, var: &str, simplify: bool) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    let d = maths::symbolic::differentiate(&e, var).map_err(map_err)?;
    if simplify {
        println!("{}", maths::simplify::simplify(&d));
    } else {
        println!("{}", d);
    }
    Ok(())
}

fn cmd_simplify(expr: &str) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    println!("{}", maths::simplify::simplify(&e));
    Ok(())
}

fn cmd_integrate(expr: &str, var: &str, a: f64, b: f64, n: Option<usize>, adaptive: bool) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    let ctx = Context::standard();
    let expr_clone = e.clone();
    let ctx_clone = ctx.clone();
    let var_owned = var.to_string();
    let f = move |x: f64| {
        let mut c = ctx_clone.clone();
        c.set(&var_owned, x);
        eval(&expr_clone, &c).unwrap_or(f64::NAN)
    };
    let v = if adaptive || n.is_none() {
        maths::calculus::integrate_adaptive(f, a, b, 1e-9, 30).map_err(map_err)?
    } else {
        let n = n.unwrap();
        maths::calculus::integrate_simpson(f, a, b, n).map_err(map_err)?
    };
    println!("{}", v);
    Ok(())
}

fn cmd_solve(
    expr: &str,
    var: &str,
    guess: f64,
    bisect: Option<Vec<f64>>,
    max_iter: usize,
    tol: f64,
) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    let opts = SolveOptions { max_iter, tol, h: 1e-6 };
    let ctx = Context::standard();
    let var_owned = var.to_string();
    let f = move |x: f64| {
        let mut c = ctx.clone();
        c.set(&var_owned, x);
        eval(&e, &c).unwrap_or(f64::NAN)
    };
    match bisect {
        Some(bracket) if bracket.len() == 2 => {
            let (r, fv) = maths::solver::bisect(f, bracket[0], bracket[1], opts).map_err(map_err)?;
            println!("root = {}  (f = {})", r, fv);
        }
        _ => {
            let (r, fv) = maths::solver::newton_central(f, guess, opts).map_err(map_err)?;
            println!("root = {}  (f = {})", r, fv);
        }
    }
    Ok(())
}

fn cmd_poly_roots(coeffs: &[f64]) -> Result<()> {
    let r = maths::solver::polynomial_roots(coeffs).map_err(map_err)?;
    if r.is_empty() {
        println!("(no real roots found)");
    } else {
        for (x, fx) in r {
            println!("x = {}  (f = {})", x, fx);
        }
    }
    Ok(())
}

fn cmd_plot(
    expr: &str,
    output: &str,
    var: &str,
    a: f64,
    b: f64,
    samples: usize,
    _height: u32,
    _width: u32,
) -> Result<()> {
    let e = Parser::parse(expr).map_err(map_err)?;
    maths::plot::plot_function(output, &e, var, a, b, samples, &format!("y = {}", expr))
        .map_err(map_err)?;
    println!("wrote {}", output);
    Ok(())
}

fn cmd_fft(samples: &[String], complex: bool, inverse: bool, magnitude: bool, power: bool) -> Result<()> {
    if complex {
        let nums: Result<Vec<f64>> = samples
            .iter()
            .map(|s| s.parse::<f64>().context("could not parse sample"))
            .collect();
        let nums = nums?;
        if nums.len() % 2 != 0 {
            return Err(anyhow::anyhow!("complex samples need an even count (re/im interleaved)"));
        }
        let mut data: Vec<Complex<f64>> = Vec::with_capacity(nums.len() / 2);
        for chunk in nums.chunks(2) {
            data.push(Complex::new(chunk[0], chunk[1]));
        }
        let out = if inverse {
            maths::fft::ifft(&data).map_err(map_err)?
        } else {
            maths::fft::fft(&data).map_err(map_err)?
        };
        for v in out {
            println!("{} {} {}", v.re, v.im, v.abs());
        }
    } else {
        let samples: Result<Vec<f64>> = samples
            .iter()
            .map(|s| s.parse::<f64>().context("could not parse sample"))
            .collect();
        let samples = samples?;
        // pad to a power of two if needed
        let n = maths::fft::next_pow2(samples.len());
        let mut padded = samples;
        padded.resize(n, 0.0);
        if inverse {
            // rebuild complex from real, run ifft
            let data: Vec<Complex<f64>> = padded.iter().map(|&x| Complex::new(x, 0.0)).collect();
            let out = maths::fft::ifft(&data).map_err(map_err)?;
            for v in out {
                println!("{}", v.re);
            }
        } else if magnitude {
            let m = maths::fft::magnitude_spectrum(&padded).map_err(map_err)?;
            for v in m {
                println!("{}", v);
            }
        } else if power {
            let p = maths::fft::power_spectrum(&padded).map_err(map_err)?;
            for v in p {
                println!("{}", v);
            }
        } else {
            let full = maths::fft::fft(
                &padded.iter().map(|&x| Complex::new(x, 0.0)).collect::<Vec<_>>(),
            )
            .map_err(map_err)?;
            for v in full {
                println!("{} {} {}", v.re, v.im, v.abs());
            }
        }
    }
    Ok(())
}

fn map_err(e: MathError) -> anyhow::Error {
    anyhow::anyhow!("{}", e)
}

// suppress unused warnings when the binary embeds but doesn't use Expr
#[allow(dead_code)]
fn _expr_used(_: &Expr) {}