//! Command-line entry point.
//!
//! Usage:
//!   mathr "sin(pi/4) + 2^3"       evaluate an expression
//!   mathr "diff x^3"              symbolic derivative
//!   mathr "solve x^2 - 4"         find roots
//!   mathr "int sin(x) 0 pi"       numerical integral
//!   mathr "gcd 48 36"             number theory
//!   mathr "is-prime 97"           primality test
//!   mathr "gamma(0.5)"            special functions (as expressions)
//!   echo "sin(pi/4)" | mathr      read from stdin
//!   mathr                         interactive REPL

use anyhow::Result;
use clap::Parser as ClapParser;
use mathr::eval::Context;

#[derive(ClapParser)]
#[command(
    name = "mathr",
    version,
    about = "Mathematical computation in your terminal — just give it an expression.",
    long_about = None
)]
struct Cli {
    /// Expression or command string (e.g. "sin(pi/4)", "diff x^3", "solve x^2-4").
    /// If omitted, reads from stdin (when piped) or starts the REPL.
    input: Option<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let input = match cli.input {
        Some(s) => s,
        None => {
            use std::io::IsTerminal;
            if std::io::stdin().is_terminal() {
                return mathr::repl::run().map_err(Into::into);
            }
            use std::io::Read;
            let mut buf = String::new();
            std::io::stdin().read_to_string(&mut buf)?;
            buf.trim().to_string()
        }
    };

    if input.is_empty() {
        return mathr::repl::run().map_err(Into::into);
    }

    let ctx = Context::standard();
    match mathr::repl::dispatch_str(&input, ctx) {
        Ok(Some(s)) => { println!("{}", s); Ok(()) }
        Ok(None) => Ok(()),
        Err(e) => Err(anyhow::anyhow!("{}", e)),
    }
}