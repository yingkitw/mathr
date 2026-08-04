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
//!   mathr notebook [file.mnb]     web notebook UI
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
    #[arg(trailing_var_arg = true, allow_hyphen_values = true)]
    input: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let input = if cli.input.is_empty() {
        use std::io::IsTerminal;
        if std::io::stdin().is_terminal() {
            return mathr::repl::run().map_err(Into::into);
        }
        use std::io::Read;
        let mut buf = String::new();
        std::io::stdin().read_to_string(&mut buf)?;
        buf.trim().to_string()
    } else {
        cli.input.join(" ")
    };

    if input.is_empty() {
        return mathr::repl::run().map_err(Into::into);
    }

    // Check for "notebook" subcommand
    if input == "notebook" || input.starts_with("notebook ") {
        let args: Vec<&str> = input.split_whitespace().collect();
        let file = args.get(1).map(|s| std::path::PathBuf::from(s));
        let port: u16 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3000);
        return run_notebook(file, port).map_err(Into::into);
    }

    let ctx = Context::standard();
    match mathr::repl::dispatch_str(&input, ctx) {
        Ok(Some(s)) => { println!("{}", s); Ok(()) }
        Ok(None) => Ok(()),
        Err(e) => Err(anyhow::anyhow!("{}", e)),
    }
}

fn run_notebook(file: Option<std::path::PathBuf>, port: u16) -> mathr::error::Result<()> {
    let notebook = match &file {
        Some(path) => {
            if path.exists() {
                mathr::notebook::Notebook::load(path)?
            } else {
                mathr::notebook::Notebook::new()
            }
        }
        None => {
            let mut nb = mathr::notebook::Notebook::new();
            nb.add_cell("sin(pi/4)");
            nb.add_cell("\\frac{1}{2} + \\frac{3}{4}");
            nb.add_cell("diff x^3");
            nb
        }
    };
    let server = mathr::server::NotebookServer::new(notebook, file);
    server.serve(port)
}