//! The MCP tool catalogue and execution.
//!
//! Every tool is a thin wrapper that builds a mathr command string and
//! runs it through [`crate::repl::dispatch_str`], so the MCP surface
//! cannot drift from the CLI.

use super::json::Json;
use crate::error::{MathError, Result};

/// A tool definition: name, description and JSON Schema for its arguments.
pub struct Tool {
    pub name: &'static str,
    pub description: &'static str,
    pub input_schema: fn() -> Json,
}

const EXPR_DESC: &str = "Math expression. Plain syntax (x^2 + 3*x) or LaTeX (\\frac{1}{2}, \\sin(x)).";
const VAR_DESC: &str =
    "Variable to operate on (default: auto-detected, usually x).";

/// The tools exposed by the server.
///
/// The set mirrors the `unsolve` TypeScript port's MCP surface, adapted to
/// mathr's command vocabulary.
pub fn tools() -> Vec<Tool> {
    vec![
        Tool {
            name: "evaluate",
            description: "Evaluate a math expression. Exact rational arithmetic for fractions (1/2 + 3/4 -> 5/4); numeric results otherwise (sin(pi/4) -> 0.7071...). Supports LaTeX input, constants pi/e/tau, and functions sin/cos/tan/exp/ln/log/sqrt/gamma and more.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variables".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("object")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Optional variable substitutions, e.g. {\"x\": 2}",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "differentiate",
            description: "Symbolic differentiation with step-by-step working (sum, product, quotient, chain, power rules). Returns the derivative in plain syntax, e.g. diff(sin(x^2)) -> cos(x^2)*2*x.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "integrate",
            description: "Integration: a symbolic antiderivative (+ C) by default, or a numeric definite integral when both 'from' and 'to' are given. Bounds may be numbers or expressions such as \"pi\" or \"pi/2\".",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Lower bound of a definite integral (needs 'to' as well)",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Upper bound of a definite integral (needs 'from' as well)",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "solve",
            description: "Solve an equation for a variable. Exact roots for linear/quadratic polynomials (rational arithmetic), rational roots + numeric roots for higher degrees, and robust numeric root-finding for non-polynomial equations.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "equation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Equation like 'x^2 - 5*x + 6 = 0' or an expression (assumed = 0)",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("equation")])),
                ])
            },
        },
        Tool {
            name: "nsolve",
            description: "Find ALL real roots of an equation by scanning an interval for sign changes and bisecting each one. Unlike 'solve', which returns a single root near an initial guess, this returns every root: nsolve(sin(x)) over [-10, 10] gives 0, +/-pi, +/-2pi, +/-3pi. Poles are rejected (a sign flip across a pole is not a root).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Lower bound of the scan (default -30)"),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Upper bound of the scan (default 30)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "simplify",
            description: "Simplify an expression with constant folding and algebraic identities (collecting like terms, 0 + x -> x, x/1 -> x, exact fraction arithmetic).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "taylor",
            description: "Symbolic Taylor series expansion of a function around a point with N terms, e.g. exp(x) around 0 -> 1 + x + 0.5*x^2 + ...",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "at".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Expansion point (default 0)"),
                                    ),
                                ]),
                            ),
                            (
                                "terms".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    (
                                        "description".into(),
                                        Json::str("Number of terms (default 5)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "limit",
            description: "Limit of an expression as a variable approaches a point. Handles direct substitution, 0/0 and infinity/infinity via recursive L'Hopital's rule, numeric probing, and limits at +/- infinity. Example: sin(x)/x at 0 -> 1.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "at".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str("Point to approach (default 0)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "poly-roots",
            description: "All roots of a polynomial with numeric coefficients - real AND complex - via the Durand-Kerner iteration. Coefficients in descending order: [1, -5, 6] solves x^2 - 5x + 6 = 0 -> roots 2, 3.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "coefficients".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("number"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Polynomial coefficients, leading (highest degree) first",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("coefficients")]),
                    ),
                ])
            },
        },
        Tool {
            name: "number-theory",
            description: "Integer number theory: gcd, lcm, factor (prime factorization), is-prime (deterministic Miller-Rabin), fib (fast doubling, exact BigInt), binom (exact binomial coefficient).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "operation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "enum".into(),
                                        Json::Array(vec![
                                            Json::str("gcd"),
                                            Json::str("lcm"),
                                            Json::str("factor"),
                                            Json::str("is-prime"),
                                            Json::str("fib"),
                                            Json::str("binom"),
                                        ]),
                                    ),
                                ]),
                            ),
                            (
                                "a".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    ("description".into(), Json::str("First integer")),
                                ]),
                            ),
                            (
                                "b".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("integer")),
                                    (
                                        "description".into(),
                                        Json::str("Second integer (gcd/lcm/binom)"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("operation"), Json::str("a")]),
                    ),
                ])
            },
        },
        Tool {
            name: "stats",
            description: "Descriptive statistics over a numeric dataset: n, sum, mean, median, min, max, range, sample variance, and sample standard deviation.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "data".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("number"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str("Numeric dataset"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("data")])),
                ])
            },
        },
        Tool {
            name: "matrix",
            description: "Matrix operations on a square/rectangular numeric matrix: det, inv, transpose, linsolve (solve A*x = b given the augmented n x (n+1) matrix [A | b]), trace, eigen (eigenvalues + eigenvectors of a real symmetric matrix).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "operation".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "enum".into(),
                                        Json::Array(vec![
                                            Json::str("det"),
                                            Json::str("inv"),
                                            Json::str("transpose"),
                                            Json::str("linsolve"),
                                            Json::str("trace"),
                                            Json::str("eigen"),
                                        ]),
                                    ),
                                ]),
                            ),
                            (
                                "matrix".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("array")),
                                    (
                                        "items".into(),
                                        Json::obj(vec![("type".into(), Json::str("array"))]),
                                    ),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Matrix as rows of numbers, e.g. [[4, 7], [2, 6]]; for linsolve, the augmented matrix [[2, 1, 5], [1, 3, 10]]",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("operation"), Json::str("matrix")]),
                    ),
                ])
            },
        },
        Tool {
            name: "minimize",
            description: "Minimize a single-variable function over a bounded interval via golden-section search (assumes the function is unimodal on the interval). Returns the argmin x and the minimum value f(x).",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                            (
                                "from".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Lower bound of the search interval (a number, or an expression such as \"pi\")",
                                        ),
                                    ),
                                ]),
                            ),
                            (
                                "to".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("number")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "Upper bound of the search interval (a number, or an expression such as \"2*pi\")",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![
                            Json::str("expression"),
                            Json::str("from"),
                            Json::str("to"),
                        ]),
                    ),
                ])
            },
        },
        Tool {
            name: "msolve",
            description: "Solve a system of 2-4 linear equations EXACTLY, using rational arithmetic (the answer stays a fraction, never a decimal). Equations are separated by ';', e.g. '2x + y = 5; x + 3y = 10' -> x = 1, y = 3. The right-hand side may be omitted (it defaults to 0). A singular system or a non-linear equation is reported as an error.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "equations".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "description".into(),
                                        Json::str(
                                            "The equations, separated by ';', e.g. \"2x + y = 5; x + 3y = 10\"",
                                        ),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("equations")])),
                ])
            },
        },
        Tool {
            name: "expand-in",
            description: "Expand an expression in ONE variable, treating every other sub-expression as an opaque exact coefficient (never rounded to a float). For example expand-in '(x + sin(y))^2' with variable 'x' gives x^2 + 2*sin(y)*x + sin(y)*sin(y). This reaches cases the ordinary expand cannot, since it only distributes when both sides are polynomials.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    (
                                        "description".into(),
                                        Json::str("The variable to expand in"),
                                    ),
                                ]),
                            ),
                        ]),
                    ),
                    (
                        "required".into(),
                        Json::Array(vec![Json::str("expression"), Json::str("variable")]),
                    ),
                ])
            },
        },
        Tool {
            name: "apart",
            description: "Partial fraction decomposition of a rational function whose denominator factors into linear factors (e.g. 1/(x^2-1) -> 1/2/(x-1) - 1/2/(x+1)). Handles improper fractions via polynomial long division.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                            (
                                "variable".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(VAR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
        Tool {
            name: "expand",
            description: "Expand a polynomial expression into a collected sum of monomials (multivariate), e.g. (x + y)^2 -> x^2 + 2*x*y + y^2.",
            input_schema: || {
                Json::obj(vec![
                    ("type".into(), Json::str("object")),
                    (
                        "properties".into(),
                        Json::obj(vec![
                            (
                                "expression".into(),
                                Json::obj(vec![
                                    ("type".into(), Json::str("string")),
                                    ("description".into(), Json::str(EXPR_DESC)),
                                ]),
                            ),
                        ]),
                    ),
                    ("required".into(), Json::Array(vec![Json::str("expression")])),
                ])
            },
        },
    ]
}

// =========================================================================
// Tool execution
// =========================================================================

/// Run a tool by name with the given arguments, returning the text result.
///
/// Errors are returned as `Err` and reported by the caller as an MCP tool
/// execution error (`isError: true`).
pub fn call_tool(name: &str, args: &Json) -> Result<String> {
    match name {
        "evaluate" => {
            let expr = require_str(args, "expression")?;
            let mut ctx = crate::eval::Context::standard();
            if let Some(Json::Object(members)) = args.get("variables") {
                for (k, v) in members {
                    let Some(n) = v.as_f64() else {
                        return Err(MathError::Eval(format!(
                            "variable '{}' must be a number",
                            k
                        )));
                    };
                    ctx.set(k, n);
                }
            }
            let e = crate::parser::Parser::parse(&expr)?;
            // Prefer an exact rational result when the expression is rational,
            // matching the REPL's plain-evaluation path.
            if let Some(r) = crate::rational::eval_rational(&e) {
                return Ok(r.to_string());
            }
            let v = crate::eval::eval(&e, &ctx)?;
            Ok(crate::repl::format_value(v))
        }
        "differentiate" => {
            let expr = require_str(args, "expression")?;
            let cmd = format!("diff {}", expr);
            run_with_var(cmd, args)
        }
        "integrate" => {
            let expr = require_str(args, "expression")?;
            // With both bounds, compute the numeric definite integral.
            if args.get("from").is_some() && args.get("to").is_some() {
                let from = bound_text(args, "from")?;
                let to = bound_text(args, "to")?;
                let cmd = format!("int {} from {} to {}", expr, from, to);
                return run_with_var(cmd, args);
            }
            // Otherwise the symbolic antiderivative.
            let cmd = format!("integrate {}", expr);
            run_with_var(cmd, args)
        }
        "solve" => {
            let expr = require_str(args, "equation")?;
            let cmd = format!("solve {}", expr);
            run_with_var(cmd, args)
        }
        "nsolve" => {
            // All real roots by scanning for sign changes.
            let expr = require_str(args, "expression")?;
            let mut cmd = format!("nsolve {}", expr);
            if let Some(from) = args.get("from").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(from)));
            }
            if let Some(to) = args.get("to").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(to)));
            }
            run_with_var(cmd, args)
        }
        "simplify" => {
            let expr = require_str(args, "expression")?;
            run_command(&format!("simplify {}", expr))
        }
        "expand" => {
            let expr = require_str(args, "expression")?;
            run_command(&format!("expand {}", expr))
        }
        "apart" => {
            let expr = require_str(args, "expression")?;
            run_with_var(format!("apart {}", expr), args)
        }
        "msolve" => {
            // Solve a system of linear equations exactly.
            let equations = require_str(args, "equations")?;
            let solution = crate::linsolve::solve_system(&equations)?;
            let parts: Vec<String> = solution
                .iter()
                .map(|(v, r)| format!("{} = {}", v, r.to_string()))
                .collect();
            Ok(parts.join(", "))
        }
        "expand-in" => {
            // Expand in one variable, others become exact coefficients.
            let expr = require_str(args, "expression")?;
            let var = require_str(args, "variable")?;
            let e = crate::parser::Parser::parse(&expr)?;
            Ok(crate::poly::expand_in(&e, &var).to_string())
        }
        "taylor" => {
            let expr = require_str(args, "expression")?;
            let at = args.get("at").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let terms = args
                .get("terms")
                .and_then(|v| v.as_f64())
                .unwrap_or(5.0)
                .clamp(1.0, 64.0) as usize;
            let mut cmd = format!("taylor {}", expr);
            if at != 0.0 {
                cmd.push_str(&format!(" {}", crate::repl::format_value(at)));
            }
            cmd.push_str(&format!(" {}", terms));
            run_with_var(cmd, args)
        }
        "limit" => {
            let expr = require_str(args, "expression")?;
            let at = args.get("at").and_then(|v| v.as_f64()).unwrap_or(0.0);
            let mut cmd = format!("limit {}", expr);
            cmd.push_str(&format!(" {}", crate::repl::format_value(at)));
            run_with_var(cmd, args)
        }
        "poly-roots" => {
            let coeffs = args
                .get("coefficients")
                .and_then(|v| v.as_array())
                .ok_or_else(|| {
                    MathError::Eval("poly-roots needs a 'coefficients' array".into())
                })?;
            let mut parts = Vec::with_capacity(coeffs.len());
            for c in coeffs {
                let n = c.as_f64().ok_or_else(|| {
                    MathError::Eval("poly-roots coefficients must be numbers".into())
                })?;
                parts.push(crate::repl::format_value(n));
            }
            run_command(&format!("poly-roots {}", parts.join(" ")))
        }
        "number-theory" => {
            let op = require_str(args, "operation")?;
            let a = require_num(args, "a")?;
            let mut cmd = format!("{} {}", op, crate::repl::format_value(a));
            if let Some(b) = args.get("b").and_then(|v| v.as_f64()) {
                cmd.push_str(&format!(" {}", crate::repl::format_value(b)));
            }
            run_command(&cmd)
        }
        "stats" => {
            let data = args
                .get("data")
                .and_then(|v| v.as_array())
                .ok_or_else(|| MathError::Eval("stats needs a 'data' array".into()))?;
            let mut parts = Vec::with_capacity(data.len());
            for d in data {
                let n = d
                    .as_f64()
                    .ok_or_else(|| MathError::Eval("stats data must be numbers".into()))?;
                parts.push(crate::repl::format_value(n));
            }
            run_command(&format!("stats {}", parts.join(" ")))
        }
        "matrix" => {
            let op = require_str(args, "operation")?;
            let rows = args
                .get("matrix")
                .and_then(|v| v.as_array())
                .ok_or_else(|| MathError::Eval("matrix needs a 'matrix' array".into()))?;
            // mathr's matrix syntax separates rows with `|` and entries with
            // commas or whitespace.
            let mut lines = Vec::with_capacity(rows.len());
            for row in rows {
                let cells = row.as_array().ok_or_else(|| {
                    MathError::Eval("matrix must be an array of rows".into())
                })?;
                let mut parts = Vec::with_capacity(cells.len());
                for c in cells {
                    let n = c.as_f64().ok_or_else(|| {
                        MathError::Eval("matrix entries must be numbers".into())
                    })?;
                    parts.push(crate::repl::format_value(n));
                }
                lines.push(parts.join(", "));
            }
            run_command(&format!("{} {}", op, lines.join(" | ")))
        }
        "minimize" => {
            let expr = require_str(args, "expression")?;
            // The bounds may be numbers or expressions ("pi", "pi/2").
            let from = bound_text(args, "from")?;
            let to = bound_text(args, "to")?;
            let cmd = format!("minimize {} from {} to {}", expr, from, to);
            run_with_var(cmd, args)
        }
        other => Err(MathError::Eval(format!("unknown tool: {}", other))),
    }
}

/// Read a bound argument as source text.
///
/// A number is formatted the REPL way; a string is passed through so
/// expressions like `pi` or `pi/2` work.
fn bound_text(args: &Json, key: &str) -> Result<String> {
    match args.get(key) {
        Some(Json::Num(n)) => Ok(crate::repl::format_value(*n)),
        Some(Json::Str(s)) => Ok(s.trim().to_string()),
        _ => Err(MathError::Eval(format!("missing required argument '{}'", key))),
    }
}

/// Fetch a required string argument.
fn require_str(args: &Json, key: &str) -> Result<String> {
    args.get(key)
        .and_then(|v| v.as_str())
        .map(|s| s.to_string())
        .ok_or_else(|| MathError::Eval(format!("missing required argument '{}'", key)))
}

/// Fetch a required numeric argument.
fn require_num(args: &Json, key: &str) -> Result<f64> {
    args.get(key)
        .and_then(|v| v.as_f64())
        .ok_or_else(|| MathError::Eval(format!("missing required argument '{}'", key)))
}

/// Run a command, optionally appending `for <variable>`.
fn run_with_var(mut cmd: String, args: &Json) -> Result<String> {
    if let Some(var) = args.get("variable").and_then(|v| v.as_str()) {
        if !var.is_empty() {
            cmd.push_str(&format!(" for {}", var));
        }
    }
    run_command(&cmd)
}

/// Run a mathr command string and return its output.
fn run_command(cmd: &str) -> Result<String> {
    let ctx = crate::eval::Context::standard();
    match crate::repl::dispatch_str(cmd, ctx)? {
        Some(out) => Ok(out),
        None => Ok(String::new()),
    }
}
