//! Number-theory command handlers.
//!
//! Handlers moved verbatim from the original monolithic `repl.rs`.

use super::*;
pub(super) fn do_numtheory(rest: &str, op: &str) -> Result<Option<String>> {
    use crate::numtheory;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    match op {
        "gcd" => {
            let nums = parse_u64_list(&tokens)?;
            if nums.len() < 2 { return Err(crate::error::MathError::Eval("gcd needs ≥2 numbers".into())); }
            let g = nums.iter().copied().reduce(numtheory::gcd).unwrap();
            Ok(Some(g.to_string()))
        }
        "lcm" => {
            let nums = parse_u64_list(&tokens)?;
            if nums.len() < 2 { return Err(crate::error::MathError::Eval("lcm needs ≥2 numbers".into())); }
            let l = nums.iter().copied().reduce(numtheory::lcm).unwrap();
            Ok(Some(l.to_string()))
        }
        "is-prime" => {
            // Auto-upgrade to BigInt when the argument does not fit in u64
            // (the same pattern `fib`/`binom`/`fact` use).
            match parse_u64_single(&tokens) {
                Ok(n) => Ok(Some(numtheory::is_prime(n).to_string())),
                Err(_) => {
                    let n = crate::bigint::parse(tokens.first().map(|s| s.trim()).unwrap_or(""))?;
                    Ok(Some(crate::bigint::is_prime(&n, 20).to_string()))
                }
            }
        }
        "factor" => {
            // Auto-upgrade to BigInt when the argument does not fit in u64.
            match parse_u64_single(&tokens) {
                Ok(n) => {
                    let factors = numtheory::prime_factors(n);
                    if factors.is_empty() {
                        Ok(Some("1".into()))
                    } else {
                        // Group repeated factors with exponents: 360 = 2^3 * 3^2 * 5
                        let mut grouped: Vec<(u64, u32)> = Vec::new();
                        for &p in &factors {
                            if let Some(last) = grouped.last_mut() {
                                if last.0 == p { last.1 += 1; continue; }
                            }
                            grouped.push((p, 1));
                        }
                        let strs: Vec<String> = grouped
                            .iter()
                            .map(|(p, e)| if *e == 1 { p.to_string() } else { format!("{}^{}", p, e) })
                            .collect();
                        Ok(Some(strs.join(" · ")))
                    }
                }
                Err(_) => {
                    let n = crate::bigint::parse(tokens.first().map(|s| s.trim()).unwrap_or(""))?;
                    let factors = crate::bigint::factorize(&n);
                    if factors.is_empty() {
                        return Ok(Some("1".into()));
                    }
                    let strs: Vec<String> = factors
                        .iter()
                        .map(|(p, e)| {
                            if *e == 1 {
                                p.to_string()
                            } else {
                                format!("{}^{}", p, e)
                            }
                        })
                        .collect();
                    Ok(Some(strs.join(" · ")))
                }
            }
        }
        "fib" => {
            let n = parse_u64_single(&tokens)?;
            // F(90) ≈ 2.88e18 fits in u64; F(91) overflows. Auto-upgrade to BigInt.
            if n <= 90 {
                Ok(Some(numtheory::fibonacci(n).to_string()))
            } else {
                Ok(Some(crate::bigint::fibonacci(n).to_string()))
            }
        }
        "binom" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("binom needs n k".into())); }
            let n: u64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
            let k: u64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad k".into()))?;
            // Try u64 first; fall back to BigInt on overflow (matching SymPy behavior).
            match numtheory::binomial(n, k) {
                Ok(r) => Ok(Some(r.to_string())),
                Err(_) => Ok(Some(crate::bigint::binomial(n, k).to_string())),
            }
        }
        "fact" => {
            let n = parse_u64_single(&tokens)?;
            // 20! fits in u64; 21! overflows. Auto-upgrade to BigInt.
            if n <= 20 {
                let r = numtheory::factorial(n)?;
                Ok(Some(r.to_string()))
            } else {
                Ok(Some(crate::bigint::factorial(n).to_string()))
            }
        }
        "mr-prime" => {
            let n = parse_u64_single(&tokens)?;
            let rounds: usize = tokens.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
            Ok(Some(numtheory::is_prime_miller_rabin(n, rounds).to_string()))
        }
        "jacobi" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("jacobi needs a n".into())); }
            let a: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
            let n: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad n".into()))?;
            Ok(Some(numtheory::jacobi_symbol(a, n).to_string()))
        }
        "cf" => {
            if tokens.len() < 2 { return Err(crate::error::MathError::Eval("cf needs p q".into())); }
            let p: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad p".into()))?;
            let q: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad q".into()))?;
            let cf = numtheory::continued_fraction(p, q)?;
            let strs: Vec<String> = cf.iter().map(|v| v.to_string()).collect();
            Ok(Some(format!("[{}]", strs.join("; "))))
        }
        "diophantine" => {
            if tokens.len() < 3 { return Err(crate::error::MathError::Eval("diophantine needs a b c".into())); }
            let a: i64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad a".into()))?;
            let b: i64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad b".into()))?;
            let c: i64 = tokens[2].parse().map_err(|_| crate::error::MathError::Eval("bad c".into()))?;
            let (x, y) = numtheory::diophantine(a, b, c)?;
            Ok(Some(format!("x = {}, y = {} (a*x + b*y = {})", x, y, a * x + b * y)))
        }
        "dlog" => {
            if tokens.len() < 3 { return Err(crate::error::MathError::Eval("dlog needs g h p".into())); }
            let g: u64 = tokens[0].parse().map_err(|_| crate::error::MathError::Eval("bad g".into()))?;
            let h: u64 = tokens[1].parse().map_err(|_| crate::error::MathError::Eval("bad h".into()))?;
            let p: u64 = tokens[2].parse().map_err(|_| crate::error::MathError::Eval("bad p".into()))?;
            match numtheory::discrete_log(g, h, p) {
                Some(x) => Ok(Some(format!("x = {} (g^x mod p = {})", x, numtheory::mod_pow(g, x, p)))),
                None => Ok(Some("(no discrete log found)".into())),
            }
        }
        _ => Err(crate::error::MathError::Eval(format!("unknown op: {}", op))),
    }
}

pub(super) fn do_fast(rest: &str) -> Result<Option<String>> {
    // Format: "func x" — e.g. "fast sin 1.5"
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "fast needs: <func> <x>  (func: sin, cos, tan, exp, log, sqrt, pow)".into(),
        ));
    }
    let func = tokens[0];
    let x: f64 = tokens[1]
        .parse()
        .map_err(|_| crate::error::MathError::Eval("bad x".into()))?;
    let (fast_val, exact_val, name) = match func {
        "sin" => (crate::fastmath::fast_sin(x), x.sin(), "sin"),
        "cos" => (crate::fastmath::fast_cos(x), x.cos(), "cos"),
        "tan" => (crate::fastmath::fast_tan(x), x.tan(), "tan"),
        "exp" => (crate::fastmath::fast_exp(x), x.exp(), "exp"),
        "log" | "ln" => (crate::fastmath::fast_log(x), x.ln(), "ln"),
        "sqrt" => (crate::fastmath::fast_sqrt(x), x.sqrt(), "sqrt"),
        "pow" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("fast pow needs: <x> <y>".into()));
            }
            let y: f64 = tokens[2]
                .parse()
                .map_err(|_| crate::error::MathError::Eval("bad y".into()))?;
            let fv = crate::fastmath::fast_pow(x, y);
            let ev = x.powf(y);
            return Ok(Some(format!(
                "fast pow({}, {}) = {}  (exact: {}, err: {:e})",
                x, y, format_value(fv), format_value(ev), (fv - ev).abs()
            )));
        }
        _ => {
            return Err(crate::error::MathError::Eval(format!(
                "unknown fast func '{}' (try: sin, cos, tan, exp, log, sqrt, pow)",
                func
            )));
        }
    };
    Ok(Some(format!(
        "fast {}({}) = {}  (exact: {}, err: {:e})",
        name,
        x,
        format_value(fast_val),
        format_value(exact_val),
        (fast_val - exact_val).abs()
    )))
}

pub(super) fn do_big(rest: &str) -> Result<Option<String>> {
    use crate::bigint;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "big needs: <op> <args>  (op: prime, factor, gcd, lcm, modpow, totient)".into(),
        ));
    }
    let op = tokens[0];
    match op {
        "prime" | "is-prime" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big prime needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            let result = bigint::is_prime(&n, 20);
            Ok(Some(format!("{} is {}", n, if result { "prime" } else { "composite" })))
        }
        "factor" | "factorize" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big factor needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            let factors = bigint::factorize(&n);
            if factors.is_empty() {
                return Ok(Some(format!("{} = 1", n)));
            }
            let parts: Vec<String> = factors
                .iter()
                .map(|(p, e)| {
                    if *e == 1 {
                        p.to_string()
                    } else {
                        format!("{}^{}", p, e)
                    }
                })
                .collect();
            Ok(Some(format!("{} = {}", n, parts.join(" · "))))
        }
        "gcd" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("big gcd needs: <a> <b>".into()));
            }
            let a = bigint::parse(tokens[1])?;
            let b = bigint::parse(tokens[2])?;
            Ok(Some(format!("gcd({}, {}) = {}", a, b, bigint::gcd(&a, &b))))
        }
        "lcm" => {
            if tokens.len() < 3 {
                return Err(crate::error::MathError::Eval("big lcm needs: <a> <b>".into()));
            }
            let a = bigint::parse(tokens[1])?;
            let b = bigint::parse(tokens[2])?;
            Ok(Some(format!("lcm({}, {}) = {}", a, b, bigint::lcm(&a, &b))))
        }
        "modpow" => {
            if tokens.len() < 4 {
                return Err(crate::error::MathError::Eval("big modpow needs: <base> <exp> <mod>".into()));
            }
            let base = bigint::parse(tokens[1])?;
            let exp = bigint::parse(tokens[2])?;
            let m = bigint::parse(tokens[3])?;
            let result = bigint::mod_pow(&base, &exp, &m)?;
            Ok(Some(format!("{}^{} ≡ {} (mod {})", base, exp, result, m)))
        }
        "totient" | "phi" => {
            if tokens.len() < 2 {
                return Err(crate::error::MathError::Eval("big totient needs: <n>".into()));
            }
            let n = bigint::parse(tokens[1])?;
            Ok(Some(format!("φ({}) = {}", n, bigint::totient(&n))))
        }
        _ => Err(crate::error::MathError::Eval(format!(
            "unknown big op '{}' (try: prime, factor, gcd, lcm, modpow, totient)",
            op
        ))),
    }
}

/// `dec <expr> [prec <n>] [with <var>=<val>,...]` — arbitrary-precision
/// decimal evaluation. `prec` is the number of significant digits (default
/// 30, max 1000); `pi`, `e`, and `tau` are computed at full precision.
pub(super) fn do_dec(rest: &str) -> Result<Option<String>> {
    use crate::bigdec;
    const DEFAULT_PREC: usize = 30;
    let mut prec = DEFAULT_PREC;
    let mut work = rest.trim().to_string();
    if let Some(idx) = work.rfind(" prec ") {
        let tail = work[idx + 6..].trim();
        let n: usize = tail.parse().map_err(|_| {
            crate::error::MathError::Eval("dec: prec must be a positive integer".into())
        })?;
        if n == 0 || n > 1000 {
            return Err(crate::error::MathError::Eval(
                "dec: prec must be between 1 and 1000".into(),
            ));
        }
        prec = n;
        work.truncate(idx);
    }
    let mut assign_part: Option<String> = None;
    if let Some(idx) = work.find(" with ") {
        assign_part = Some(work[idx + 6..].to_string());
        work.truncate(idx);
    }
    let expr_src = work.trim();
    if expr_src.is_empty() {
        return Err(crate::error::MathError::Eval(
            "dec needs: <expr> [prec <n>] [with <var>=<val>,...]".into(),
        ));
    }
    let e = Parser::parse(expr_src)?;
    let mut vars: std::collections::HashMap<String, bigdecimal::BigDecimal> =
        std::collections::HashMap::new();
    if let Some(a) = assign_part {
        for assignment in split_assignments(&a) {
            let kv: Vec<&str> = assignment.splitn(2, '=').collect();
            if kv.len() != 2 {
                return Err(crate::error::MathError::Eval(format!(
                    "bad assignment `{}` (expected: var=value)",
                    assignment
                )));
            }
            let name = kv[0].trim().to_string();
            let val = bigdec::parse(kv[1].trim())?;
            vars.insert(name, val);
        }
    }
    let result = bigdec::eval_decimal_rounded(&e, &vars, prec)?;
    Ok(Some(result.to_string()))
}

/// `limit <expr> <var> <point>` or `limit <expr> <point>` — compute the
/// limit as the variable approaches the point (finite, `inf`, or `-inf`).
/// The variable may be omitted when the expression has exactly one.
pub(super) fn do_base(rest: &str) -> Result<Option<String>> {
    // `base <n> <target>` — convert <n> to the target base.
    // `base <n> from <base> to <base>` — convert between two bases.
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.is_empty() {
        return Err(crate::error::MathError::Eval(
            "base needs: <n> <target>   or   <n> from <base> to <base>".into(),
        ));
    }
    use crate::numbase;

    // The `from … to …` form: convert between two bases.
    let lower: Vec<String> = tokens.iter().map(|t| t.to_ascii_lowercase()).collect();
    if let Some(from_pos) = lower.iter().rposition(|t| t == "from") {
        let to_pos = lower
            .iter()
            .rposition(|t| t == "to")
            .ok_or_else(|| crate::error::MathError::Eval("base needs: <n> from <base> to <base>".into()))?;
        if to_pos <= from_pos + 1 || to_pos + 1 >= tokens.len() || from_pos == 0 {
            return Err(crate::error::MathError::Eval(
                "base needs: <n> from <base> to <base>".into(),
            ));
        }
        let value = tokens[0];
        let from = numbase::parse_base(tokens[from_pos + 1])?;
        let to = numbase::parse_base(tokens[to_pos + 1])?;
        let n = numbase::from_base(value, from)?;
        let converted = numbase::to_base(&n, to)?;
        return Ok(Some(format!("{}_{}", converted, to)));
    }

    // The simple form: `base <n> <target>`.
    if tokens.len() < 2 {
        return Err(crate::error::MathError::Eval(
            "base needs: <n> <target>   or   <n> from <base> to <base>".into(),
        ));
    }
    let value = tokens[0];
    let target = numbase::parse_base(tokens[1])?;
    // The value is read in the base implied by its own digits: a bare decimal
    // number is decimal, but a string like "FF" or "0b1011" is not.
    let n = numbase::from_base(value, 10)?;
    let converted = numbase::to_base(&n, target)?;
    Ok(Some(format!("{}_{}", converted, target)))
}

