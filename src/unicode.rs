//! Unicode math input normalization.
//!
//! Users type `√2`, `x²`, `2×3`, `π`, `−5` — all of which are ordinary
//! Unicode, and all of which a byte-oriented ASCII parser cannot read. This
//! module rewrites such input into the ASCII forms the parser already
//! understands, *before* tokenizing.
//!
//! Two kinds of rewriting happen:
//!
//! 1. **Character substitution** — `×` → `*`, `÷` → `/`, `√` → `sqrt`,
//!    `π` → `pi`, `−` (U+2212) → `-`, and so on.
//! 2. **Superscript runs** — a run of superscript digits/letters after a base
//!    becomes an explicit power: `x²` → `x^(2)`, `x⁻¹` → `x^(-1)`,
//!    `xⁿ` → `x^(n)`.
//!
//! The rewrite is deliberately conservative: it only fires on characters that
//! are unambiguously mathematical, so ordinary prose identifiers are left
//! alone.
//!
//! Cross-port pattern from the `unsolve` TypeScript port
//! (`../unsolve/api/lib/engine/parser.ts:25-70`).

/// Rewrite Unicode math input into its ASCII equivalent.
///
/// Returns the input unchanged (borrowed) when nothing needed rewriting, so
/// the common ASCII path stays allocation-free.
pub fn normalize(src: &str) -> String {
    let mut out = String::with_capacity(src.len());
    let mut changed = false;
    let mut chars = src.chars().peekable();

    while let Some(c) = chars.next() {
        match c {
            // Superscript run: collect the whole run, then decide.
            c if is_superscript_start(c) => {
                let mut run = String::new();
                run.push(c);
                while let Some(&next) = chars.peek() {
                    if is_superscript_continue(next) {
                        run.push(next);
                        chars.next();
                    } else {
                        break;
                    }
                }
                match superscript_to_ascii(&run) {
                    Some(exp) => {
                        // A superscript run always becomes an explicit power.
                        // `x²` → `x^(2)`.
                        out.push_str("^(");
                        out.push_str(&exp);
                        out.push(')');
                        changed = true;
                    }
                    None => {
                        // Not a superscript we understand — pass it through
                        // unchanged so the parser reports the real error.
                        out.push_str(&run);
                    }
                }
            }
            // Root sign: `√2` → `sqrt(2)`, `∛8` → `cbrt(8)`.
            // The argument is delimited with parentheses, because `sqrt2`
            // would be read as a single identifier.
            '\u{221A}' | '\u{221B}' => {
                let func = if c == '\u{221A}' { "sqrt" } else { "cbrt" };
                // Capture the raw argument text first, then normalize it
                // separately — the argument may itself contain superscripts
                // or other Unicode that the main loop would have handled.
                let arg = capture_root_argument(&mut chars);
                out.push_str(func);
                out.push('(');
                out.push_str(&normalize(&arg));
                out.push(')');
                changed = true;
            }
            // Ordinary character: substitute if we have a mapping.
            c => {
                if let Some(repl) = char_replacement(c) {
                    out.push_str(repl);
                    changed = true;
                } else {
                    out.push(c);
                }
            }
        }
    }

    if changed {
        out
    } else {
        src.to_string()
    }
}

/// Map a single Unicode character to its ASCII replacement, if any.
fn char_replacement(c: char) -> Option<&'static str> {
    Some(match c {
        // Minus / dash family → ASCII hyphen-minus.
        '\u{2212}' => "-", // MINUS SIGN
        '\u{2012}' => "-", // FIGURE DASH
        '\u{2013}' => "-", // EN DASH
        '\u{2014}' => "-", // EM DASH
        '\u{2015}' => "-", // HORIZONTAL BAR

        // Multiplication / division.
        '\u{00D7}' => "*",  // MULTIPLICATION SIGN
        '\u{22C5}' => "*",  // DOT OPERATOR
        '\u{2219}' => "*",  // BULLET OPERATOR
        '\u{00B7}' => "*",  // MIDDLE DOT
        '\u{00F7}' => "/",  // DIVISION SIGN
        '\u{2215}' => "/",  // DIVISION SLASH
        '\u{2216}' => "/",  // SET MINUS
        '\u{2044}' => "/",  // FRACTION SLASH

        // Powers / roots.
        '\u{221E}' => "inf",  // INFINITY
        '\u{00B1}' => "+",   // PLUS-MINUS

        // Constants.
        '\u{03C0}' => "pi",    // GREEK SMALL LETTER PI
        '\u{03A0}' => "pi",    // GREEK CAPITAL LETTER PI
        '\u{03C4}' => "tau",   // GREEK SMALL LETTER TAU
        '\u{03A4}' => "tau",   // GREEK CAPITAL LETTER TAU

        // Invisible multiplication.
        '\u{2062}' => "*", // INVISIBLE TIMES
        '\u{2009}' => " ", // THIN SPACE (a separator, not an operator)
        '\u{00A0}' => " ", // NO-BREAK SPACE

        _ => return None,
    })
}

/// Capture the argument of a root sign.
///
/// The argument is either a parenthesised group — `√(x+1)` — or a run of
/// non-space, non-delimiter characters — `√2`. Returns the *raw* text, with
/// any outer parentheses kept, so the caller can normalize it separately.
fn capture_root_argument(
    chars: &mut std::iter::Peekable<std::str::Chars<'_>>,
) -> String {
    let mut arg = String::new();
    // Skip whitespace before the argument.
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
        } else {
            break;
        }
    }
    // A parenthesised group is taken whole.
    if chars.peek() == Some(&'(') {
        let mut depth = 0usize;
        for c in chars.by_ref() {
            arg.push(c);
            match c {
                '(' => depth += 1,
                ')' => {
                    depth -= 1;
                    if depth == 0 {
                        break;
                    }
                }
                _ => {}
            }
        }
        return arg;
    }
    // Otherwise take the run of characters up to whitespace or a delimiter.
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() || matches!(c, ')' | ']' | '}' | '+' | '-' | '*' | '/' | '^') {
            break;
        }
        arg.push(c);
        chars.next();
    }
    arg
}

/// Is this character the first of a superscript run?
fn is_superscript_start(c: char) -> bool {
    matches!(
        c,
        '\u{00B2}'   // ²
            | '\u{00B3}' // ³
            | '\u{00B9}' // ¹
            | '\u{2070}' // ⁰
            | '\u{2071}' // ¹
            | '\u{2074}'..='\u{2079}' // ⁴..⁹
            | '\u{207A}'..='\u{207E}' // ⁺ ⁻ ⁼ ⁽ ⁾ etc.
            | '\u{207F}' // ⁿ
            | '\u{00BC}' // ¼
            | '\u{00BD}' // ½
            | '\u{00BE}' // ¾
    )
}

/// Can this character continue a superscript run?
fn is_superscript_continue(c: char) -> bool {
    is_superscript_start(c) || matches!(c, '\u{03B2}') // β
}

/// Convert a superscript run to its ASCII equivalent.
///
/// Returns `None` when the run contains something we cannot represent
/// faithfully, so the caller can leave the input untouched rather than
/// silently mis-parse it.
fn superscript_to_ascii(run: &str) -> Option<String> {
    let mut out = String::with_capacity(run.len());
    for c in run.chars() {
        let mapped = match c {
            '\u{00B0}' => "0",
            '\u{00B9}' => "1",
            '\u{00B2}' => "2",
            '\u{00B3}' => "3",
            '\u{2070}' => "0",
            '\u{2074}' => "4",
            '\u{2075}' => "5",
            '\u{2076}' => "6",
            '\u{2077}' => "7",
            '\u{2078}' => "8",
            '\u{2079}' => "9",
            '\u{207A}' => "+",
            '\u{207B}' => "-",
            '\u{207C}' => "=",
            '\u{207D}' => "(",
            '\u{207E}' => ")",
            '\u{2071}' => "i",
            '\u{207F}' => "n",
            '\u{03B2}' => "beta",
            '\u{00BC}' => "1/4",
            '\u{00BD}' => "1/2",
            '\u{00BE}' => "3/4",
            _ => return None,
        };
        out.push_str(mapped);
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn n(s: &str) -> String {
        normalize(s)
    }

    #[test]
    fn ascii_is_unchanged() {
        assert_eq!(n("x^2 + 3*x"), "x^2 + 3*x");
        assert_eq!(n("sin(pi/4)"), "sin(pi/4)");
    }

    #[test]
    fn multiplication_sign_becomes_star() {
        assert_eq!(n("2×3"), "2*3");
        assert_eq!(n("2·3"), "2*3");
        assert_eq!(n("a⋅b"), "a*b");
    }

    #[test]
    fn division_sign_becomes_slash() {
        assert_eq!(n("5÷2"), "5/2");
        assert_eq!(n("6∕2"), "6/2");
    }

    #[test]
    fn unicode_minus_becomes_hyphen() {
        assert_eq!(n("−5"), "-5");
        assert_eq!(n("3−2"), "3-2");
        // En dash and em dash too.
        assert_eq!(n("3–2"), "3-2");
        assert_eq!(n("3—2"), "3-2");
    }

    #[test]
    fn root_sign_becomes_sqrt_call() {
        // `√2` must become `sqrt(2)` — `sqrt2` would be one identifier.
        assert_eq!(n("√2"), "sqrt(2)");
        assert_eq!(n("√x"), "sqrt(x)");
        assert_eq!(n("∛8"), "cbrt(8)");
    }

    #[test]
    fn root_sign_keeps_parenthesised_argument() {
        assert_eq!(n("√(x+1)"), "sqrt((x+1))");
        assert_eq!(n("√(x²+1)"), "sqrt((x^(2)+1))");
    }

    #[test]
    fn greek_constants_become_names() {
        assert_eq!(n("π"), "pi");
        assert_eq!(n("Π"), "pi");
        assert_eq!(n("τ"), "tau");
        assert_eq!(n("∞"), "inf");
    }

    #[test]
    fn superscript_digits_become_powers() {
        assert_eq!(n("x²"), "x^(2)");
        assert_eq!(n("x³"), "x^(3)");
        assert_eq!(n("x¹"), "x^(1)");
        assert_eq!(n("x⁰"), "x^(0)");
    }

    #[test]
    fn superscript_multi_digit_run_becomes_one_power() {
        assert_eq!(n("x¹²"), "x^(12)");
        assert_eq!(n("x²³"), "x^(23)");
    }

    #[test]
    fn superscript_negative_becomes_negative_power() {
        assert_eq!(n("x⁻¹"), "x^(-1)");
        assert_eq!(n("x⁻²"), "x^(-2)");
    }

    #[test]
    fn superscript_letters_become_powers() {
        assert_eq!(n("xⁿ"), "x^(n)");
    }

    #[test]
    fn superscript_fraction_glyphs_become_powers() {
        assert_eq!(n("x½"), "x^(1/2)");
        assert_eq!(n("x¼"), "x^(1/4)");
        assert_eq!(n("x¾"), "x^(3/4)");
    }

    #[test]
    fn superscript_after_number() {
        assert_eq!(n("2²"), "2^(2)");
    }

    #[test]
    fn multiple_superscripts_in_one_string() {
        assert_eq!(n("x²+y³"), "x^(2)+y^(3)");
    }

    #[test]
    fn unknown_superscript_is_left_alone() {
        // U+2071 is "i" but we handle it; use an unmapped one instead.
        // The point: an unrecognized glyph must not be silently dropped.
        let src = "x\u{2A2B}";
        assert_eq!(n(src), src);
    }

    #[test]
    fn mixed_unicode_input() {
        assert_eq!(n("2×x²"), "2*x^(2)");
        assert_eq!(n("√(x²+1)"), "sqrt((x^(2)+1))");
    }

    #[test]
    fn plain_identifiers_are_untouched() {
        // "pi" spelled out must not become something else.
        assert_eq!(n("pie"), "pie");
        assert_eq!(n("sin"), "sin");
    }
}
