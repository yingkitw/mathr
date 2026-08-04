use mathr::laurent;

fn main() {
    println!("=== Laurent Series Expansion Demo ===\n");

    // 1. Simple pole: f(x) = 1/x
    let ls = laurent::laurent_series_str("1/x", "x", 0.0, 1, 3).unwrap();
    println!("f(x) = 1/x  (pole of order 1 at x=0)");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-1 = {:.6}", ls.coeff(-1));
    println!("  c_0  = {:.6}", ls.coeff(0));
    println!();

    // 2. Double pole: f(x) = 1/x^2
    let ls = laurent::laurent_series_str("1/x^2", "x", 0.0, 2, 3).unwrap();
    println!("f(x) = 1/x^2  (pole of order 2 at x=0)");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-2 = {:.6}", ls.coeff(-2));
    println!();

    // 3. Simple pole with analytic part: f(x) = 1/x + 1 + x
    let ls = laurent::laurent_series_str("1/x + 1 + x", "x", 0.0, 1, 3).unwrap();
    println!("f(x) = 1/x + 1 + x");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-1 = {:.6}, c_0 = {:.6}, c_1 = {:.6}",
        ls.coeff(-1), ls.coeff(0), ls.coeff(1));
    println!();

    // 4. exp(x)/x: 1/x + 1 + x/2 + x^2/6 + ...
    let ls = laurent::laurent_series_str("exp(x)/x", "x", 0.0, 1, 5).unwrap();
    println!("f(x) = exp(x)/x");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-1 = {:.6}, c_0 = {:.6}, c_1 = {:.6}, c_2 = {:.6}",
        ls.coeff(-1), ls.coeff(0), ls.coeff(1), ls.coeff(2));
    let val = ls.eval(1.0);
    println!("  eval at x=1: {:.6}  (exact: e = {:.6})", val, std::f64::consts::E);
    println!();

    // 5. Rational function: 1/(x(1-x)) = 1/x + 1 + x + x^2 + ...
    let ls = laurent::laurent_series_str("1/(x*(1-x))", "x", 0.0, 1, 5).unwrap();
    println!("f(x) = 1/(x(1-x))");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-1 = {:.6}, c_0 = {:.6}, c_1 = {:.6}, c_2 = {:.6}",
        ls.coeff(-1), ls.coeff(0), ls.coeff(1), ls.coeff(2));
    println!();

    // 6. Expansion around a non-zero center: f(x) = 1/(x-2), pole at x=2
    let ls = laurent::laurent_series_str("1/(x-2)", "x", 2.0, 1, 3).unwrap();
    println!("f(x) = 1/(x-2)  (pole of order 1 at x=2)");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_-1 = {:.6}", ls.coeff(-1));
    let val = ls.eval(3.0);
    println!("  eval at x=3: {:.6}  (exact: 1.0)", val);
    println!();

    // 7. No pole (Taylor-like): f(x) = cos(x)
    let ls = laurent::laurent_series_str("cos(x)", "x", 0.0, 0, 6).unwrap();
    println!("f(x) = cos(x)  (no pole)");
    println!("  Laurent series: {}", ls.to_string());
    println!("  c_0 = {:.6}, c_2 = {:.6}", ls.coeff(0), ls.coeff(2));
    let val = ls.eval(0.5);
    println!("  eval at x=0.5: {:.6}  (exact: {:.6})", val, 0.5_f64.cos());

    println!("\n=== Done ===");
}
