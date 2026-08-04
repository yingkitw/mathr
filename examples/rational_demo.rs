use mathr::rational::{parse_rational, Rational};

fn main() {
    println!("=== Rational Arithmetic Demo ===\n");

    // 1. Construction and reduction
    let a = Rational::new(6, 8).unwrap();
    println!("Rational::new(6, 8) = {}  (reduced to {}/{})", a, a.num(), a.den());

    let b = Rational::new(3, -4).unwrap();
    println!("Rational::new(3, -4) = {}  (denominator normalised positive)", b);
    println!();

    // 2. Arithmetic
    let half = parse_rational("1/2").unwrap();
    let third = parse_rational("1/3").unwrap();
    let quarter = parse_rational("1/4").unwrap();

    println!("1/2 + 1/3 = {}", half + third);
    println!("1/2 - 1/3 = {}", half - third);
    println!("1/2 * 1/3 = {}", half * third);
    println!("1/2 / 1/3 = {}", half / third);
    println!();

    // 3. Chained: (1/2 + 1/3) * 1/4 = 5/24
    let result = (half + third) * quarter;
    println!("(1/2 + 1/3) * 1/4 = {}", result);
    println!();

    // 4. Powers
    let two_thirds = Rational::new(2, 3).unwrap();
    println!("(2/3)^3 = {}", two_thirds.powi(3));
    println!("(2/3)^(-2) = {}", two_thirds.powi(-2));
    println!();

    // 5. Parsing: integers, fractions, decimals
    let inputs = ["42", "-7", "3/4", "-3/4", "0.5", "-1.25", "0.125"];
    for s in inputs {
        let r = parse_rational(s).unwrap();
        println!("parse_rational(\"{}\") = {}  (= {:.6} as f64)", s, r, r.to_f64());
    }
    println!();

    // 6. Comparison and ordering
    let r1 = Rational::new(1, 3).unwrap();
    let r2 = Rational::new(1, 2).unwrap();
    println!("1/3 < 1/2: {}", r1 < r2);
    println!("1/3 == 2/6: {}", r1 == Rational::new(2, 6).unwrap());
    println!();

    // 7. Large denominators (i128 intermediate arithmetic)
    let big1 = Rational::new(1, 1_000_000_000).unwrap();
    let big2 = Rational::new(1, 1_000_000_001).unwrap();
    println!("1/1e9 > 1/(1e9+1): {}", big1 > big2);
    let sum = big1 + big2;
    println!("1/1e9 + 1/(1e9+1) = {}", sum);
    println!();

    // 8. Absolute value and negation
    let neg = Rational::new(-3, 4).unwrap();
    println!("-3/4: abs = {}, neg = {}", neg.abs(), -neg);

    println!("\n=== Done ===");
}
