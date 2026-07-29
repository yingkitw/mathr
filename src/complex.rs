use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// A small generic complex number type used by the FFT module.
/// Implemented from scratch (no external crates) to keep the dependency
/// footprint minimal and to demonstrate the FFT pipeline end-to-end.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex<T> {
    pub re: T,
    pub im: T,
}

impl<T> Complex<T> {
    pub const fn new(re: T, im: T) -> Self {
        Self { re, im }
    }
}

impl<T: Copy + Neg<Output = T>> Complex<T> {
    pub fn conj(&self) -> Self {
        Self::new(self.re, -self.im)
    }
}

impl<T: Copy + Neg<Output = T>> Neg for Complex<T> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.re, -self.im)
    }
}

impl<T: Copy + Add<Output = T>> Add for Complex<T> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl<T: Copy + Sub<Output = T>> Sub for Complex<T> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T>> Mul for Complex<T> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl<T: Copy + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T>> Div
    for Complex<T>
{
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

impl<T: fmt::Display> fmt::Display for Complex<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} + {}i", self.re, self.im)
    }
}

impl<T: Default> Default for Complex<T> {
    fn default() -> Self {
        Self {
            re: T::default(),
            im: T::default(),
        }
    }
}

// Numerical traits for f64 Complex
impl Complex<f64> {
    pub const ZERO: Self = Self::new(0.0, 0.0);
    pub const ONE: Self = Self::new(1.0, 0.0);

    pub fn abs(&self) -> f64 {
        (self.re * self.re + self.im * self.im).sqrt()
    }

    pub fn arg(&self) -> f64 {
        self.im.atan2(self.re)
    }

    pub fn exp(&self) -> Self {
        let r = self.re.exp();
        Self::new(r * self.im.cos(), r * self.im.sin())
    }

    pub fn from_polar(r: f64, theta: f64) -> Self {
        Self::new(r * theta.cos(), r * theta.sin())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn arithmetic() {
        let a = Complex::new(1.0, 2.0);
        let b = Complex::new(3.0, 4.0);
        assert_eq!(a + b, Complex::new(4.0, 6.0));
        assert_eq!(b - a, Complex::new(2.0, 2.0));
        assert_eq!(a * b, Complex::new(-5.0, 10.0));
        let q = a / b;
        assert!(approx(q.re, 0.44));
        assert!(approx(q.im, 0.08));
    }

    #[test]
    fn exp_and_polar() {
        // z = iπ/2 (a complex number whose real part is 0, imaginary part is π/2)
        let z = Complex::new(0.0, std::f64::consts::FRAC_PI_2);
        // exp(iπ/2) = cos(π/2) + i sin(π/2) = 0 + i·1 = i
        let w = z.exp();
        // cos(pi/2) is ~6e-17 due to fp rounding
        assert!(w.re.abs() < 1e-9, "expected ~0, got {}", w.re);
        assert!(approx(w.im, 1.0));

        // round-trip: from_polar then exp gives e^{r·e^{iθ}}, not e^{iθ}.
        // Verify the magnitude.
        let r = 1.5;
        let theta = 0.7;
        let z2 = Complex::from_polar(r, theta); // z2 = r·e^{iθ}
        let w2 = z2.exp();                       // w2 = e^{z2} = e^{r·cos(θ)} · e^{i·r·sin(θ)}
        let expected_re = (r * theta.cos()).exp() * (r * theta.sin()).cos();
        let expected_im = (r * theta.cos()).exp() * (r * theta.sin()).sin();
        assert!(approx(w2.re, expected_re));
        assert!(approx(w2.im, expected_im));
    }

    #[test]
    fn conj_and_neg() {
        let z = Complex::new(1.5, -2.5);
        assert_eq!(z.conj(), Complex::new(1.5, 2.5));
        assert_eq!(-z, Complex::new(-1.5, 2.5));
    }

    #[test]
    fn default_works() {
        let z: Complex<f64> = Complex::default();
        assert_eq!(z, Complex::new(0.0, 0.0));
    }
}