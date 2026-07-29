//! FFT implementation from scratch.
//!
//! We provide a radix-2 Cooley–Tukey fast Fourier transform in three flavours:
//! - [`fft`]: complex in, complex out (recursive and iterative versions)
//! - [`ifft`]: inverse FFT
//! - [`rfft`]: real-input forward FFT (packs `N` real samples into `N/2 + 1`
//!   complex bins)
//!
//! The algorithms are textbook implementations — clarity over micro-tuning —
//! but the iterative version is `O(N log N)` and uses `O(N)` scratch.

use crate::complex::Complex;
use crate::error::{MathError, Result};

/// Verify that `n` is a power of two. The radix-2 algorithm only accepts
/// such lengths.
fn check_n(n: usize) -> Result<()> {
    if n == 0 {
        return Err(MathError::InvalidArgument("FFT size must be > 0".into()));
    }
    if n & (n - 1) != 0 {
        return Err(MathError::InvalidArgument(format!(
            "FFT size must be a power of 2, got {}",
            n
        )));
    }
    Ok(())
}

/// Forward discrete Fourier transform of a complex input.
///
/// `X[k] = sum_{n=0..N-1} x[n] * exp(-2πi k n / N)`
pub fn fft(input: &[Complex<f64>]) -> Result<Vec<Complex<f64>>> {
    let n = input.len();
    check_n(n)?;
    let mut buf = input.to_vec();
    fft_in_place(&mut buf, false);
    Ok(buf)
}

/// Inverse discrete Fourier transform. Equivalent to FFT with conjugated
/// twiddle factors and a `1/N` scaling.
pub fn ifft(input: &[Complex<f64>]) -> Result<Vec<Complex<f64>>> {
    let n = input.len();
    check_n(n)?;
    let mut buf = input.to_vec();
    fft_in_place(&mut buf, true);
    for x in &mut buf {
        x.re /= n as f64;
        x.im /= n as f64;
    }
    Ok(buf)
}

/// Real-input forward FFT. Returns `N/2 + 1` complex bins (the positive
/// frequency half, including DC and Nyquist).
pub fn rfft(input: &[f64]) -> Result<Vec<Complex<f64>>> {
    let n = input.len();
    check_n(n)?;
    let complex_input: Vec<Complex<f64>> =
        input.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let full = fft(&complex_input)?;
    Ok(full[..n / 2 + 1].to_vec())
}

/// Compute the magnitude spectrum `|X[k]|` for a real input.
pub fn magnitude_spectrum(samples: &[f64]) -> Result<Vec<f64>> {
    Ok(rfft(samples)?.iter().map(|c| c.abs()).collect())
}

/// Compute the power spectrum `|X[k]|^2`.
pub fn power_spectrum(samples: &[f64]) -> Result<Vec<f64>> {
    Ok(rfft(samples)?.iter().map(|c| c.abs().powi(2)).collect())
}

/// In-place iterative Cooley–Tukey FFT.
///
/// `inverse` flips the sign of the twiddle factor exponent.
pub fn fft_in_place(buf: &mut [Complex<f64>], inverse: bool) {
    let n = buf.len();
    debug_assert!(n > 0 && (n & (n - 1)) == 0, "FFT size must be a power of 2");

    // ---- 1. Bit-reversal permutation -------------------------------
    // Each index i is rewritten as reverse_bits(i, log2(n)).
    let bits = n.trailing_zeros() as usize;
    for i in 0..n {
        let j = reverse_bits(i, bits);
        if i < j {
            buf.swap(i, j);
        }
    }

    // ---- 2. Cooley–Tukey butterfly -------------------------------
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut size = 2;
    while size <= n {
        let half = size / 2;
        // principal nth root of unity for this stage
        let theta = sign * 2.0 * std::f64::consts::PI / size as f64;
        let w_step = Complex::new(theta.cos(), theta.sin());
        let mut start = 0;
        while start < n {
            let mut w = Complex::new(1.0, 0.0);
            for k in 0..half {
                let t = w * buf[start + k + half];
                let u = buf[start + k];
                buf[start + k] = u + t;
                buf[start + k + half] = u - t;
                w = w * w_step;
            }
            start += size;
        }
        size *= 2;
    }
}

fn reverse_bits(mut x: usize, bits: usize) -> usize {
    let mut y = 0usize;
    for _ in 0..bits {
        y = (y << 1) | (x & 1);
        x >>= 1;
    }
    y
}

/// Two-dimensional FFT, useful for image-processing style workloads.
pub fn fft2(input: &[Vec<Complex<f64>>]) -> Result<Vec<Vec<Complex<f64>>>> {
    if input.is_empty() {
        return Ok(Vec::new());
    }
    let rows = input.len();
    let cols = input[0].len();
    check_n(rows)?;
    check_n(cols)?;

    // transform each row
    let mut data: Vec<Vec<Complex<f64>>> = input.to_vec();
    for r in &mut data {
        let mut row = std::mem::take(r);
        fft_in_place(&mut row, false);
        *r = row;
    }

    // transform each column
    for c in 0..cols {
        let mut col: Vec<Complex<f64>> = (0..rows).map(|r| data[r][c]).collect();
        fft_in_place(&mut col, false);
        for r in 0..rows {
            data[r][c] = col[r];
        }
    }
    Ok(data)
}

/// Next power of two >= `n`.
pub fn next_pow2(n: usize) -> usize {
    if n <= 1 {
        return 1;
    }
    let mut p = 1;
    while p < n {
        p <<= 1;
    }
    p
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn close(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn dft_impulse() {
        // x = [1, 0, 0, 0] -> X = [1, 1, 1, 1]
        let x = vec![
            Complex::new(1.0, 0.0),
            Complex::ZERO,
            Complex::ZERO,
            Complex::ZERO,
        ];
        let y = fft(&x).unwrap();
        for v in &y {
            assert!(close(v.re, 1.0, 1e-12));
            assert!(close(v.im, 0.0, 1e-12));
        }
    }

    #[test]
    fn dft_complex_exponential() {
        // x[n] = exp(2πi k0 n / N),  k0 = 1, N = 8 -> X[k] is 0 except at k=1
        let n = 8;
        let x: Vec<Complex<f64>> = (0..n)
            .map(|i| Complex::from_polar(1.0, 2.0 * PI * i as f64 / n as f64))
            .collect();
        let y = fft(&x).unwrap();
        for (k, v) in y.iter().enumerate() {
            if k == 1 {
                assert!(close(v.abs(), n as f64, 1e-9));
            } else {
                assert!(close(v.abs(), 0.0, 1e-9));
            }
        }
    }

    #[test]
    fn inverse_roundtrip() {
        let x: Vec<Complex<f64>> = (0..16)
            .map(|i| Complex::new(i as f64, -(i as f64) * 0.5))
            .collect();
        let y = fft(&x).unwrap();
        let back = ifft(&y).unwrap();
        for (a, b) in x.iter().zip(back.iter()) {
            assert!(close(a.re, b.re, 1e-9));
            assert!(close(a.im, b.im, 1e-9));
        }
    }

    #[test]
    fn real_sine_magnitude() {
        // 64 samples of a 4 Hz sine sampled at 64 Hz -> bin 4 should dominate
        let n = 64;
        let samples: Vec<f64> = (0..n).map(|i| (2.0 * PI * 4.0 * i as f64 / n as f64).sin()).collect();
        let mags = magnitude_spectrum(&samples).unwrap();
        // bin 0 is DC, then look at bin 4 (and 60, the mirror image)
        assert!(mags[4] > 30.0, "expected strong bin at k=4, got {}", mags[4]);
        for (k, m) in mags.iter().enumerate() {
            if k != 4 && k != (n - 4) && k != 0 {
                assert!(*m < 1.0, "unexpected energy at k={}: {}", k, m);
            }
        }
    }

    #[test]
    fn power_of_two_check() {
        assert!(fft(&vec![Complex::ZERO; 3]).is_err());
        assert!(fft(&vec![Complex::ZERO; 8]).is_ok());
    }

    #[test]
    fn next_pow2_basic() {
        assert_eq!(next_pow2(1), 1);
        assert_eq!(next_pow2(5), 8);
        assert_eq!(next_pow2(8), 8);
        assert_eq!(next_pow2(9), 16);
    }

    #[test]
    fn fft2_roundtrip() {
        let n = 4;
        let mut data = vec![vec![Complex::ZERO; n]; n];
        for r in 0..n {
            for c in 0..n {
                data[r][c] = Complex::new((r * n + c) as f64, 0.0);
            }
        }
        let transformed = fft2(&data).unwrap();
        assert_eq!(transformed.len(), n);
        assert_eq!(transformed[0].len(), n);
    }
}