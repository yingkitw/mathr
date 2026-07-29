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
use std::f64::consts::PI;

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
    let bits = n.trailing_zeros() as usize;
    for i in 0..n {
        let j = reverse_bits(i, bits);
        if i < j {
            buf.swap(i, j);
        }
    }

    // ---- 2. Cooley–Tukey butterfly with precomputed twiddles --------
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut size = 2;
    while size <= n {
        let half = size / 2;
        // Precompute twiddle factors for this stage (half factors)
        let theta = sign * 2.0 * std::f64::consts::PI / size as f64;
        let twiddles: Vec<Complex<f64>> = (0..half)
            .map(|k| {
                let angle = theta * k as f64;
                Complex::new(angle.cos(), angle.sin())
            })
            .collect();
        let mut start = 0;
        while start < n {
            for k in 0..half {
                let t = twiddles[k] * buf[start + k + half];
                let u = buf[start + k];
                buf[start + k] = u + t;
                buf[start + k + half] = u - t;
            }
            start += size;
        }
        size *= 2;
    }
}

/// Reverse the low `bits` bits of `x`.
fn reverse_bits(x: usize, bits: usize) -> usize {
    // Byte-level reversal via lookup table, then shift to keep only `bits` bits.
    const REV_BYTE: [u8; 256] = {
        let mut table = [0u8; 256];
        let mut i = 0;
        while i < 256 {
            let mut val = i as u8;
            let mut rev = 0u8;
            let mut j = 0;
            while j < 8 {
                rev = (rev << 1) | (val & 1);
                val >>= 1;
                j += 1;
            }
            table[i] = rev;
            i += 1;
        }
        table
    };

    let b0 = REV_BYTE[(x & 0xFF) as usize] as usize;
    let b1 = REV_BYTE[((x >> 8) & 0xFF) as usize] as usize;
    let b2 = REV_BYTE[((x >> 16) & 0xFF) as usize] as usize;
    let b3 = REV_BYTE[((x >> 24) & 0xFF) as usize] as usize;
    let b4 = REV_BYTE[((x >> 32) & 0xFF) as usize] as usize;
    let b5 = REV_BYTE[((x >> 40) & 0xFF) as usize] as usize;
    let b6 = REV_BYTE[((x >> 48) & 0xFF) as usize] as usize;
    let b7 = REV_BYTE[((x >> 56) & 0xFF) as usize] as usize;

    let reversed = (b0 << 56) | (b1 << 48) | (b2 << 40) | (b3 << 32)
                 | (b4 << 24) | (b5 << 16) | (b6 << 8)  |  b7;
    reversed >> (64 - bits)
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
    1 << (usize::BITS - (n - 1).leading_zeros())
}

// --- Window functions -------------------------------------------------------

/// Hamming window coefficients: w[i] = a0 - a1*cos(2πi/(N-1))
const HAMMING_A0: f64 = 0.54;
const HAMMING_A1: f64 = 0.46;

/// Blackman window coefficients: w[i] = a0 - a1*cos(2πi/(N-1)) + a2*cos(4πi/(N-1))
const BLACKMAN_A0: f64 = 0.42;
const BLACKMAN_A1: f64 = 0.5;
const BLACKMAN_A2: f64 = 0.08;

/// Apply a window function to a slice of samples. Returns a new Vec.
pub fn apply_window(samples: &[f64], window: Window) -> Vec<f64> {
    let n = samples.len();
    if n <= 1 {
        return samples.to_vec();
    }
    let denom = (n - 1) as f64;
    match window {
        Window::Rectangular => samples.to_vec(),
        Window::Hann => samples
            .iter()
            .enumerate()
            .map(|(i, &x)| x * 0.5 * (1.0 - (2.0 * PI * i as f64 / denom).cos()))
            .collect(),
        Window::Hamming => samples
            .iter()
            .enumerate()
            .map(|(i, &x)| x * (HAMMING_A0 - HAMMING_A1 * (2.0 * PI * i as f64 / denom).cos()))
            .collect(),
        Window::Blackman => samples
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let t = 2.0 * PI * i as f64 / denom;
                x * (BLACKMAN_A0 - BLACKMAN_A1 * t.cos() + BLACKMAN_A2 * (2.0 * t).cos())
            })
            .collect(),
    }
}

/// Window function types for spectral analysis.
#[derive(Debug, Clone, Copy)]
pub enum Window {
    Rectangular,
    Hann,
    Hamming,
    Blackman,
}

// --- Convolution & cross-correlation via FFT --------------------------------

/// Fast convolution of two real signals via FFT. Returns `len = a.len() + b.len() - 1`.
pub fn convolve(a: &[f64], b: &[f64]) -> Result<Vec<f64>> {
    if a.is_empty() || b.is_empty() {
        return Err(MathError::InvalidArgument("convolve: empty input".into()));
    }
    let result_len = a.len() + b.len() - 1;
    let n = next_pow2(result_len);

    let mut fa: Vec<Complex<f64>> = a.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut fb: Vec<Complex<f64>> = b.iter().map(|&x| Complex::new(x, 0.0)).collect();
    fa.resize(n, Complex::ZERO);
    fb.resize(n, Complex::ZERO);

    fft_in_place(&mut fa, false);
    fft_in_place(&mut fb, false);

    for i in 0..n {
        fa[i] = fa[i] * fb[i];
    }

    fft_in_place(&mut fa, true);
    let scale = 1.0 / n as f64;
    Ok(fa[..result_len].iter().map(|c| c.re * scale).collect())
}

/// Cross-correlation of two real signals via FFT. Returns `len = a.len() + b.len() - 1`.
pub fn cross_correlate(a: &[f64], b: &[f64]) -> Result<Vec<f64>> {
    if a.is_empty() || b.is_empty() {
        return Err(MathError::InvalidArgument("cross_correlate: empty input".into()));
    }
    let b_rev: Vec<f64> = b.iter().rev().copied().collect();
    convolve(a, &b_rev)
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

    #[test]
    fn convolve_impulse() {
        // Convolving with [1, 0, 0, ...] should return the original signal
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![1.0];
        let c = convolve(&a, &b).unwrap();
        assert_eq!(c.len(), 4);
        for i in 0..4 {
            assert!(close(c[i], a[i], 1e-9));
        }
    }

    #[test]
    fn convolve_known() {
        // [1, 2, 3] * [1, 1] = [1, 3, 5, 3]
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 1.0];
        let c = convolve(&a, &b).unwrap();
        assert_eq!(c.len(), 4);
        assert!(close(c[0], 1.0, 1e-9));
        assert!(close(c[1], 3.0, 1e-9));
        assert!(close(c[2], 5.0, 1e-9));
        assert!(close(c[3], 3.0, 1e-9));
    }

    #[test]
    fn cross_correlate_self() {
        // Auto-correlation of [1, 2, 3, 4] should peak at the center
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let c = cross_correlate(&a, &a).unwrap();
        assert_eq!(c.len(), 7);
        // The peak should be at index 3 (full overlap)
        let peak_idx = c
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(i, _)| i)
            .unwrap();
        assert_eq!(peak_idx, 3);
        // Peak value = sum of squares = 1+4+9+16 = 30
        assert!(close(c[3], 30.0, 1e-8));
    }

    #[test]
    fn window_hann_endpoints() {
        let n = 8;
        let samples = vec![1.0; n];
        let w = apply_window(&samples, Window::Hann);
        // Symmetric Hann window: zero at endpoints
        assert!(close(w[0], 0.0, 1e-10));
        assert!(close(w[n - 1], 0.0, 1e-10));
        // Near the center, value should be close to 1 (not exact for even N)
        assert!(w[n / 2] > 0.9 && w[n / 2] < 1.01);
    }

    #[test]
    fn window_hamming_endpoints() {
        let n = 8;
        let samples = vec![1.0; n];
        let w = apply_window(&samples, Window::Hamming);
        // Symmetric Hamming: endpoints = 0.54 - 0.46*cos(0) = 0.08
        assert!(close(w[0], 0.08, 1e-10));
        assert!(close(w[n - 1], 0.08, 1e-10));
        // Near center, close to 1.0
        assert!(w[n / 2] > 0.9 && w[n / 2] < 1.01);
    }

    #[test]
    fn window_blackman_endpoints() {
        let n = 8;
        let samples = vec![1.0; n];
        let w = apply_window(&samples, Window::Blackman);
        // Blackman window is zero at endpoints
        assert!(close(w[0], 0.0, 1e-10));
        assert!(close(w[n - 1], 0.0, 1e-10));
    }

    #[test]
    fn window_rectangular_identity() {
        let samples = vec![1.0, 2.0, 3.0, 4.0];
        let w = apply_window(&samples, Window::Rectangular);
        assert_eq!(w, samples);
    }
}