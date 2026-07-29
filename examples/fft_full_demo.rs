use mathr::complex::Complex;
use mathr::fft;
use mathr::stats;

fn main() {
    println!("=== FFT & Signal Processing Demo ===\n");

    // --- Basic FFT ---
    println!("--- FFT of a simple signal ---");
    let signal: Vec<f64> = (0..8)
        .map(|i| (2.0 * std::f64::consts::PI * i as f64 / 8.0).sin())
        .collect();
    let complex_signal: Vec<Complex<f64>> = signal.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let spectrum = fft::fft(&complex_signal).unwrap();
    println!("  Input: 8 samples of sin(2πk/8)");
    for (i, c) in spectrum.iter().enumerate() {
        println!("  X[{}] = {:.4} + {:.4}j  |X| = {:.4}", i, c.re, c.im, c.abs());
    }

    // --- Magnitude spectrum ---
    println!("\n--- Magnitude spectrum of a 64-point signal ---");
    let n = 64;
    let samples: Vec<f64> = (0..n)
        .map(|i| (2.0 * std::f64::consts::PI * 8.0 * i as f64 / n as f64).sin())
        .collect();
    let mags = fft::magnitude_spectrum(&samples).unwrap();
    let peak_bin = mags.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    println!("  Peak at bin {} (expected bin 8), magnitude = {:.2}", peak_bin, mags[peak_bin]);

    // --- Power spectrum ---
    println!("\n--- Power spectrum ---");
    let power = fft::power_spectrum(&samples).unwrap();
    println!("  Power at peak bin: {:.2}", power[peak_bin]);
    let total_power: f64 = power.iter().sum();
    println!("  Total power: {:.2}", total_power);

    // --- IFFT roundtrip ---
    println!("\n--- FFT → IFFT roundtrip ---");
    let original: Vec<Complex<f64>> = (0..16)
        .map(|i| Complex::new((i as f64 * 0.5).sin(), (i as f64 * 0.3).cos()))
        .collect();
    let freq = fft::fft(&original).unwrap();
    let recovered = fft::ifft(&freq).unwrap();
    let max_error = original.iter().zip(recovered.iter())
        .map(|(a, b)| (a.re - b.re).abs().max((a.im - b.im).abs()))
        .fold(0.0, f64::max);
    println!("  16-point roundtrip max error: {:.2e}", max_error);

    // --- Convolution ---
    println!("\n--- Convolution via FFT ---");
    let a = vec![1.0, 2.0, 3.0, 4.0];
    let b = vec![1.0, 1.0, 1.0, 1.0];
    let conv = fft::convolve(&a, &b).unwrap();
    println!("  conv([1,2,3,4], [1,1,1,1]) = {:?}", conv);

    // --- Cross-correlation ---
    println!("\n--- Cross-correlation ---");
    let x = vec![1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0];
    let y = vec![1.0, 2.0, 3.0, 4.0];
    let xcorr = fft::cross_correlate(&x, &y).unwrap();
    println!("  xcorr result length: {}", xcorr.len());
    let peak = xcorr.iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();
    println!("  Peak correlation at lag {} = {:.4}", peak.0, peak.1);

    // --- Window functions ---
    println!("\n--- Window functions ---");
    let signal = vec![1.0; 32];
    for (name, window) in [
        ("Hann", fft::apply_window(&signal, fft::Window::Hann)),
        ("Hamming", fft::apply_window(&signal, fft::Window::Hamming)),
        ("Blackman", fft::apply_window(&signal, fft::Window::Blackman)),
    ] {
        let coherent_gain: f64 = window.iter().sum::<f64>() / 32.0;
        println!("  {} window: coherent gain = {:.4}, first/last = {:.4}/{:.4}",
            name, coherent_gain, window[0], window[31]);
    }

    // --- Statistics on signal ---
    println!("\n--- Statistics on signal ---");
    let signal: Vec<f64> = (0..100)
        .map(|i| (2.0 * std::f64::consts::PI * 5.0 * i as f64 / 100.0).sin())
        .collect();
    let summary = stats::summary(&signal).unwrap();
    println!("  100 samples of sin(2π·5·k/100):");
    println!("    mean   = {:.6} (expected 0)", summary.mean);
    println!("    stddev = {:.6} (expected ~0.707)", summary.stddev);
    println!("    min    = {:.6}", summary.min);
    println!("    max    = {:.6}", summary.max);

    // --- Correlation ---
    println!("\n--- Linear correlation ---");
    let x: Vec<f64> = (0..20).map(|i| i as f64).collect();
    let y: Vec<f64> = x.iter().map(|v| 2.0 * v + 1.0).collect();
    let corr = stats::correlation(&x, &y).unwrap();
    let (slope, intercept) = stats::linear_regression(&x, &y).unwrap();
    println!("  y = 2x + 1: correlation = {:.6}, slope = {:.4}, intercept = {:.4}", corr, slope, intercept);
}
