use mathr::complex::Complex;
use mathr::fft;

fn main() {
    // --- FFT of a pure sinusoid ---
    // 16 samples of a 2 Hz sine wave sampled at 16 Hz
    let n = 16;
    let freq = 2.0;
    let samples: Vec<f64> = (0..n)
        .map(|i| (2.0 * std::f64::consts::PI * freq * i as f64 / n as f64).sin())
        .collect();

    let mags = fft::magnitude_spectrum(&samples).unwrap();
    println!("Magnitude spectrum (16-point FFT of 2 Hz sine):");
    for (k, m) in mags.iter().enumerate() {
        if *m > 1e-6 {
            println!("  bin {:2}: {:.4}", k, m);
        }
    }

    // --- Convolution via FFT ---
    let a = vec![1.0, 2.0, 3.0, 4.0];
    let b = vec![1.0, 0.0, -1.0];
    let c = fft::convolve(&a, &b).unwrap();
    println!("\nConvolution of [1,2,3,4] * [1,0,-1]:");
    println!("  {:?} * {:?} = {:?}", a, b, c);

    // --- Cross-correlation ---
    let x = vec![1.0, 2.0, 3.0, 4.0];
    let corr = fft::cross_correlate(&x, &x).unwrap();
    println!("\nAuto-correlation of [1,2,3,4]:");
    println!("  {:?}", corr);
    let peak_idx = corr
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    println!("  peak at index {} (value {:.1})", peak_idx, corr[peak_idx]);

    // --- Window functions ---
    let signal = vec![1.0; 8];
    let windowed = fft::apply_window(&signal, fft::Window::Hann);
    println!("\nHann window applied to 8 ones:");
    println!("  {:?}", windowed);

    // --- Inverse FFT round-trip ---
    let input: Vec<Complex<f64>> = (0..8)
        .map(|i| Complex::new((i as f64).sin(), (i as f64).cos()))
        .collect();
    let freq_data = fft::fft(&input).unwrap();
    let recovered = fft::ifft(&freq_data).unwrap();
    let max_err = input
        .iter()
        .zip(recovered.iter())
        .map(|(a, b)| (a.re - b.re).abs() + (a.im - b.im).abs())
        .fold(0.0, f64::max);
    println!("\nFFT → IFFT round-trip max error: {:.2e}", max_err);
}
