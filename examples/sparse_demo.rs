//! Sparse matrices: CSR/CSC storage and conjugate-gradient solving.
//!
//! Run with: cargo run --example sparse_demo

use mathr::sparse::{conjugate_gradient, Csr, Triplet};

fn main() {
    // 1. Build a sparse matrix from triplets: 1-D Poisson stencil
    //    A = tridiag(-1, 2, -1), n = 8 interior points.
    let n = 8usize;
    let mut triplets: Vec<Triplet> = Vec::new();
    for i in 0..n {
        triplets.push((i, i, 2.0));
        if i > 0 {
            triplets.push((i, i - 1, -1.0));
        }
        if i + 1 < n {
            triplets.push((i, i + 1, -1.0));
        }
    }
    let a = Csr::from_triplets(n, n, &triplets).unwrap();
    println!("1-D Poisson matrix: {}x{}, {} nonzeros (dense would store {})",
        n, n, a.nnz(), n * n);

    // 2. Solve A x = b with conjugate gradient.
    let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let res = conjugate_gradient(&a, &b, 1e-12, 200).unwrap();
    println!("CG solve: {} iterations, residual {:.2e}", res.iterations, res.residual);
    let xs: Vec<String> = res.x.iter().map(|v| format!("{:.6}", v)).collect();
    println!("x = [{}]", xs.join(", "));

    // 3. matvec and sparse x sparse multiplication.
    let x = vec![1.0; n];
    let y = a.matvec(&x).unwrap();
    println!("\nA · [1, 1, ..., 1] = [{}, {}, ..., {}] (interior rows are 0)", y[0], y[1], y[n / 2]);

    let aa = a.multiply(&a).unwrap();
    println!("A·A: {}x{}, {} nonzeros (pentadiagonal)", aa.rows, aa.cols, aa.nnz());

    // 4. Transpose to CSC and back.
    let csc = a.transpose();
    let back = csc.to_csr();
    let d1 = a.to_dense();
    let d2 = back.to_dense();
    let max_diff: f64 = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .map(|(i, j)| (d1[(i, j)] - d2[(i, j)]).abs())
        .fold(0.0, f64::max);
    println!("\nCSC round-trip max error: {:.2e}", max_diff);

    // 5. Jacobi-preconditioned CG is invariant under diagonal scaling.
    let s: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let mut scaled: Vec<Triplet> = Vec::new();
    for i in 0..n {
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            scaled.push((i, j, s[i] * a.values[k] * s[j]));
        }
    }
    let a2 = Csr::from_triplets(n, n, &scaled).unwrap();
    let b2: Vec<f64> = b.iter().zip(&s).map(|(bi, si)| bi * si).collect();
    let plain = conjugate_gradient(&a, &b, 1e-12, 200).unwrap();
    let prec = mathr::sparse::conjugate_gradient_jacobi(&a2, &b2, 1e-12, 200).unwrap();
    println!("\nDiagonal scaling S = diag(1..{n}): plain CG on (A, b) vs Jacobi PCG on (S·A·S, S·b):");
    println!("  iterations: {} vs {}", plain.iterations, prec.iterations);

    // 6. BiCGStab handles nonsymmetric systems (convection-diffusion);
    //    ILU(0) preconditioning accelerates it vs Jacobi.
    let ih = (n + 1) as f64;
    let mut cd: Vec<Triplet> = Vec::new();
    for i in 0..n {
        cd.push((i, i, 2.0 * ih * ih));
        if i > 0 {
            cd.push((i, i - 1, -ih * ih - ih / 2.0));
        }
        if i + 1 < n {
            cd.push((i, i + 1, -ih * ih + ih / 2.0));
        }
    }
    let acd = Csr::from_triplets(n, n, &cd).unwrap();
    let jacobi_res = mathr::sparse::bicgstab(&acd, &b, 1e-12, 200).unwrap();
    let ilu = mathr::sparse::Ilu0::factorize(&acd).unwrap();
    let ilu_res = mathr::sparse::bicgstab_ilu(&acd, &b, 1e-12, 200, &ilu).unwrap();
    println!("\nBiCGStab on nonsymmetric convection-diffusion ({}x{}, {} nnz):", n, n, acd.nnz());
    println!("  Jacobi-preconditioned : {} iterations, residual {:.2e}", jacobi_res.iterations, jacobi_res.residual);
    println!("  ILU(0)-preconditioned : {} iterations, residual {:.2e}", ilu_res.iterations, ilu_res.residual);
}
