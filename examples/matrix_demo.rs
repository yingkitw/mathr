use mathr::matrix::Matrix;

fn main() {
    // --- Create matrices ---
    let a = Matrix::from_rows(&[
        vec![1.0, 2.0, 3.0],
        vec![4.0, 5.0, 6.0],
        vec![7.0, 8.0, 10.0],
    ]).unwrap();
    println!("Matrix A:");
    print_matrix(&a);

    // --- Determinant ---
    let det = a.determinant().unwrap();
    println!("det(A) = {:.6}\n", det);

    // --- Transpose ---
    let at = a.transpose();
    println!("Aᵀ:");
    print_matrix(&at);

    // --- Matrix multiplication ---
    let b = Matrix::from_rows(&[
        vec![1.0, 0.0],
        vec![0.0, 1.0],
        vec![1.0, 1.0],
    ]).unwrap();
    let c = (&a * &b).unwrap();
    println!("A × B (3×3 × 3×2):");
    print_matrix(&c);

    // --- Inverse ---
    let inv = a.inverse().unwrap();
    println!("A⁻¹:");
    print_matrix(&inv);

    // --- Verify A × A⁻¹ = I ---
    let identity = (&a * &inv).unwrap();
    println!("A × A⁻¹ (should be identity):");
    print_matrix(&identity);

    // --- Solve linear system Ax = b ---
    let b_vec = vec![6.0, 15.0, 25.0];
    let x = a.solve(&b_vec).unwrap();
    println!("Solve Ax = {:?}:  x = {:?}", b_vec, x);

    // Verify
    let check = a.mul_vec(&x).unwrap();
    println!("  Verification: Ax = {:?}", check);

    // --- Trace ---
    let tr = a.trace().unwrap();
    println!("\ntrace(A) = {:.1}", tr);
}

fn print_matrix(m: &Matrix) {
    for i in 0..m.rows {
        let row: Vec<String> = (0..m.cols).map(|j| format!("{:8.4}", m.get(i, j))).collect();
        println!("  [{}]", row.join(" "));
    }
    println!();
}
