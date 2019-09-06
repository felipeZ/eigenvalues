//! ### Davidson implementation

extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// Structure with the configuration data
pub struct Eigen {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64>,
}

/// API
pub fn davidson(h: DMatrix<f64>, nvalues: usize, max_iters: usize) -> Result<Eigen, &'static str> {
    // Numerical tolerance of the algorithm
    const TOLERANCE: f64 = 1e-8;

    // Method to compute the correction
    let method = "DPR";

    // Data used for the algorithm
    let max_dim_sub = nvalues * 10;

    // Initial subpace
    let mut dim_sub = nvalues * 2;
    // 1. Select the initial ortogonal subspace based on lowest elements
    let mut basis = generate_subspace(&h.diagonal(), max_dim_sub);

    // Outer loop block Davidson schema
    let mut i = 0;
    loop {
        // 2. Generate subpace matrix problem by projecting into V
        let subspace = basis.columns(0, dim_sub);
        let matrix_proj = subspace.transpose() * (&h * subspace);

        // 3. compute the eigenvalues and their corresponding ritz_vectors
        let eig = na::linalg::SymmetricEigen::new(matrix_proj);

        // 4. Check for convergence
        // 4.1 Compute the residues
        let ritz_vectors = subspace * eig.eigenvectors.columns(0, dim_sub);
        let residues = compute_residues(&h, &eig.eigenvalues, &ritz_vectors);

        // 4.2 Check Converge for each pair eigenvalue/eigenvector
        let mut errors =
            DVector::<f64>::from_iterator(nvalues, residues.column_iter().map(|col| col.norm()));

        // 4.3 Check if all eigenvalues/eigenvectors have converged
        if errors.iter().all(|&x| x < TOLERANCE) {
            break;
        }
        // 5. Update subspace basis set
        // 5.1 Add the correction vectors to the current basis
        if dim_sub <= max_dim_sub {
            let correction = compute_correction(&h, &residues, &eig.eigenvalues, &method);
            update_subspace(&mut basis, correction, dim_sub, dim_sub * 2);

            // 6. Orthogonalize the subspace
            basis = orthogonalize_subspace(basis);
            // update counter
            dim_sub *= 2;

        // 5.2 Otherwise reduce the basis of the subspace to the current correction
        } else {
            dim_sub = nvalues * 2;
            basis.fill(0.0);
            update_subspace(&mut basis, ritz_vectors, 0, dim_sub);
        }
        // Check number of iterations
        if i > max_iters {
            break;
        }
    }
    let eigenvalues = DVector::<f64>::zeros(3);
    let eigenvectors = DMatrix::<f64>::identity(3, 3);
    Ok(Eigen {
        eigenvalues,
        eigenvectors,
    })
}

/// Update the subpace with new vectors
fn update_subspace(basis: &mut DMatrix<f64>, vectors: DMatrix<f64>, start: usize, end: usize) {
    for k in start..end {
        basis.set_column(k, &vectors.column(k));
    }
}

/// Orthogonalize the subpsace using the QR method
fn orthogonalize_subspace(basis: DMatrix<f64>) -> DMatrix<f64> {
    let qr = na::linalg::QR::new(basis);
    qr.q()
}

/// Residue vectors
fn compute_residues(
    h: &DMatrix<f64>,
    eigenvalues: &DVector<f64>,
    ritz_vectors: &DMatrix<f64>,
) -> DMatrix<f64> {
    let dim_sub = eigenvalues.nrows();
    let mut residues = DMatrix::<f64>::zeros(h.nrows(), dim_sub);
    for k in 0..dim_sub {
        let guess = eigenvalues[k] * ritz_vectors.column(k);
        let vs = h * ritz_vectors.column(k);
        residues.set_column(k, &(vs - guess));
    }
    residues
}

/// compute the correction vectors using either DPR or GJD
fn compute_correction(
    h: &DMatrix<f64>,
    residues: &DMatrix<f64>,
    eigenvalues: &DVector<f64>,
    method: &str,
) -> DMatrix<f64> {
    match method.to_uppercase().as_ref() {
        "DPR" => compute_dpr_correction(&h, &residues, &eigenvalues),
        _ => panic!("Method {} has not been implemented", method),
    }
}

/// Use the Diagonal-Preconditioned-Residue (DPR) method to compute the correction
fn compute_dpr_correction(
    h: &DMatrix<f64>,
    residues: &DMatrix<f64>,
    eigenvalues: &DVector<f64>,
) -> DMatrix<f64> {
    let d = h.diagonal();
    let mut correction = DMatrix::<f64>::zeros(h.nrows(), residues.ncols());
    for (k, r) in eigenvalues.iter().enumerate() {
        let rs = DVector::<f64>::repeat(h.nrows(), *r);
        let x = (rs - &d) * residues.column(k);
        correction.set_column(k, &x);
    }
    correction
}

/// Generate initial orthonormal subspace
fn generate_subspace(diag: &DVector<f64>, max_dim_sub: usize) -> DMatrix<f64> {
    // TODO implement the case when the diagonal is not sorted
    DMatrix::<f64>::identity(diag.nrows(), max_dim_sub)
}

#[cfg(test)]
mod test {

    #[test]
    fn test_something() {
        let universe = 42;
        assert_eq!(universe, 42);
    }
}
