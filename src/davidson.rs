//! ### Davidson implementation

extern crate nalgebra as na;
use na::linalg::SymmetricEigen;
use na::{DMatrix, DVector, Dynamic};

/// Structure with the configuration data
pub struct EigenDavidson {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64>,
}

impl EigenDavidson {
    /// The new static method takes the following arguments:
    /// * `h` - A highly diagonal symmetric matrix
    /// * `nvalues` - the number of eigenvalues/eigenvectors pair to compute

    pub fn new(h: DMatrix<f64>, nvalues: usize)  -> Result<EigenDavidson, &'static str> {

    // Maximum number of iterations
    let max_iters = 100;
    
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
    let mut result = Err("Algorithm didn't converge!");
    for i in 0..max_iters {
        // 2. Generate subpace matrix problem by projecting into the basis
        let subspace = basis.columns(0, dim_sub);
        let matrix_proj = subspace.transpose() * (&h * subspace);

        // 3. compute the eigenvalues and their corresponding ritz_vectors
        let eig = SymmetricEigen::new(matrix_proj);

        // 4. Check for convergence
        // 4.1 Compute the residues
        let ritz_vectors = subspace * eig.eigenvectors.columns(0, dim_sub);
        let residues = compute_residues(&h, &eig.eigenvalues, &ritz_vectors);

        // 4.2 Check Converge for each pair eigenvalue/eigenvector
        let errors =
            DVector::<f64>::from_iterator(nvalues, residues.column_iter().map(|col| col.norm()));

        // 4.3 Check if all eigenvalues/eigenvectors have converged
        if errors.iter().all(|&x| x < TOLERANCE) {
            result = Ok(create_results(&eig.eigenvalues, &ritz_vectors, nvalues));
            break;
        }
        // 5. Update subspace basis set
        // 5.1 Add the correction vectors to the current basis
        if dim_sub <= max_dim_sub {
            let correction = compute_correction(&h, residues, eig, &method);
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
    result
    }
}

/// Extract the requested eigenvalues/eigenvectors pairs
fn create_results(
    subspace_eigenvalues: &DVector<f64>,
    ritz_vectors: &DMatrix<f64>,
    nvalues: usize,
) -> EigenDavidson {
    let eigenvectors = DMatrix::<f64>::from_iterator(
        ritz_vectors.nrows(),
        nvalues,
        subspace_eigenvalues.columns(0, nvalues).iter().cloned(),
    );
    let eigenvalues =
        DVector::<f64>::from_iterator(nvalues, ritz_vectors.rows(0, nvalues).iter().cloned());
    EigenDavidson {
        eigenvalues,
        eigenvectors,
    }
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
    residues: DMatrix<f64>,
    eigenpairs: SymmetricEigen<f64, Dynamic>,
    method: &str,
) -> DMatrix<f64> {
    match method.to_uppercase().as_ref() {
        "DPR" => compute_dpr_correction(&h, residues, &eigenpairs.eigenvalues),
        "GJD" => compute_gjd_correction(&h, residues, &eigenpairs),
        _ => panic!("Method {} has not been implemented", method),
    }
}

/// Use the Diagonal-Preconditioned-Residue (DPR) method to compute the correction
fn compute_dpr_correction(
    h: &DMatrix<f64>,
    residues: DMatrix<f64>,
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

/// Use the Generalized Jacobi Davidson (GJD) to compute the correction
fn compute_gjd_correction(
    h: &DMatrix<f64>,
    residues: DMatrix<f64>,
    eigenpairs: &SymmetricEigen<f64, Dynamic>,
) -> DMatrix<f64> {
    let dimx = h.nrows();
    let dimy = eigenpairs.eigenvalues.nrows();
    let id = DMatrix::<f64>::identity(dimx, dimx);
    let ones = DVector::<f64>::repeat(dimx, 1.0);
    let mut correction = DMatrix::<f64>::zeros(dimx, dimy);
    for (k, r) in eigenpairs.eigenvalues.iter().enumerate() {
        // Create the components of the linear system
        let x = eigenpairs.eigenvectors.column(k);
        let t1 = &id - x * x.transpose();
        let mut t2 = h.clone();
        t2.set_diagonal(&(*r * &ones));
        let arr = &t1 * t2 * &t1;
        // Solve the linear system
        let decomp = arr.lu(); 
        let mut b = -residues.column(k);
        decomp.solve_mut(&mut b);
        correction.set_column(k, &b);
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
    extern crate nalgebra as na;
    use na::{DMatrix, DVector};

    #[test]
    fn test_update_subspace() {
        let mut arr = DMatrix::<f64>::repeat(3, 3, 1.);
        let brr = DMatrix::<f64>::zeros(3, 2);
        super::update_subspace(&mut arr, brr, 0, 2);
        assert_eq!(arr.column(1).sum(), 0.);
        assert_eq!(arr.column(2).sum(), 3.);
    }

}
