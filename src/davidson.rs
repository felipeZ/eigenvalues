//! ### Davidson implementation

extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// Structure with the configuration data
pub struct Eigen {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64> 
}

/// API
pub fn davidson(
    h: DMatrix<f64>,
    nvalues: usize,
    max__iters: usize,
) -> Eigen {
    // Numerical tolerance of the algorithm
    const TOLERANCE: f64 = 1e-8;

    // Data used for the algorithm
    let max_dim_sub = nvalues * 10;

    // Initial subpace
    let mut dim_sub = nvalues * 2;
    // 1. Select the initial ortogonal subspace based on lowest elements
    let mut V = generate_subspace(&h.diagonal(), max_dim_sub);

    // Outer loop block Davidson schema
    // for i in 0..max__iters {
    for i in 0..1 {
        // 2. Generate subpace matrix problem by projecting into V
        let subspace = V.columns(0, dim_sub);
        let matrix_proj = subspace.transpose() * (&h * subspace);

        // 3. compute the eigenvalues and their corresponding ritz_vectors
        let eig = na::linalg::SymmetricEigen::new(matrix_proj);

        // 4. Check for convergence
        let ritz_vectors = subspace * eig.eigenvectors.columns(0, nvalues);

        let mut errors = DVector::<f64>::zeros(nvalues);
        let mut residues = DMatrix::<f64>::zeros(h.ncols(), nvalues);

        // 4.1 Check Converge for each pair eigenvalue/eigenvector
        for k in 0..nvalues {
            let guess = eig.eigenvalues[k] * ritz_vectors.column(k);
            let vs = &h * ritz_vectors.column(k);
            residues.set_column(k, &(vs - guess));
            let err = residues.column(k).norm();
            errors[k] = err;
        }

        // 4.2 Check if all eigenvalues/eigenvectors have converged
        if errors.iter().all(|&x| x < TOLERANCE) {
            break;
        }
        // 5. Add the correction vectors to the current basis
        // if dim_sub <= max_dim {
        //     correction = compute_correction(&h, &rs);

        //     dim_sub *= 2;
        // } else {
        //     dim_sub = nvalues * 2;
        // }
    }
    let eigenvalues = DVector::<f64>::zeros(3);
    let eigenvectors = DMatrix::<f64>::identity(3, 3);
    Eigen{eigenvalues, eigenvectors}
}

/// Generate initial orthonormal subspace
pub fn generate_subspace(diag: &DVector<f64>, max_dim_sub: usize) -> DMatrix<f64> {
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
