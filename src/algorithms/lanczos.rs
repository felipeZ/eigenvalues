/*!

# Hermitian Lanczos algorithm

The [Hermitian Lanczos](https://en.wikipedia.org/wiki/Lanczos_algorithm) is an algorithm to compute the lowest/highest
eigenvalues of an hermitian matrix using a [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace)

*/
extern crate nalgebra as na;
use super::SpectrumTarget;
use crate::matrix_operations::MatrixOperations;
use na::{DMatrix, DVector};

pub struct HermitianLanczos {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64>,
}

impl HermitianLanczos {
    /// The new static method takes the following arguments:
    /// * `h` - A highly diagonal symmetric matrix
    /// * `nvalues` - the number of eigenvalues/eigenvectors pair to compute
    /// * `spectrum_target` Lowest or Highest part of the spectrum
    /// * `tolerance` numerical tolerance.

    pub fn new<M: MatrixOperations>(
        h: M,
        nvalues: usize,
        spectrum_target: SpectrumTarget,
        tolerance: f64,
    ) -> Result<Self, &'static str> {
        let eigenvalues = DVector::<f64>::zeros(h.nrows());
        let eigenvectors = DMatrix::<f64>::zeros(h.nrows(), h.ncols());
        let max_iters = 50;

        // Off-diagonal elements
        let mut betas = DVector::<f64>::zeros(max_iters);
        // Diagonal elements
        let mut alphas = DVector::<f64>::zeros(max_iters);

        // Initial guess of the eigenvector
        let mut vs = DVector::<f64>::zeros(h.ncols());
	vs[1] = 0.0;

        let mut residues = DVector::<f64>::zeros(h.ncols());
        // Compute the elements of the tridiagonal matrix
        for i in 0..max_iters {
            let tmp = &h.matrix_vector_prod(&vs);
            alphas[i] = tmp.dot(&vs);
	    vs.copy_from(&(tmp));
            betas[i] = residues.norm();
            if (betas[i]) < tolerance {
                break;
            } else {
                vs.copy_from(&residues / betas[i]);
            }
        }

        Ok(HermitianLanczos {
            eigenvalues,
            eigenvectors,
        })
    }
}
