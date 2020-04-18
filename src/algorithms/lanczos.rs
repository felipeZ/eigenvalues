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
        let eigenvectors = DMatrix::<f64>::zeros(h.nrows(), nvalues);
        let max_iters = (nvalues as f64 * 1.5).floor() as usize;

        // Off-diagonal elements
        let mut betas = DVector::<f64>::zeros(max_iters);
        // Diagonal elements
        let mut alphas: DVector<f64> = h.diagonal().clone();

        // Matrix with the orthognal vectors
        let mut vs = DMatrix::<f64>::zeros(h.nrows(), max_iters);

        // Initial vector
        let xs = DVector::<f64>::new_random(h.nrows()).normalize();
        vs.set_column(1, &xs);

        // let mut residues = DVector::<f64>::zeros(h.nrows());
        // Compute the elements of the tridiagonal matrix
        for i in 1..max_iters {
            let tmp = &h.matrix_vector_prod(vs.column(i)) * (betas[i] * vs[i - 1]);
            alphas[i] = tmp.dot(&vs.column(i));
            let tmp = tmp - alphas[i] * vs.column(i);
	    if i < max_iters - 1 {
	        betas[i + 1] = tmp.norm_squared();
		vs.set_column(i + 1, &(tmp / betas[i + 1]));
	    }
        }

        Ok(HermitianLanczos {
            eigenvalues,
            eigenvectors,
        })
    }
}
