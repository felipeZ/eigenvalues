/*!

# Hermitian Lanczos algorithm

The [Hermitian Lanczos](https://en.wikipedia.org/wiki/Lanczos_algorithm) is an algorithm to compute the lowest/highest
eigenvalues of an hermitian matrix using a [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace)

*/
extern crate nalgebra as na;
use super::SpectrumTarget;
use crate::matrix_operations::MatrixOperations;
use crate::utils;
use na::linalg::SymmetricEigen;
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

    pub fn new<M: MatrixOperations>(
        h: M,
        nvalues: usize,
        spectrum_target: SpectrumTarget,
    ) -> Result<Self, &'static str> {
        let max_iters = (nvalues as f64 * 3.0).floor() as usize;

        // Off-diagonal elements
        let mut betas = DVector::<f64>::zeros(max_iters - 1);
        // Diagonal elements
        let mut alphas: DVector<f64> = DVector::<f64>::zeros(max_iters);

        // Matrix with the orthognal vectors
        let mut vs = DMatrix::<f64>::zeros(h.nrows(), max_iters);

        // Initial vector
        let xs = DVector::<f64>::new_random(h.nrows()).normalize();
        vs.set_column(0, &xs);

        // Compute the elements of the tridiagonal matrix
        for i in 0..max_iters {
            let tmp: DVector<f64> = h.matrix_vector_prod(vs.column(i));
            alphas[i] = tmp.dot(&vs.column(i));
            let tmp = {
                if i == 0 {
                    &tmp - alphas[0] * vs.column(0)
                } else {
                    &tmp - alphas[i] * vs.column(i) - betas[i - 1] * vs.column(i - 1)
                }
            };
            if i < max_iters - 1 {
                betas[i] = tmp.norm();
                if betas[i] > 1.0e-8 {
                    vs.set_column(i + 1, &(tmp / betas[i]));
                } else {
                    vs.set_column(i + 1, &tmp);
                }
            }
        }

        let tridiagonal = Self::construct_tridiagonal(&alphas, &betas);
        let ord_sort = match spectrum_target {
            SpectrumTarget::Highest => false,
            _ => true,
        };
        let eig = utils::sort_eigenpairs(SymmetricEigen::new(tridiagonal), ord_sort);
        let eigenvalues = eig.eigenvalues;
        let eigenvectors = vs * eig.eigenvectors; // Ritz vectors

        Ok(HermitianLanczos {
            eigenvalues,
            eigenvectors,
        })
    }

    fn construct_tridiagonal(alphas: &DVector<f64>, betas: &DVector<f64>) -> DMatrix<f64> {
        let dim = alphas.len();
        let lambda = |i, j| {
            if i == j {
                alphas[i]
            } else if i == j + 1 {
                betas[j]
            } else if j == i + 1 {
                betas[i]
            } else {
                0.0
            }
        };
        DMatrix::<f64>::from_fn(dim, dim, lambda)
    }
}
