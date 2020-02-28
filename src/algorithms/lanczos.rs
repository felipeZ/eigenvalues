/*!

# Hermitian Lanczos algorithm

The [Hermitian Lanczos](https://en.wikipedia.org/wiki/Lanczos_algorithm) is an algorithm to compute the lowest/highest
eigenvalues of an hermitian matrix using a [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace)

*/
extern crate nalgebra as na;

Pub struct HermitianLanczos {
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
    }
}
