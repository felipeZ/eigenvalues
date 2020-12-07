/*!

## Common matrix operations for all the matrix representations.

### Other matrix representations
Currently the algorithms are implemented for the `nalgebra` **DMatrix** type.
You can use the algorithms for other matrix representations (e.g. matrix-free)
by providing your own implementation of the **Matrixoperations** trait.

*/
use nalgebra::{DMatrix, DMatrixSlice, DVector, DVectorSlice};
use std::clone::Clone;

/// Trait containing the matrix free operations
pub trait MatrixOperations: Clone  {
    /// Matrix vector multiplication
    fn matrix_vector_prod(&self, vs: DVectorSlice<f64>) -> DVector<f64>;
    /// Matrix matrix multiplication
    fn matrix_matrix_prod(&self, mtx: DMatrixSlice<f64>) -> DMatrix<f64>;
    /// Get the matrix diagonal
    fn diagonal(&self) -> DVector<f64>;
    /// Set the matrix diagonal
    fn set_diagonal(&mut self, diag: &DVector<f64>);
    /// Get the number of columns
    fn ncols(&self) -> usize;
    /// Get the number of rows
    fn nrows(&self) -> usize;
}

impl MatrixOperations for DMatrix<f64> {
    fn matrix_vector_prod(&self, vs: DVectorSlice<f64>) -> DVector<f64> {
        self * vs
    }
    fn matrix_matrix_prod(&self, mtx: DMatrixSlice<f64>) -> DMatrix<f64> {
        self * mtx
    }
    fn diagonal(&self) -> DVector<f64> {
        self.diagonal()
    }
    fn set_diagonal(&mut self, diag: &DVector<f64>) {
        self.set_diagonal(diag);
    }

    fn ncols(&self) -> usize {
        self.ncols()
    }
    fn nrows(&self) -> usize {
        self.nrows()
    }
}


