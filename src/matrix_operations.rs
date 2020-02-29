/*!

## Common matrix operations for all the matrix representations.

### Other matrix representations
Currently the algorithms are implemented for the `nalgebra` **DMatrix** type.
You can use the algorithms for other matrix representations (e.g. matrix-free)
by providing your own implementation of the **Matrixoperations** trait.

*/
extern crate nalgebra as na;
use na::{DMatrix, DMatrixSlice, DVector, DVectorSlice};
use std::clone::Clone;

/// Trait containing the matrix free operations
pub trait MatrixOperations: Clone  {
    fn matrix_vector_prod(&self, vs: DVectorSlice<f64>) -> DVector<f64>;
    fn matrix_matrix_prod(&self, mtx: DMatrixSlice<f64>) -> DMatrix<f64>;
    fn diagonal(&self) -> DVector<f64>;
    fn set_diagonal(&mut self, diag: &DVector<f64>);
    fn ncols(&self) -> usize;
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


