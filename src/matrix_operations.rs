/*!

## Common matrix operations for all the matrix representations.


*/
extern crate nalgebra as na;
use na::{DMatrix, DMatrixSlice, DVector, DVectorSlice};

/// Trait containing the matrix free operations
pub trait MatrixOperations {
    fn matrix_vector_prod(&self, vs: DVectorSlice<f64>) -> DVector<f64>;
    fn matrix_matrix_prod(&self, mtx: DMatrixSlice<f64>) -> DMatrix<f64>;
    fn diagonal(&self) -> DVector<f64>;
    fn set_diagonal(&mut self, diag: &DVector<f64>);
    fn cols(&self) -> usize;
    fn rows(&self) -> usize;
    fn clone(&self) -> Self;
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

    fn cols(&self) -> usize {
        self.ncols()
    }
    fn rows(&self) -> usize {
        self.nrows()
    }
    fn clone(&self) -> DMatrix<f64> {
        std::clone::Clone::clone(&self)
    }
}
