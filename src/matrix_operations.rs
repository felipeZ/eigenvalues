/*!

## Common matrix operations for all the matrix representations.


*/
extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// Trait containing the matrix free operations
pub trait MatrixOperations {
    fn matrix_vec_prod(&self, vs: DVector<f64>) -> DVector<f64>;
    fn matrix_matrix_prod(&self, mtx: DMatrix<f64>) -> DMatrix<f64>;
    fn diagonal(&self) -> DVector<f64>;
}

impl MatrixOperations for DMatrix<f64> {
    fn matrix_vec_prod(&self, vs: DVector<f64>) -> DVector<f64> {
        self * vs
    }
    fn matrix_matrix_prod(&self, mtx: DMatrix<f64>) -> DMatrix<f64> {
        self * mtx
    }
    fn diagonal(&self) -> DVector<f64> {
        self.diagonal()
    }
}
