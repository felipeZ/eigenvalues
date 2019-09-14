extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// Generate a random highly diagonal symmetric matrix
fn generate_diagonal_dominant(dim: usize, sparsity: f64) -> DMatrix<f64> {
    let xs = 1..(dim + 1);
    let it = xs.map(|x: usize| x as f64);
    let vs = DVector::<f64>::from_iterator(dim, it);
    let mut arr = DMatrix::<f64>::new_random(dim, dim);
    arr += &arr.transpose();
    arr *= sparsity;
    arr.set_diagonal(&vs);
    arr
}
