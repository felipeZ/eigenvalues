extern crate nalgebra as na;
use na::linalg::SymmetricEigen;
use na::Dynamic;
use na::{DMatrix, DVector};

/// Generate a random highly diagonal symmetric matrix
pub fn generate_diagonal_dominant(dim: usize, sparsity: f64) -> DMatrix<f64> {
    let xs = 1..(dim + 1);
    let it = xs.map(|x: usize| x as f64);
    let vs = DVector::<f64>::from_iterator(dim, it);
    let mut arr = DMatrix::<f64>::new_random(dim, dim);
    arr += &arr.transpose();
    arr *= sparsity;
    arr.set_diagonal(&vs);
    arr
}

/// Sort the eigenvalues and their corresponding eigenvectors in ascending order
pub fn sort_eigenpairs(eig: SymmetricEigen<f64, Dynamic>) -> SymmetricEigen<f64, Dynamic> {
    // Sort the eigenvalues
    let mut vs: Vec<(f64, usize)> = eig
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(idx, &x)| (x, idx))
        .collect();
    vs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    // Sorted eigenvalues
    let eigenvalues = DVector::<f64>::from_iterator(vs.len(), vs.iter().map(|t| t.0));

    // Indices of the sorted eigenvalues
    let indices: Vec<_> = vs.iter().map(|t| t.1).collect();

    // Create sorted eigenvectors
    let dim_rows = eig.eigenvectors.nrows();
    let dim_cols = eig.eigenvectors.ncols();
    let mut eigenvectors = DMatrix::<f64>::zeros(dim_rows, dim_cols);

    for i in 0..dim_cols {
        eigenvectors.set_column(i, &eig.eigenvectors.column(indices[i]));
    }
    SymmetricEigen {
        eigenvalues,
        eigenvectors,
    }
}
