/*!

# Eigenvalues decomposition

This crate contains implementations of several algorithm to
diagonalize symmetric matrices.

## Usage Example
```
extern crate eigenvalues;
extern crate nalgebra as na;

// Use the Davidson method
use eigenvalues::davidson::EigenDavidson;
use na::{DMatrix, DVector};

// Generate random symmetric matrix
let brr = eigenvalues::utils::generate_diagonal_dominant(10, 0.005);

// Compute the first 2 eigenvalues/eigenvectors
let eig = EigenDavidson::new (brr, 2).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);
```
*/

pub mod algorithms;
pub mod matrix_operations;
pub mod utils;
pub use algorithms::davidson;
