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
use eigenvalues::SpectrumTarget;
use na::{DMatrix, DVector};

// Generate random symmetric matrix
let brr = eigenvalues::utils::generate_diagonal_dominant(20, 0.005);
let tolerance = 1e-6;

// Compute the first 2 lowest eigenvalues/eigenvectors using the DPR method
let eig = EigenDavidson::new (brr.clone(), 2, "DPR", SpectrumTarget::Lowest, tolerance).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);

// Compute the first 2 highest eigenvalues/eigenvectors using the GJD method
let eig = EigenDavidson::new (brr, 2, "GJD", SpectrumTarget::Highest, tolerance).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);
```

*/

pub mod algorithms;
pub mod matrix_operations;
pub mod modified_gram_schmidt;
pub mod utils;
pub use algorithms::{davidson, SpectrumTarget};
pub use modified_gram_schmidt::MGS;
