/*!

# Eigenvalues decomposition

This crate contains implementations of several algorithm to
diagonalize symmetric matrices.

## Davidson Usage Example
```
// Use the Davidson method
use eigenvalues::{Davidson, DavidsonCorrection, SpectrumTarget};

// Generate random symmetric matrix
let matrix = eigenvalues::utils::generate_diagonal_dominant(20, 0.005);
let tolerance = 1e-4;

// Compute the first 2 lowest eigenvalues/eigenvectors using the DPR method
let eig = Davidson::new(
    matrix.clone(), 2, DavidsonCorrection::DPR, SpectrumTarget::Lowest, tolerance).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);

// Compute the first 2 highest eigenvalues/eigenvectors using the GJD method
let eig = Davidson::new(
    matrix, 2, DavidsonCorrection::GJD, SpectrumTarget::Highest, tolerance).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);
```

## Lanczos Usage Example
```
use nalgebra as na;

use eigenvalues::algorithms::lanczos::HermitianLanczos;
use eigenvalues::utils::{generate_random_sparse_symmetric, sort_eigenpairs};
use eigenvalues::SpectrumTarget;

// Generate sparse matrix
let matrix = generate_random_sparse_symmetric(100, 5, 0.5);

// Use 20 iterations to approximate the highest part of the spectrum
let spectrum_target = SpectrumTarget::Highest;
let lanczos = HermitianLanczos::new(matrix.clone(), 20, spectrum_target).unwrap();
let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(matrix.clone()), false);
// Compare against reference
println!("Computed eigenvalues:\n{}", lanczos.eigenvalues.rows(0, 3));
println!("Expected eigenvalues:\n{}", eig.eigenvalues.rows(0, 3));
```

*/

pub mod algorithms;
pub mod matrix_operations;
pub mod modified_gram_schmidt;
pub mod utils;
pub use algorithms::{davidson::Davidson, lanczos, DavidsonCorrection, SpectrumTarget};
pub use modified_gram_schmidt::MGS;
