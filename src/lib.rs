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

// Compute the first 2 eigenvalues/eigenvectors using the DPR method
let eig = EigenDavidson::new (brr, 2, "DPR", None).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);
```
## Highest eigenvalues
By default the library compute the Lowest eigenvalues of the 
spectrum. It is also, possible to compute the highest eigenvalues like:
```
extern crate eigenvalues;
extern crate nalgebra as na;

// Use the Davidson method
use eigenvalues::SpectrumTarget;
use eigenvalues::davidson::EigenDavidson;
use na::{DMatrix, DVector};

// Generate random symmetric matrix
let brr = eigenvalues::utils::generate_diagonal_dominant(10, 0.005);

// Compute the first 2 eigenvalues/eigenvectors using the DPR method
let eig = EigenDavidson::new (brr, 2, "DPR", Some(SpectrumTarget::Highest)).unwrap();
println!("eigenvalues:{}", eig.eigenvalues);
println!("eigenvectors:{}", eig.eigenvectors);
```


*/

pub mod algorithms;
pub mod matrix_operations;
pub mod utils;
pub use algorithms::{SpectrumTarget, davidson};
