/*!

# Eigenvalues decomposition

This crate contains implementations of several algorithm to
diagonalize symmetric matrices.

*/

pub mod algorithms;
pub mod matrix_operations;
pub mod utils;
pub use algorithms::davidson;
