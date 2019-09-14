/*!

# Eigenvalues decomposition

This crate contains implementations of several algorithm to
diagonalize symmetric matrices.

*/

pub mod algorithms;
pub use algorithms::davidson::EigenDavidson;
mod utils;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
