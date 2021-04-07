/*!

## Algorithms to compute (some) eigenvalues/eigenvectors for symmetric matrices.

*/
pub mod davidson;
pub mod lanczos;

/// Option to compute the lowest, highest or somewhere in the middle part of the
/// spectrum
#[derive(Clone, PartialEq)]
pub enum SpectrumTarget {
    Lowest,
    Highest,
    Target(f64),
}

/// Correction method for the Davidson algorithm
#[derive(Debug, Copy, Clone)]
pub enum DavidsonCorrection {
    DPR,
    GJD,
}
