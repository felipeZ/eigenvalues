extern crate approx;
extern crate eigenvalues;
extern crate nalgebra as na;

use eigenvalues::algorithms::lanczos::HermitianLanczos;
use eigenvalues::utils::{generate_random_sparse_symmetric, sort_eigenpairs, test_eigenpairs};
use eigenvalues::SpectrumTarget;

#[test]
fn test_lanczos() {
    let matrix = generate_random_sparse_symmetric(100, 5, 0.5);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(matrix.clone()), false);
    let spectrum_target = SpectrumTarget::Highest;
    let lanczos = HermitianLanczos::new(matrix.clone(), 40, spectrum_target).unwrap();

    println!("Computed eigenvalues:\n{}", lanczos.eigenvalues[0]);
    println!("Expected eigenvalues:\n{}", eig.eigenvalues[0]);
    test_eigenpairs(&eig, (lanczos.eigenvalues, lanczos.eigenvectors), 1);
}
