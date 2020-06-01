extern crate approx;
extern crate eigenvalues;
extern crate nalgebra as na;

use eigenvalues::algorithms::lanczos::HermitianLanczos;
use eigenvalues::SpectrumTarget;
use eigenvalues::utils::generate_diagonal_dominant;
use eigenvalues::utils::{sort_eigenpairs, test_eigenpairs};


#[test]
fn test_lanczos() {
    let matrix = generate_diagonal_dominant(10, 1.0);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(matrix.clone()), false);
    let spectrum_target = SpectrumTarget::Highest;

    let lanczos = HermitianLanczos::new(matrix.clone(), 1, spectrum_target).unwrap();

    println!("Computed eigenvalues:\n{}", lanczos.eigenvalues[1]);
    println!("Expected eigenvalues:\n{}", eig.eigenvalues[1]);
    test_eigenpairs(&eig, (lanczos.eigenvalues, lanczos.eigenvectors), 2);

}
