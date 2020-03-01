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
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(matrix.clone()), true);
    let spectrum_target = SpectrumTarget::Lowest;
    let tolerance = 1.0e-4;

    let lanczos = HermitianLanczos::new(matrix.clone(), 2, spectrum_target, tolerance).unwrap();

    test_eigenpairs(&eig, (lanczos.eigenvalues, lanczos.eigenvectors), 2);

}
