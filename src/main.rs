use nalgebra as na;

use eigenvalues::algorithms::lanczos::HermitianLanczos;
use eigenvalues::utils::{generate_random_sparse_symmetric, sort_eigenpairs};
use eigenvalues::SpectrumTarget;

fn main() {
    let matrix = generate_random_sparse_symmetric(200, 5, 0.5);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(matrix.clone()), false);
    let spectrum_target = SpectrumTarget::Highest;
    let lanczos = HermitianLanczos::new(matrix.clone(), 50, spectrum_target).unwrap();
    println!("Computed eigenvalues:\n{}", lanczos.eigenvalues.rows(0, 3));
    println!("Expected eigenvalues:\n{}", eig.eigenvalues.rows(0, 3));
    let x = eig.eigenvectors.column(0);
    let y = lanczos.eigenvectors.column(0);
    println!("parallel:{}", x.dot(&y));
}
