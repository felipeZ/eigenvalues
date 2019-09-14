extern crate approx;
extern crate eigenvalues;
extern crate nalgebra as na;

use approx::relative_eq;
use eigenvalues::algorithms::davidson::EigenDavidson;
use eigenvalues::utils::generate_diagonal_dominant;
use eigenvalues::utils::sort_eigenpairs;

#[test]
fn davidson_eigenvalues() {
    let arr = generate_diagonal_dominant(10, 0.005);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(arr.clone()));
    let dav_eig = EigenDavidson::new(arr, 2).unwrap();
    for i in 0..2 {
        // Test Eigenvalues
        relative_eq!(eig.eigenvalues[i], dav_eig.eigenvalues[i]);
        // Test Eigenvectors
        let x = eig.eigenvectors.column(i);
        let y = dav_eig.eigenvectors.column(i);
        // The autovectors may different in their sign
        // They should be either parallel or antiparallel
        let dot = x.dot(&y).abs();
        relative_eq!(dot, 1.0);
    }
}
