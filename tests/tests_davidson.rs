extern crate nalgebra as na;

use eigenvalues::algorithms::davidson::Davidson;
use eigenvalues::utils::generate_diagonal_dominant;
use eigenvalues::utils::{sort_eigenpairs, test_eigenpairs};
use eigenvalues::{DavidsonCorrection, SpectrumTarget};

#[test]
fn test_davidson_lowest() {
    let arr = generate_diagonal_dominant(10, 0.005);
    let eig = sort_eigenpairs(nalgebra::linalg::SymmetricEigen::new(arr.clone()), true);
    let spectrum_target = SpectrumTarget::Lowest;
    let tolerance = 1.0e-4;

    let dav = Davidson::new(
        arr.clone(),
        2,
        DavidsonCorrection::DPR,
        spectrum_target.clone(),
        tolerance,
    )
    .unwrap();
    test_eigenpairs(&eig, (dav.eigenvalues, dav.eigenvectors), 2);
    let dav = Davidson::new(
        arr.clone(),
        2,
        DavidsonCorrection::GJD,
        spectrum_target,
        tolerance,
    )
    .unwrap();
    test_eigenpairs(&eig, (dav.eigenvalues, dav.eigenvectors), 2);
}

#[test]
fn test_davidson_unsorted() {
    // Test the algorithm when the diagonal is unsorted
    let mut arr = generate_diagonal_dominant(8, 0.005);
    let tolerance = 1.0e-6;
    let vs = nalgebra::DVector::<f64>::from_vec(vec![3.0, 2.0, 4.0, 1.0, 5.0, 6.0, 7.0, 8.0]);
    arr.set_diagonal(&vs);
    let eig = sort_eigenpairs(nalgebra::linalg::SymmetricEigen::new(arr.clone()), true);
    let dav = Davidson::new(
        arr,
        1,
        DavidsonCorrection::DPR,
        SpectrumTarget::Lowest,
        tolerance,
    )
    .unwrap();
    test_eigenpairs(&eig, (dav.eigenvalues, dav.eigenvectors), 1);
}

#[test]
fn test_davidson_highest() {
    // Test the compution of the highest eigenvalues
    let dim = 20;
    let nvalues = 2;
    let tolerance = 1.0e-4;
    let arr = generate_diagonal_dominant(dim, 0.005);
    let eig = sort_eigenpairs(nalgebra::linalg::SymmetricEigen::new(arr.clone()), false);

    let target = SpectrumTarget::Highest;
    println!("running DPR");
    let dav = Davidson::new(
        arr.clone(),
        nvalues,
        DavidsonCorrection::DPR,
        target.clone(),
        tolerance,
    )
    .unwrap();
    test_eigenpairs(&eig, (dav.eigenvalues, dav.eigenvectors), nvalues);
    println!("running GJD");
    let dav = Davidson::new(
        arr.clone(),
        nvalues,
        DavidsonCorrection::GJD,
        target,
        tolerance,
    )
    .unwrap();
    test_eigenpairs(&eig, (dav.eigenvalues, dav.eigenvectors), nvalues);
}
