extern crate approx;
extern crate eigenvalues;
extern crate nalgebra as na;

use approx::relative_eq;
use eigenvalues::algorithms::davidson::EigenDavidson;
use eigenvalues::utils::generate_diagonal_dominant;
use eigenvalues::utils::sort_eigenpairs;
use eigenvalues::SpectrumTarget;

#[test]
fn test_davidson_lowest() {
    let arr = generate_diagonal_dominant(10, 0.005);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(arr.clone()), true);

    let spectrum_target = SpectrumTarget::Lowest;

    let dav_eig = EigenDavidson::new(arr.clone(), 2, "DPR", spectrum_target.clone()).unwrap();
    test_eigenpairs(&eig, dav_eig, 2);
    let dav_eig = EigenDavidson::new(arr.clone(), 2, "GJD", spectrum_target).unwrap();
    test_eigenpairs(&eig, dav_eig, 2);
}

#[test]
fn test_davidson_unsorted() {
    // Test the algorithm when the diagonal is unsorted
    let mut arr = generate_diagonal_dominant(8, 0.005);
    let vs = na::DVector::<f64>::from_vec(vec![3.0, 2.0, 4.0, 1.0, 5.0, 6.0, 7.0, 8.0]);
    arr.set_diagonal(&vs);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(arr.clone()), true);
    let dav_eig = EigenDavidson::new(arr, 1, "DPR", SpectrumTarget::Lowest).unwrap();
    test_eigenpairs(&eig, dav_eig, 1);
}

#[test]
fn test_davidson_highest() {
    // Test the compution of the highest eigenvalues
    let dim = 20;
    let nvalues = 2;
    let arr = generate_diagonal_dominant(dim, 0.005);
    let eig = sort_eigenpairs(na::linalg::SymmetricEigen::new(arr.clone()), false);

    let target = SpectrumTarget::Highest;
    println!("running DPR");
    let dav_eig = EigenDavidson::new(arr.clone(), nvalues, "DPR", target.clone()).unwrap();
    test_eigenpairs(&eig, dav_eig, nvalues);
    println!("running GJD");
    let dav_eig = EigenDavidson::new(arr.clone(), nvalues, "GJD", target).unwrap();
    test_eigenpairs(&eig, dav_eig, nvalues);
}

fn test_eigenpairs(
    reference: &na::linalg::SymmetricEigen<f64, na::Dynamic>,
    dav_eig: EigenDavidson,
    number: usize,
) {
    for i in 0..number {
        // Test Eigenvalues
        assert!(relative_eq!(
            reference.eigenvalues[i],
            dav_eig.eigenvalues[i],
            epsilon = 1e-6
        ));
        // Test Eigenvectors
        let x = reference.eigenvectors.column(i);
        let y = dav_eig.eigenvectors.column(i);
        // The autovectors may different in their sign
        // They should be either parallel or antiparallel
        let dot = x.dot(&y).abs();
        assert!(relative_eq!(dot, 1.0, epsilon = 1e-8));
    }
}
