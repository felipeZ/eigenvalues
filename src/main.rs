extern crate eigenvalues;
extern crate nalgebra as na;

use eigenvalues::davidson::Davidson;
use eigenvalues::SpectrumTarget;

fn main() {
    let brr = eigenvalues::utils::generate_diagonal_dominant(20, 0.05);
    let tolerance = 1e-6;
    let eig = Davidson::new(brr, 2, "GJD", SpectrumTarget::Highest, tolerance).unwrap();
    println!("eigenvalues:{}", eig.eigenvalues);
    println!("eigenvectors:{}", eig.eigenvectors);
}
