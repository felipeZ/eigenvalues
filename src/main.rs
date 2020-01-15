extern crate eigenvalues;
extern crate nalgebra as na;

use eigenvalues::davidson::EigenDavidson;

fn main() {
  let brr = eigenvalues::utils::generate_diagonal_dominant(10, 0.005);
  let eig = EigenDavidson::new (brr, 2).unwrap();
  println!("eigenvalues:{}", eig.eigenvalues);
  println!("eigenvectors:{}", eig.eigenvectors);
}
