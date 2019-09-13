extern crate davidson;
extern crate nalgebra as na;

use davidson::EigenDavidson;
use na::{DMatrix, DVector};


fn main() {
    let xs = 1..11;
    let dim = 10;
    let ys = xs.map(|x: i32| x as f64);
    let vs = DVector::<f64>::from_iterator(dim, ys);
    let mut brr = DMatrix::<f64>::new_random(dim, dim);
    brr += &brr.transpose();
    brr *= 0.005;
    brr.set_diagonal(&vs);
    println!("brr:{}", brr);
    // let mut eig = na::linalg::SymmetricEigen::new(brr);
    let eig = EigenDavidson::new(brr, 2).unwrap();
    println!("eigenvalues:{}", eig.eigenvalues);
    println!("eigenvectors:{}", eig.eigenvectors);

}
