extern crate davidson;
extern crate nalgebra as na;

use davidson::EigenDavidson;
use na::{DMatrix, DVector};


fn main() {
    let vs = DVector::<f64>::from_iterator(5, [1., 2., 3., 4., 5.].iter().cloned());
    let mut brr = DMatrix::<f64>::new_random(5, 5);
    brr += &brr.transpose();
    brr *= 0.005;
    brr.set_diagonal(&vs);
    println!("brr:{}", brr);
    let mut eig = na::linalg::SymmetricEigen::new(brr);
    // let eig = EigenDavidson::new(brr, 2);
}
