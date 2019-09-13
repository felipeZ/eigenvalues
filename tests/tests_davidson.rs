extern crate davidson;
extern crate nalgebra as na;

use davidson::EigenDavidson;
use na::linalg::SymmetricEigen;

#[test]
fn the_answer(){
    let answer = 42;
    assert_eq!(answer, 42);
}
