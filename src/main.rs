extern crate davidson;
extern crate nalgebra as na;

use na::{DMatrix, DVector};

fn main() {
  let vs = DVector::<f64>::from_iterator(3, [ 1., 3., 5. ].iter().cloned());
  let mut brr = DMatrix::<f64>::new_random(3, 3);
  brr += &brr.transpose();
  brr *= 0.005;
  brr.set_diagonal(&vs);
}
