/*!

# Modified Gram-Schmidt (MGS)

The Gram-Schmidt method is a method for orthonormalising a set of vectors. see:
[Gram-Schmidt process](https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process)
The MGS method improves the orthogonality loss due to the finite numerical precision
on computers.
 */

extern crate nalgebra as na;
use na::DMatrix;

pub struct MGS {
    pub basis: DMatrix<f64>,
}

impl MGS {
    /// The new static method takes a single argument:
    /// * `vectors` to diagonalize as columns of the matrix
    pub fn new(vectors: DMatrix<f64>) -> Result<Self, &'static str> {
        let mut result = Err("Something when wrong!");
        result
    }
}

#[cfg(test)]
mod test {
    extern crate nalgebra as na;
    use approx::relative_eq;
    use na::DMatrix;

    #[test]
    fn test_gram_schmidt() {
        let dim = 10;
        let vectors = DMatrix::<f64>::new_random(dim, dim);
        let mgs_result = super::MGS::new(vectors);
        let basis: DMatrix<f64> = match mgs_result {
            Ok(ortho) => ortho.basis,
            Err(message) => panic!(message),
        };

        let result = basis.transpose() * &basis;
        let diag = result.sum();
        assert!(relative_eq!(diag, dim as f64, epsilon = 1e-8));
    }
}
