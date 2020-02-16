/*!

# Modified Gram-Schmidt (MGS)

The Gram-Schmidt method is a method for orthonormalising a set of vectors. see:
[Gram-Schmidt process](https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process)
The MGS method improves the orthogonality loss due to the finite numerical precision
on computers.
 */

extern crate nalgebra as na;
use na::{DMatrix, DVector, DVectorSlice};

pub struct MGS {
    pub basis: DMatrix<f64>,
}

impl MGS {
    /// The new static method takes three argument:
    /// * `vectors` to diagonalize as columns of the matrix
    /// * `start` index of the column to start orthogonalizing
    /// * `end` last index of the column to diagonalize (non-inclusive)
    pub fn new(vectors: DMatrix<f64>, start: usize, end: usize) -> Result<Self, &'static str> {
        let mut basis = vectors.clone();
        for i in start..end {
            basis.set_column(i, &vectors.column(i));
            for j in 0..i {
                let proj = MGS::project(&basis.column(j), &basis.column(i));
                basis.set_column(i, &(basis.column(i) - proj));
            }
            basis.set_column(i, &basis.column(i).normalize());
        }
        Ok(MGS { basis })
    }

    // Project
    fn project(v1: &DVectorSlice<f64>, v2: &DVectorSlice<f64>) -> DVector<f64> {
        let magnitud = v1.dot(&v2) / v1.dot(&v1);
        return v1 * magnitud;
    }
}

#[cfg(test)]
mod test {
    extern crate nalgebra as na;
    use na::DMatrix;

    #[test]
    fn test_gram_schmidt() {
        let dim = 10;
        let vectors = DMatrix::<f64>::new_random(dim, dim);
        fun_test(vectors, 0);
    }

    #[test]
    fn test_start_gram_schmidt() {
        let arr = DMatrix::<f64>::from_vec(3, 3, vec![1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 2.0]);
        fun_test(arr, 1);
    }

    fn fun_test(vectors: DMatrix<f64>, start: usize) {
        let end = vectors.ncols();
        let mgs_result = super::MGS::new(vectors, start, end);
        let basis = match mgs_result {
            Ok(ortho) => ortho.basis,
            Err(message) => panic!(message),
        };
        let result = basis.transpose() * &basis;
        assert!(result.is_identity(1e-8));
    }
}
