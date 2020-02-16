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
    /// The new static method takes a single argument:
    /// * `vectors` to diagonalize as columns of the matrix
    pub fn new(vectors: DMatrix<f64>) -> Result<Self, &'static str> {
        if !vectors.is_square() {
            return Err("vectors must be an square matrix");
        }
        let dim = vectors.nrows();
        let mut basis = DMatrix::<f64>::zeros(dim, dim);
        // first vector of the basis
        let first = vectors.column(0) / vectors.column(0).norm();
        basis.set_column(0, &first);
        // Iterate over the rest of the columns
        for i in 1..basis.ncols() {
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
        let mgs_result = super::MGS::new(vectors);
        let basis: DMatrix<f64> = match mgs_result {
            Ok(ortho) => ortho.basis,
            Err(message) => panic!(message),
        };

        let result = basis.transpose() * &basis;
        assert!(result.is_identity(1e-8));
    }
}
