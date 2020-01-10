/*!

# Davidson Diagonalization

The Davidson method is suitable for diagonal-dominant symmetric matrices,
that are quite common in certain scientific problems like [electronic structure](https://en.wikipedia.org/wiki/Electronic_structure).
The Davidson method could be not practical for other kind of symmetric matrices.

The current implementation uses a general davidson algorithm, meaning
that it compute all the requested eigenvalues simultaneusly using a variable size
 block approach.
The family of Davidson algorithm only differ in the way that the correction
vector is computed.

Available correction methods are:
 * **DPR**: Diagonal-Preconditioned-Residue
 * **GJD**: Generalized Jacobi Davidson

*/

extern crate nalgebra as na;
use crate::matrix_operations::MatrixOperations;
use crate::utils;
use na::linalg::SymmetricEigen;
use na::{DMatrix, DVector, DVectorSlice, Dynamic};
use std::f64;

/// Structure containing the initial configuration data
struct Config {
    method: String,
    tolerance: f64,
    max_iters: usize,
    max_search_space: usize,
    init_dim: usize,
}
impl Config {
    /// Choose sensible default values for the davidson algorithm, where:
    /// * `nvalues` - Number of eigenvalue/eigenvector pairs to compute
    /// * `dim` - dimension of the matrix to diagonalize
    fn new(nvalues: usize, dim: usize) -> Self {
        let max_search_space = if nvalues * 10 < dim {
            nvalues * 10
        } else {
            dim
        };
        Config {
            method: String::from("DPR"),
            tolerance: 1e-8,
            max_iters: 100,
            max_search_space: max_search_space,
            init_dim: nvalues * 2,
        }
    }
}

/// Structure with the configuration data
pub struct EigenDavidson {
    pub eigenvalues: DVector<f64>,
    pub eigenvectors: DMatrix<f64>,
}

impl EigenDavidson {
    /// The new static method takes the following arguments:
    /// * `h` - A highly diagonal symmetric matrix
    /// * `nvalues` - the number of eigenvalues/eigenvectors pair to compute

    pub fn new<M: MatrixOperations>(h: M, nvalues: usize) -> Result<Self, &'static str> {
        // Initial configuration
        let conf = Config::new(nvalues, h.rows());

        // Initial subpace
        let mut dim_sub = conf.init_dim;
        // 1.1 Select the initial ortogonal subspace based on lowest elements
        let mut basis = generate_subspace(&h.diagonal(), conf.max_search_space);

        // 1.2 Vector containing the indices of th converged and deflated eigenpairs
        // let mut deflated = Vec::new();

        // Outer loop block Davidson schema
        let mut result = Err("Algorithm didn't converge!");
        for i in 0..conf.max_iters {
            // 2. Generate subpace matrix problem by projecting into the basis
            let subspace = basis.columns(0, dim_sub);
            let matrix_proj = subspace.transpose() * &h.matrix_matrix_prod(subspace); // (&h * subspace);

            // 3. compute the eigenvalues and their corresponding ritz_vectors
            let eig = utils::sort_eigenpairs(SymmetricEigen::new(matrix_proj));

            // 4. Check for convergence
            // 4.1 Compute the residues
            let ritz_vectors = subspace * eig.eigenvectors.columns(0, dim_sub);
            let residues = compute_residues(&h, &eig.eigenvalues, &ritz_vectors);

            // 4.2 Check Converge for each pair eigenvalue/eigenvector
            let errors = DVector::<f64>::from_iterator(
                nvalues,
                residues
                    .columns(0, nvalues)
                    .column_iter()
                    .map(|col| col.norm()),
            );

            // 4.3 Check if all eigenvalues/eigenvectors have converged
            if errors.iter().all(|&x| x < conf.tolerance) {
                result = Ok(create_results(&eig.eigenvalues, &ritz_vectors, nvalues));
                break;
            }

            // 4.4 Deflates the converged eigenvalues
            // println!("Iter:{}", i);
            // println!("errors:{}", errors);
            // for (k, err) in errors.iter().enumerate() {
            //     if (*err < conf.tolerance) & (!deflated.contains(&k)) {
            //         println!("deflating:{}", k);
            //         deflate_converged(&mut h, eig.eigenvalues[k], &eig.eigenvectors.column(k));
            //         deflated.push(i);
            //     }
            // }

            // 5. Update subspace basis set
            // 5.1 Add the correction vectors to the current basis
            if 2 * dim_sub <= conf.max_search_space {
                let correction = compute_correction(&h, residues, eig, &conf.method);
                update_subspace(&mut basis, correction, dim_sub, dim_sub * 2);

                // 6. Orthogonalize the subspace
                basis = orthogonalize_subspace(basis);
                // update counter
                dim_sub *= 2;

            // 5.2 Otherwise reduce the basis of the subspace to the current correction
            } else {
                dim_sub = conf.init_dim;
                basis.fill(0.0);
                update_subspace(&mut basis, ritz_vectors, 0, dim_sub);
            }
            // Check number of iterations
            if i > conf.max_iters {
                break;
            }
        }
        result
    }
}

/// Deflates the converged eigenvalue/eigenvector
fn deflate_converged(h: &mut DMatrix<f64>, lambda: f64, vs: &DVectorSlice<f64>) {
    let dim = h.nrows();
    let xs = DVector::<f64>::from_iterator(dim, vs.iter().cloned());
    let arr = lambda * xs.transpose() * xs;
    h.zip_apply(&arr, |x, y| x - y);
}

/// Extract the requested eigenvalues/eigenvectors pairs
fn create_results(
    subspace_eigenvalues: &DVector<f64>,
    ritz_vectors: &DMatrix<f64>,
    nvalues: usize,
) -> EigenDavidson {
    let eigenvectors = DMatrix::<f64>::from_iterator(
        ritz_vectors.nrows(),
        nvalues,
        ritz_vectors.columns(0, nvalues).iter().cloned(),
    );
    let eigenvalues = DVector::<f64>::from_iterator(
        nvalues,
        subspace_eigenvalues.rows(0, nvalues).iter().cloned(),
    );
    EigenDavidson {
        eigenvalues,
        eigenvectors,
    }
}

/// Update the subpace with new vectors
fn update_subspace(basis: &mut DMatrix<f64>, vectors: DMatrix<f64>, start: usize, end: usize) {
    let mut i = 0; // indices for the new vector to add
    for k in start..end {
        basis.set_column(k, &vectors.column(i));
        i += 1;
    }
}

/// Orthogonalize the subpsace using the QR method
fn orthogonalize_subspace(basis: DMatrix<f64>) -> DMatrix<f64> {
    let qr = na::linalg::QR::new(basis);
    qr.q()
}

/// Residue vectors
fn compute_residues<M: MatrixOperations>(
    h: &M,
    eigenvalues: &DVector<f64>,
    ritz_vectors: &DMatrix<f64>,
) -> DMatrix<f64> {
    let dim_sub = eigenvalues.nrows();
    let mut residues = DMatrix::<f64>::zeros(h.rows(), dim_sub);
    for k in 0..dim_sub {
        let guess = eigenvalues[k] * ritz_vectors.column(k);
        let vs = h.matrix_vector_prod(ritz_vectors.column(k));
        residues.set_column(k, &(vs - guess));
    }
    residues
}

/// compute the correction vectors using either DPR or GJD
fn compute_correction<M: MatrixOperations>(
    h: &M,
    residues: DMatrix<f64>,
    eigenpairs: SymmetricEigen<f64, Dynamic>,
    method: &str,
) -> DMatrix<f64> {
    match method.to_uppercase().as_ref() {
        "DPR" => compute_dpr_correction(h, residues, &eigenpairs.eigenvalues),
        "GJD" => compute_gjd_correction(h, residues, &eigenpairs),
        _ => panic!("Method {} has not been implemented", method),
    }
}

/// Use the Diagonal-Preconditioned-Residue (DPR) method to compute the correction
fn compute_dpr_correction<M: MatrixOperations>(
    h: &M,
    residues: DMatrix<f64>,
    eigenvalues: &DVector<f64>,
) -> DMatrix<f64> {
    let d = h.diagonal();
    let mut correction = DMatrix::<f64>::zeros(h.rows(), residues.ncols());
    for (k, r) in eigenvalues.iter().enumerate() {
        let rs = DVector::<f64>::repeat(h.rows(), *r);
        let x = residues.column(k).component_mul(&(rs - &d));
        correction.set_column(k, &x);
    }
    correction
}

/// Use the Generalized Jacobi Davidson (GJD) to compute the correction
fn compute_gjd_correction<M: MatrixOperations>(
    h: &M,
    residues: DMatrix<f64>,
    eigenpairs: &SymmetricEigen<f64, Dynamic>,
) -> DMatrix<f64> {
    // let dimx = h.rows();
    // let dimy = eigenpairs.eigenvalues.nrows();
    // let id = DMatrix::<f64>::identity(dimx, dimx);
    // let ones = DVector::<f64>::repeat(dimx, 1.0);
    // let mut correction = DMatrix::<f64>::zeros(dimx, dimy);
    // for (k, r) in eigenpairs.eigenvalues.iter().enumerate() {
    //     // Create the components of the linear system
    //     let x = eigenpairs.eigenvectors.column(k);
    //     let t1 = &id - x * x.transpose();
    //     let mut t2 = h.clone();
    //     t2.set_diagonal(&(*r * &ones));
    //     let arr = &t1 * t2 * &t1;
    //     // Solve the linear system
    //     let decomp = arr.lu();
    //     let mut b = -residues.column(k);
    //     decomp.solve_mut(&mut b);
    //     correction.set_column(k, &b);
    // }
    let correction = DMatrix::<f64>::zeros(2, 2);
    correction
}

/// Generate initial orthonormal subspace
fn generate_subspace(diag: &DVector<f64>, max_search_space: usize) -> DMatrix<f64> {
    if is_sorted(diag) {
        DMatrix::<f64>::identity(diag.nrows(), max_search_space)
    } else {
        // TODO implement the case when the diagonal is not sorted
        panic!("Matrix diagonal elements are not sorted")
    }
}

/// Check if a vector is sorted in ascending order
fn is_sorted(xs: &DVector<f64>) -> bool {
    let mut d: Vec<f64> = xs.iter().cloned().collect();
    d.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let vs: DVector<f64> = DVector::<f64>::from_vec(d);
    let r = xs - vs;
    r.norm() < f64::EPSILON
}

#[cfg(test)]
mod test {
    extern crate nalgebra as na;
    use na::{DMatrix, DVector};

    #[test]
    fn test_update_subspace() {
        let mut arr = DMatrix::<f64>::repeat(3, 3, 1.);
        let brr = DMatrix::<f64>::zeros(3, 2);
        super::update_subspace(&mut arr, brr, 0, 2);
        assert_eq!(arr.column(1).sum(), 0.);
        assert_eq!(arr.column(2).sum(), 3.);
    }

}
