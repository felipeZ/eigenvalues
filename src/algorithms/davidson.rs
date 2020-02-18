/*!

# Davidson Diagonalization

The Davidson method is suitable for diagonal-dominant symmetric matrices,
that are quite common in certain scientific problems like [electronic
structure](https://en.wikipedia.org/wiki/Electronic_structure). The Davidson
method could be not practical for other kind of symmetric matrices.

The current implementation uses a general davidson algorithm, meaning
that it compute all the requested eigenvalues simultaneusly using a variable
size block approach. The family of Davidson algorithm only differ in the way
that the correction vector is computed.

Available correction methods are:
 * **DPR**: Diagonal-Preconditioned-Residue
 * **GJD**: Generalized Jacobi Davidson

*/

extern crate nalgebra as na;
use super::SpectrumTarget;
use crate::matrix_operations::MatrixOperations;
use crate::utils;
use crate::MGS;
use na::linalg::SymmetricEigen;
use na::{DMatrix, DVector};
use std::f64;
use std::ops::Not;

/// Structure containing the initial configuration data
struct Config {
    method: String,
    spectrum_target: SpectrumTarget,
    tolerance: f64,
    max_iters: usize,
    max_search_space: usize,
    init_dim: usize,
}
impl Config {
    /// Choose sensible default values for the davidson algorithm, where:
    /// * `nvalues` - Number of eigenvalue/eigenvector pairs to compute
    /// * `dim` - dimension of the matrix to diagonalize
    /// * `method` - Either DPR or GJD
    /// * `target` Lowest, highest or somewhere in the middle portion of the spectrum
    fn new(nvalues: usize, dim: usize, method: &str, target: SpectrumTarget) -> Self {
        let max_search_space = if nvalues * 10 < dim {
            nvalues * 10
        } else {
            dim
        };
        // Check that the correction method requested by the user is available
        let available_methods = ["DPR", "GJD"];
        if available_methods
            .iter()
            .any(|&x| x == method.to_uppercase())
            .not()
        {
            panic!("Method {} is not available!", method);
        }
        Config {
            method: String::from(method),
            spectrum_target: target,
            tolerance: 1e-4,
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

    pub fn new<M: MatrixOperations>(
        h: M,
        nvalues: usize,
        method: &str,
        spectrum_target: SpectrumTarget,
    ) -> Result<Self, &'static str> {
        // Initial configuration
        let conf = Config::new(nvalues, h.rows(), method, spectrum_target);

        // Initial subpace
        let mut dim_sub = conf.init_dim;
        // 1.1 Select the initial ortogonal subspace based on lowest elements
        let mut basis = generate_subspace(&h.diagonal(), &conf);

        // 1.2 Select the correction to use
        let corrector = CorrectionMethod::<M>::new(&h, &conf.method);

        // Outer loop block Davidson schema
        let mut result = Err("Algorithm didn't converge!");
        for i in 0..conf.max_iters {
            // 2. Generate subpace matrix problem by projecting into the basis
            let subspace = basis.columns(0, dim_sub);
            let matrix_proj = subspace.transpose() * &h.matrix_matrix_prod(subspace);

            // 3. compute the eigenvalues and their corresponding ritz_vectors
            let ord_sort = match conf.spectrum_target {
                SpectrumTarget::Highest => false,
                _ => true,
            };
            let eig = utils::sort_eigenpairs(SymmetricEigen::new(matrix_proj), ord_sort);

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
            // 5. Update subspace basis set
            // 5.1 Add the correction vectors to the current basis
            if dim_sub + nvalues <= conf.max_search_space {
                let correction =
                    corrector.compute_correction(residues, &eig.eigenvalues, &ritz_vectors);
                update_subspace(&mut basis, correction, dim_sub, dim_sub + nvalues);

                // 6. Orthogonalize the subspace
                basis = orthogonalize_subspace(basis, 0);
                // update counter
                dim_sub += nvalues;

            // 5.2 Otherwise reduce the basis of the subspace to the current
            // correction
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

/// Structure containing the correction methods
struct CorrectionMethod<'a, M>
where
    M: MatrixOperations,
{
    /// The initial target matrix
    target: &'a M,
    /// Method used to compute the correction
    method: String,
}

impl<'a, M> CorrectionMethod<'a, M>
where
    M: MatrixOperations,
{
    fn new(target: &'a M, method: &str) -> CorrectionMethod<'a, M> {
        CorrectionMethod {
            target: target,
            method: String::from(method),
        }
    }

    /// compute the correction vectors using either DPR or GJD
    fn compute_correction(
        &self,
        residues: DMatrix<f64>,
        eigenvalues: &DVector<f64>,
        ritz_vectors: &DMatrix<f64>,
    ) -> DMatrix<f64> {
        match self.method.to_uppercase().as_ref() {
            "DPR" => self.compute_dpr_correction(residues, eigenvalues),
            "GJD" => self.compute_gjd_correction(residues, eigenvalues, ritz_vectors),
            _ => panic!("Method {} has not been implemented", self.method),
        }
    }

    /// Use the Diagonal-Preconditioned-Residue (DPR) method to compute the correction
    fn compute_dpr_correction(
        &self,
        residues: DMatrix<f64>,
        eigenvalues: &DVector<f64>,
    ) -> DMatrix<f64> {
        let d = self.target.diagonal();
        let mut correction = DMatrix::<f64>::zeros(self.target.rows(), residues.ncols());
        for (k, lambda) in eigenvalues.iter().enumerate() {
            let tmp = DVector::<f64>::repeat(self.target.rows(), *lambda) - &d;
            let rs = residues.column(k).component_div(&tmp);
            correction.set_column(k, &rs);
        }
        correction
    }

    /// Use the Generalized Jacobi Davidson (GJD) to compute the correction
    fn compute_gjd_correction(
        &self,
        residues: DMatrix<f64>,
        eigenvalues: &DVector<f64>,
        ritz_vectors: &DMatrix<f64>,
    ) -> DMatrix<f64> {
        let dimx = self.target.rows();
        let dimy = residues.ncols();
        let id = DMatrix::<f64>::identity(dimx, dimx);
        let ones = DVector::<f64>::repeat(dimx, 1.0);
        let mut correction = DMatrix::<f64>::zeros(dimx, dimy);
        let diag = self.target.diagonal();
        for (k, r) in ritz_vectors.column_iter().enumerate() {
            // Create the components of the linear system
            let t1 = &id - r * r.transpose();
            let mut t2 = self.target.clone();
            let val = &diag - &(eigenvalues[k] * &ones);
            t2.set_diagonal(&val);
            let arr = &t1 * &t2.matrix_matrix_prod(t1.rows(0, dimx));
            // Solve the linear system
            let decomp = arr.lu();
            let mut b = -residues.column(k);
            decomp.solve_mut(&mut b);
            correction.set_column(k, &b);
        }
        correction
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
fn orthogonalize_subspace(vectors: DMatrix<f64>, start: usize) -> DMatrix<f64> {
    let mgs = MGS::new(vectors, start);
    match mgs {
        Ok(result) => result.basis,
        Err(msg) => panic!("Error orthonormalising the basis:{}", msg),
    }
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

/// Generate initial orthonormal subspace
fn generate_subspace(diag: &DVector<f64>, conf: &Config) -> DMatrix<f64> {
    if is_sorted(diag) {
        DMatrix::<f64>::identity(diag.nrows(), conf.max_search_space)
    } else {
        let xs = diag.as_slice().to_vec();
        let mut rs = xs.clone();

        match conf.spectrum_target {
            SpectrumTarget::Lowest => utils::sort_vector(&mut rs, true),
            SpectrumTarget::Highest => utils::sort_vector(&mut rs, false),
            _ => panic!("Not implemented error!"),
        }

        // update the matrix according to the spectrumtarget
        let mut mtx = DMatrix::<f64>::zeros(diag.nrows(), conf.max_search_space);
        for i in 0..conf.max_search_space {
            let index = rs.iter().position(|&x| x == xs[i]).unwrap();
            mtx[(i, index)] = 1.0;
        }
        mtx
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
    use na::DMatrix;

    #[test]
    fn test_update_subspace() {
        let mut arr = DMatrix::<f64>::repeat(3, 3, 1.);
        let brr = DMatrix::<f64>::zeros(3, 2);
        super::update_subspace(&mut arr, brr, 0, 2);
        assert_eq!(arr.column(1).sum(), 0.);
        assert_eq!(arr.column(2).sum(), 3.);
    }
}
