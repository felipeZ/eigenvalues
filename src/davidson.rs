//! ### Davidson implementation
extern crate nalgebra as na;
use na::{DMatrix, DVector};

/// Structure with the configuration data
struct Data {
    nvalues: usize,
    max__iters: usize,
    max_dim_sub: usize,
}

/// API
pub fn davidson(
    H: DMatrix<f64>,
    nvalues: usize,
    max__iters: usize,
) -> (DVector<f64>, DMatrix<f64>) {
    // Data used for the algorithm
    let max_dim_sub = nvalues * 10;
    let data = Data {
        nvalues,
        max__iters,
        max_dim_sub,
    };

    // Initial subpace
    let dim_sub = nvalues * 2;
    let V = generate_subspace(&H, dim_sub);

    // Outer loop block Davidson schema
    for i in 1..max__iters {
        println!("wheee, {}", i);
    }

    let es = DVector::<f64>::zeros(3);
    let vs = DMatrix::<f64>::identity(3, 3);
    (es, vs)
}

/// Generate initial orthonormal subspace
fn generate_subspace(mtx: &DMatrix<f64>, max_dim_sub: usize) -> DMatrix<f64> {
    let d = mtx.diagonal();
    DMatrix::<f64>::identity(mtx.nrows(), max_dim_sub)
}

#[cfg(test)]
mod test {

    #[test]
    fn test_something() {
        let universe = 42;
        assert_eq!(universe, 42);
    }

}
