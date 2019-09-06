/*!

# Davidson Diagonalization

This crate contains an implementation of the Davidson's algorithm to
diagonalize symmetric matrices.

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

pub mod davidson;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
