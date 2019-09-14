Eigenvalue Decomposition
========================
This package contains some implementations for computing the eigenvalues of a symmetric matrix,
implemented in [Rust](https://www.rust-lang.org/).

Available Algorithms:
 * **Davidson** using either  Diagonal-Preconditioned-Residue (**DPR**) or Generalized Jacobi Davidson (**GJD**).


### Note:
The Davidson method is suitable for **diagonal-dominant symmetric matrices**, that are quite common
in certain scientific problems like [electronic structure](https://en.wikipedia.org/wiki/Electronic_structure). The Davidson method could be not practical
for other kind of symmetric matrices.
