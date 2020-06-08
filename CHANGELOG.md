## [0.2.1] 08/06/2020

### New
* Add Lanczos algorithm (#12)

### Fixed
* Only orthogonalize the new vectors added to the subpace.
* Add Eigendavidson argument documentation.

### Changed
* Renamed struc `Eigendavidson` => `Davidson`.

## [0.2.0] 28-02-2029

## Changed
* Used [modified Gram-Schmidt orthogonalization](https://github.com/felipeZ/eigenvalues/issues/8)
* [Reduce scaling of the search space](https://github.com/felipeZ/eigenvalues/issues/10)

### Fixed
* Bug in DPR correction.

## [0.1.4]

### Changed

* The Davidson method requires a **SpectrumTarget** `enum` argument.

## [0.1.3] 17-01-2020

### Fixed
* Fixed bug generating the subpsace when the diagonal is not sorted.

### Added
* Allow to compute the lowest or highest eigenvalues of the spectrum.


## [0.1.2] 17-01-2020

### Added
* Generate the orthonormal subspace different from the identity matrix when the diagonal values are not sorted.

## [0.1.1] 16-01-2020

### Changed
* Declared **Clone** a [Super Trait](https://doc.rust-lang.org/reference/items/traits.html#supertraits) **Matrixoperations**.

## [0.1.0] 16-01-2020

### Added
* First implementation using `f64`.
* Used [nalgebra](https://github.com/rustsim/nalgebra) for the linear algebra operations. 
* Both **GDJ** and **DPR** correction methods working.
* Use a Trait to collect the methods that a matrix representation should have to compute the eigenvalue/eigenvectors.
