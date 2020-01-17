## [Unreleased]

## [0.1.3] 17-01-2020

### Fixed
* Fixed bug generating the subpsace when the diagonal is not sorted.

### Added

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
