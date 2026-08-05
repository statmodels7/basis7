# Linearly Transformed Basis

A basis obtained from another by a fixed linear map of its functions,
\\\tilde{B}(x) = B(x)\\T\\. Constructed by
[`orthonorm_basis`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`constrain_basis`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
or
[`dr_basis`](https://statmodels7.github.io/basis7/reference/dr_basis.md).

## Usage

``` r
TransformedBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list(),
  parent_basis = basis(),
  transform = integer(0)
)
```

## Arguments

- basis_name:

  A short name for the family, used when printing.

- dimension:

  The number of functions in the basis, a positive integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several.

- basis_params:

  A named list of whatever else the subclass needs.

- parent_basis:

  The basis being transformed.

- transform:

  The matrix \\T\\, with one row per parent function.

## Value

An object of class `TransformedBasis`.

## Details

Orthonormalizing a basis, restricting it to satisfy a linear constraint,
and rebuilding it so that it diagonalizes an inner product are the same
operation with different matrices, so they share one class. Derivatives
and integrals transform by the same \\T\\, because differentiation and
integration are linear and \\T\\ does not depend on \\x\\; the Gram
matrix transforms by congruence, \\T^\top G\\ T\\, so a parent with an
exact Gram matrix passes its exactness on.

\\T\\ may have fewer columns than rows, which is how a constraint
reduces the dimension.

Transforms compose by multiplication rather than by nesting:
transforming a `TransformedBasis` again produces one object holding the
product of the two matrices, so a chain of transforms costs one matrix
multiplication per evaluation however long it is.

## See also

[`orthonorm_basis`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`constrain_basis`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
[`dr_basis`](https://statmodels7.github.io/basis7/reference/dr_basis.md)

## Examples

``` r
o <- orthonorm_basis(bspline_basis(dimension = 6))
S7::S7_inherits(o, TransformedBasis)
#> [1] TRUE
dim(o@transform)
#> [1] 6 6
```
