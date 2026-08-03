# Restrict a Basis to a Linear Constraint

Returns the basis whose coefficient vectors are exactly those satisfying
\\C\beta = 0\\, with the dimension reduced by the rank of \\C\\.

## Usage

``` r
constrain_basis(basis, constraint, tol = 1e-10)
```

## Arguments

- basis:

  The basis to restrict.

- constraint:

  A matrix with one column per basis function, or a vector for a single
  constraint.

- tol:

  The relative tolerance below which a singular value counts as zero
  when determining the rank.

## Value

An object of class
[`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md).

## Details

The transform is an orthonormal basis of the null space of \\C\\,
extracted from its singular value decomposition, so the constrained
basis spans precisely the admissible functions and no reparametrisation
of the constraint changes the space it produces.

What the package supplies is the mechanics. Which constraint a model
term should carry, whether a sum-to-zero condition for identifiability
or orthogonality to a linear part, is a decision that needs to know what
the term means, and belongs to the layer that does.

## See also

[`dr_basis`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
[`orthonorm_basis`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)

## Examples

``` r
b <- bspline_basis(dimension = 6)

# sum to zero over a grid: the usual identifiability constraint
x <- seq(0, 1, length.out = 200)
cs <- constrain_basis(b, colSums(basis_eval(b, x)))
cs@dimension
#> [1] 5
max(abs(colSums(basis_eval(cs, x))))
#> [1] 2.034527e-14
```
