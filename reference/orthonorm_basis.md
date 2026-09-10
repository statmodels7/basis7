# Orthonormalize a Basis

Returns a basis spanning the same functions whose Gram matrix is the
identity, so that the design matrix has orthogonal columns of unit norm
in \\L^2\\ and a least-squares fit against it is perfectly conditioned.
The transform is read off the Gram matrix, so for the shipped families
the orthonormalization is exact: there is no grid, no number of points
to choose, and no scale factor to correct.

## Usage

``` r
orthonorm_basis(basis, order = 0L)
```

## Arguments

- basis:

  The basis to orthonormalize, any object inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- order:

  The derivative order whose inner products are made the identity, a
  single non-negative whole number, default `0`. Any value above `0`
  throws for the shipped families, their higher-order Gram matrices
  being singular.

## Value

An object of class
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
with `basis_name` `orthonorm(<parent>)` and column names `on1`, `on2`,
and so on. Its dimension is the parent's.

## The construction

Write \\G = R^\top R\\ for the Cholesky factorization of the Gram
matrix. The basis \\B R^{-1}\\ then has Gram matrix \\R^{-\top} R^\top
R\\ R^{-1} = I\\, so \\T = R^{-1}\\. Measured on a cubic B-spline of six
functions,
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
of the result differs from the identity by 4.4e-16.

The span is unchanged, \\R^{-1}\\ being invertible: a function the
parent can represent is fitted by the orthonormalized basis to rounding.

## Orthonormal in which inner product

`order` chooses it. At `0`, the default, the functions themselves are
orthonormal. Above that the `order`-th derivatives would be, and the
Gram matrix there is singular for every family, the constant
differentiating away, so the factorization fails and an error says so.
Orthonormalizing a derivative therefore needs a basis whose constant has
already been removed by
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md).

## Composing

Orthonormalizing an already orthonormal basis returns it unchanged up to
rounding, and the two transforms collapse into one matrix; only
`@basis_name` records the second pass.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
which supplies \\G\\ and which the result makes the identity;
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
to remove a direction first;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
for a rotation that diagonalizes two matrices at once.

## Examples

``` r
b <- bspline_basis(dimension = 6)
o <- orthonorm_basis(b)

# The Gram matrix is the identity, exactly.
round(basis_gram(o), 12)
#>     on1 on2 on3 on4 on5 on6
#> on1   1   0   0   0   0   0
#> on2   0   1   0   0   0   0
#> on3   0   0   1   0   0   0
#> on4   0   0   0   1   0   0
#> on5   0   0   0   0   1   0
#> on6   0   0   0   0   0   1
max(abs(basis_gram(o) - diag(6)))
#> [1] 4.440892e-16

# The span is unchanged: a function the parent represents is fitted exactly.
set.seed(1)
x <- seq(0, 1, length.out = 100)
f <- drop(basis_eval(b, x) %*% rnorm(6))
max(abs(lm.fit(basis_eval(o, x), f)$residuals))
#> [1] 5.533747e-16

# Orthonormalizing again changes nothing, and keeps one matrix.
max(abs(basis_gram(orthonorm_basis(o)) - diag(6)))
#> [1] 4.440892e-16
class(orthonorm_basis(o)@parent_basis)
#> [1] "basis7::BsplineBasis" "basis7::basis"        "S7_object"           

# A higher order is refused: that Gram matrix is singular.
try(orthonorm_basis(b, order = 2))
#> Error : The Gram matrix is singular, so the basis functions are linearly dependent and cannot be orthonormalized. Reduce 'dimension', or orthonormalize at order 0.
```
