# Orthonormalise a Basis

Returns the basis whose functions span the same space and are
orthonormal in \\L^2\\ on the basis interval.

## Usage

``` r
orthonorm_basis(basis, order = 0L)
```

## Arguments

- basis:

  The basis to orthonormalise.

- order:

  The derivative order whose inner products are made the identity. Zero,
  the default, orthonormalises the functions themselves.

## Value

An object of class
[`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md).

## Details

The transform is read off the Gram matrix rather than estimated on a
grid. Writing \\G = R^\top R\\ for the Cholesky factorisation of the
Gram matrix, the basis \\B R^{-1}\\ has Gram matrix \\R^{-\top} R^\top R
R^{-1} = I\\. Because the parent's Gram matrix is exact for the families
that ship with the package, so is the orthonormalisation: there is no
grid, no number of points to choose, and no scale factor to correct.

Orthonormalising an already orthonormal basis returns it unchanged, up
to rounding, and the transforms collapse rather than nesting.

## See also

[`basis_gram`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
[`constrain_basis`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)

## Examples

``` r
o <- orthonorm_basis(bspline_basis(dimension = 6))
round(basis_gram(o), 12) # the identity, exactly
#>     on1 on2 on3 on4 on5 on6
#> on1   1   0   0   0   0   0
#> on2   0   1   0   0   0   0
#> on3   0   0   1   0   0   0
#> on4   0   0   0   1   0   0
#> on5   0   0   0   0   1   0
#> on6   0   0   0   0   0   1
```
