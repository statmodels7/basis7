# Gram Matrix of a Basis

Returns the matrix of inner products of the `order`-th derivatives of
the basis functions, \$\$G\_{ab} = \int\_{\ell}^{u} B_a^{(d)}(t)\\
B_b^{(d)}(t)\\ \mathrm{d}t.\$\$

## Usage

``` r
basis_gram(basis, order = 0L, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order, a non-negative integer. Zero, the default, gives
  the inner products of the basis functions themselves.

- ...:

  Passed to methods.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

The Gram matrix is what a roughness penalty integrates: the quadratic
form \\\beta^\top G_2 \beta\\ is \\\int (f'')^2\\, for \\f\\ the
function the coefficients describe. It is basis mathematics, an inner
product, and says nothing about which combination a model should be
shrunk by.

It is symmetric and positive semidefinite by construction, and singular
whenever the order-\\d\\ derivatives are linearly dependent, which for
\\d \ge 1\\ they always are: constants differentiate to zero.

## See also

[`basis_deriv`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)

## Examples

``` r
b <- fourier_basis(dimension = 5)
round(basis_gram(b), 6) # diagonal, by orthogonality
#>       const sin1 cos1 sin2 cos2
#> const     1  0.0  0.0  0.0  0.0
#> sin1      0  0.5  0.0  0.0  0.0
#> cos1      0  0.0  0.5  0.0  0.0
#> sin2      0  0.0  0.0  0.5  0.0
#> cos2      0  0.0  0.0  0.0  0.5
round(basis_gram(b, order = 1), 4)
#>       const    sin1    cos1    sin2    cos2
#> const     0  0.0000  0.0000  0.0000  0.0000
#> sin1      0 19.7392  0.0000  0.0000  0.0000
#> cos1      0  0.0000 19.7392  0.0000  0.0000
#> sin2      0  0.0000  0.0000 78.9568  0.0000
#> cos2      0  0.0000  0.0000  0.0000 78.9568
```
