# Gram Matrix of a Basis

Returns the matrix of inner products of the `order`-th derivatives of
the basis functions, \$\$G\_{ab} = \int\_{\ell}^{u} B_a^{(d)}(t)\\
B_b^{(d)}(t)\\ \mathrm{d}t.\$\$

## Usage

``` r
basis_gram(basis, order = 0L, at = NULL, weight = NULL, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order, a non-negative integer. Zero, the default, gives
  the inner products of the basis functions themselves.

- at:

  An optional numeric vector of points. When given, the inner products
  are taken against the empirical measure of those points rather than
  against Lebesgue measure.

- weight:

  An optional function of one numeric vector, a density to weight the
  integral by.

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

The inner product is taken against a measure, and which measure matters.
The default is Lebesgue on the basis interval, which is what a roughness
penalty integrates. Supplying `at` takes the empirical measure of those
points instead, \\B^\top B / n\\, which is the matrix a design matrix
actually produces and the one a basis is diagonalized against when the
construction is meant to depend on where the data lie. Supplying
`weight` takes a weighted Lebesgue measure.

Both alternatives are handled in the body of the generic, before
dispatch, so a method never sees them and never has to implement them;
it always returns the plain Lebesgue matrix. A method must still name
the arguments in its signature, because S7 requires a method's formals
to contain the generic's.

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

# against the empirical measure of a sample instead
set.seed(1)
round(basis_gram(b, at = runif(2000)), 3)
#>        const   sin1   cos1   sin2   cos2
#> const  1.000  0.021  0.019 -0.021  0.011
#> sin1   0.021  0.494 -0.010  0.010  0.010
#> cos1   0.019 -0.010  0.506  0.031  0.009
#> sin2  -0.021  0.010  0.031  0.503 -0.010
#> cos2   0.011  0.010  0.009 -0.010  0.497
```
