# Construct a Legendre Polynomial Basis

The polynomials of degree below `dimension` on \\\[\ell, u\]\\,
expressed in the Legendre basis rather than in raw powers.

## Usage

``` r
poly_basis(lower = 0, upper = 1, dimension = 4)
```

## Arguments

- lower, upper:

  The endpoints of the interval.

- dimension:

  The number of polynomials, so the highest degree is `dimension - 1`.

## Value

An object of class
[`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md).

## Details

The Legendre polynomials span exactly the same space as \\1, x, x^2,
\ldots\\, and the choice between them is numerical. Raw powers give a
Gram matrix that is a Hilbert matrix, whose condition number grows so
fast that a basis of ten of them is already close to singular in double
precision. The Legendre polynomials are orthogonal, so their Gram matrix
is diagonal and perfectly conditioned.

## See also

[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)

## Examples

``` r
p <- poly_basis(dimension = 5)
p
#> Basis: legendre
#> Functions: 5   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree  4
#> Numerical: none

# orthogonal, so the Gram matrix is diagonal
round(basis_gram(p), 10)
#>    P0        P1  P2        P3        P4
#> P0  1 0.0000000 0.0 0.0000000 0.0000000
#> P1  0 0.3333333 0.0 0.0000000 0.0000000
#> P2  0 0.0000000 0.2 0.0000000 0.0000000
#> P3  0 0.0000000 0.0 0.1428571 0.0000000
#> P4  0 0.0000000 0.0 0.0000000 0.1111111

# and it spans the same space as the raw powers
x <- seq(0, 1, length.out = 40)
fitted <- lm.fit(basis_eval(p, x), x^3)$fitted.values
max(abs(fitted - x^3))
#> [1] 2.270252e-16
```
