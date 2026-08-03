# Derivatives of the Legendre Polynomials

The `order`-th derivative of \\P_0, \ldots, P\_{K-1}\\ with respect to
`x`, on the basis interval.

## Usage

``` r
legendre_deriv(basis, x, order)
```

## Arguments

- basis:

  A
  [`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

Differentiating the identity \\P\_{n+1}' - P\_{n-1}' = (2n+1) P_n\\
repeatedly gives \\P\_{n+1}^{(d)} = P\_{n-1}^{(d)} + (2n+1)
P_n^{(d-1)}\\, which builds each order from the one below it in closed
form, with no differencing anywhere. The chain rule for the shift to
\\\[\ell, u\]\\ contributes \\(2/(u-\ell))^{d}\\.
