# Derivatives of a Fourier Basis

Exact derivatives of any order, from the shift identity of
[`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md).
The constant differentiates to zero.

## Arguments

- basis:

  A
  [`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
