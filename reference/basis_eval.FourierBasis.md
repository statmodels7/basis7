# Evaluate a Fourier Basis

The constant and the sine-cosine pairs, from the shift identity of
[`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
at order zero.

## Arguments

- basis:

  A
  [`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
