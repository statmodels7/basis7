# Evaluate a Legendre Basis

The three-term recurrence, on the shifted variable.

## Arguments

- basis:

  A
  [`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
