# Derivatives of a Legendre Basis

Exact derivatives of any order, from the derivative of the recurrence.
An order above the highest degree gives the zero matrix.

## Arguments

- basis:

  A
  [`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
