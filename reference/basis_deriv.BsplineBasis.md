# Derivatives of a B-Spline Basis

Exact derivatives from the B-spline recurrence. An order above the
degree gives the zero matrix, which is the value of that derivative.

## Arguments

- basis:

  A
  [`BsplineBasis`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
