# Integral of a Transformed Basis

The parent's integral, multiplied by the transform. The anchor at the
lower endpoint survives, since a linear combination of columns that are
all zero there is zero there.

## Arguments

- basis:

  A
  [`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows.
