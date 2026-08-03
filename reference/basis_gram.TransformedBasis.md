# Gram Matrix of a Transformed Basis

The congruence \\T^\top G\\ T\\ of the parent's Gram matrix, so a parent
whose inner products are exact passes that on rather than falling back
to quadrature.

## Arguments

- basis:

  A
  [`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- order:

  The derivative order.

- at, weight:

  Handled by the generic before dispatch; unused here.

- ...:

  Passed to the parent's method.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.
