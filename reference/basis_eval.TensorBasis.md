# Evaluate a Tensor Product Basis

The row-wise Kronecker product of the marginal evaluations.

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- ...:

  Unused.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns.
