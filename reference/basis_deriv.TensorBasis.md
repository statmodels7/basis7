# Partial Derivatives of a Tensor Product Basis

The mixed partial derivative given by a multi-index, one order per
variable.

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- order:

  An integer vector with one entry per variable.

- ...:

  Unused.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns.

## Details

The product separates, so the derivative differentiates each marginal to
its own order and multiplies the results. An order beyond what a
marginal supports makes that factor, and so the whole product, zero.
