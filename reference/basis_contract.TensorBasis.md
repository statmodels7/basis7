# Contract a Tensor Product Basis Against Coefficients

The value of the function the coefficients describe, computed from the
marginal evaluations without forming the tensor design matrix.

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- coef:

  An array with one dimension per marginal, or a list of factor matrices
  in canonical polyadic form.

- block:

  The number of rows processed at once when `coef` is an array. It
  bounds the peak memory, which is otherwise what forming the design
  matrix would cost.

- ...:

  Unused.

## Value

A numeric vector with one value per row of `x`.
