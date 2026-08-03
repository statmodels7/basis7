# Integral of a Tensor Product Basis

The integral over the box from the lower corner to each point.

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

## Details

The integrand separates, so the multiple integral is the product of the
marginal integrals, and the anchor survives: a product in which every
factor is zero at the corner is zero at the corner.
