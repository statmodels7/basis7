# Gram Matrix of a Tensor Product Basis

The Kronecker product of the marginal Gram matrices.

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- order:

  An integer vector with one entry per variable.

- at, weight:

  Handled by the generic before dispatch; unused here.

- ...:

  Passed to the marginals.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

The integral over the box of a product of separable functions factorises
into a product of one-dimensional integrals, so the matrix is separable
and costs one marginal Gram matrix per variable rather than one
integration over the box. A tensor of exactly integrated marginals is
therefore exact at any number of variables, where a quadrature over the
box would not be.
