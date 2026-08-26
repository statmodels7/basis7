# Gram Matrix of a Tensor Product Basis

Returns the Kronecker product of the marginal Gram matrices, which is
the integral over the box exactly, with no quadrature. A product of
exactly integrated margins is therefore exact at any number of
variables, where a rule over the box would cost a node count exponential
in the number of variables and still carry an error.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- order:

  The derivative order, an integer vector with one entry per variable,
  or a single `0`. Each margin is asked for its own entry.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- ...:

  Passed to each margin's method.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
the product column names on both margins.

## Details

The integral over the box of a product of separable functions factorizes
into one-dimensional integrals, so \\G = G_1 \otimes \cdots \otimes
G_D\\ with \\G_j\\ the margin's own matrix at that margin's own order.
It costs one marginal Gram matrix per variable and no integration over
the box at all.

The result is symmetrized as `(G + t(G))/2`, a Kronecker product of
symmetric matrices being symmetric only up to the order its entries were
formed in. It is singular whenever any margin's is, so an `order` with
any non-zero entry gives a singular matrix.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures;
[`base::kronecker()`](https://rdrr.io/r/base/kronecker.html), whose
column order this follows.
