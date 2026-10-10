# Partial Derivatives of a Tensor Product Basis

Returns the mixed partial derivative named by a multi-index, one order
per variable: `c(2, 0)` is \\\partial^2/\partial x_1^2\\ and `c(1, 1)`
is \\\partial^2/\partial x_1 \partial x_2\\. Exact wherever the margins
are. It is the only method of the package that computes a mixed partial,
the numerical fallback signalling an error for one; a transformed tensor
basis obtains it from this method.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- order:

  An integer vector with one entry per variable, or a single `0`. A
  single non-zero order signals an error, having two readings; see
  [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns.

## Details

The product separates, so the derivative differentiates each margin to
its own order and multiplies the results. The derivative has no cross
term and needs no stencil in the plane, so the accuracy of a mixed
partial is the accuracy of the margins' own derivatives.

An order beyond what a margin carries makes that factor zero, so the
whole product is zero: `c(4, 0)` on a cubic B-spline margin is the exact
zero matrix whatever the other margins do.

## See also

[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md),
which does the work;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic and the multi-index rule.
