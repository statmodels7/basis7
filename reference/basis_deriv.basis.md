# Numerical Derivatives of a Basis

The default derivative method, applying one finite-difference stencil to
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md).
It is what a basis that registers no derivative method of its own uses,
so that implementing the evaluation is enough to have a complete basis.

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

Numerical Derivatives of a Basis

See
[`numerical_deriv_matrix`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
for the stencil and the step, and
[`basis_is_numerical`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
for asking an object whether its derivatives come from here.
