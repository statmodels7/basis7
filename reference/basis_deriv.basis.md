# Numerical Derivatives of a Basis

The derivative method every basis inherits from the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
one finite-difference stencil of the requested order, applied to
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md).
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
is therefore available for a subclass that supplies its evaluation
alone, and a closed form registered later takes over through dispatch
with no change to calling code.

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  Evaluation points inside the basis interval: a numeric vector for a
  basis of one variable, or a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns for a basis of several.

- order:

  The derivative order, already checked by the generic: a single
  non-negative whole number, or one per variable. More than one non-zero
  entry signals an error.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md).

## Details

Numerical Derivatives of a Basis

## Accuracy

One stencil of the requested order is used, never a chain of first
differences. See
[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
for the stencil, the step, the endpoint rule and the measured accuracy.

## A basis of several variables

The stencil differentiates along one coordinate at a time, replacing
that column of the points and holding the others. A mixed partial such
as `c(1, 1)` therefore signals an error: a stencil in the plane carries
the product of two errors, and
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
the family that needs mixed partials, computes them exactly from its
margins.

## See also

[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md),
which does the work;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md),
which reports whether the derivatives of an object come from here;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic.
