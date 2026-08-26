# Numerical Derivatives of a Basis

The derivative method every basis inherits from the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
one finite-difference stencil of the order asked for, applied to
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md).
A subclass supplying its evaluation alone therefore answers
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
immediately, and registering a closed form later takes over through
dispatch with no change to calling code.

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

  The derivative order, a single non-negative whole number, or one per
  variable. More than one non-zero entry throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md).

## Details

Numerical Derivatives of a Basis

## Accuracy

One stencil of the order wanted, never a chain of first differences.
Measured against exact Legendre derivatives at 21 interior points of
\\\[0, 1\]\\, relative to the scale of the answer: 3.9e-10 at order 1,
1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08 at order 4. See
[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
for the stencil, the step and the endpoint rule.

## A basis of several variables

The stencil differentiates along one coordinate at a time, replacing
that column of the points and holding the others. A mixed partial such
as `c(1, 1)` therefore throws: a stencil in the plane carries the
product of two errors, and the one family that needs mixed partials,
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
computes them exactly from its margins.

## See also

[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md),
which does the work;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
to ask an object whether its derivatives come from here;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic.
