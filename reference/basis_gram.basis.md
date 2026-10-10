# Numerical Gram Matrix of a Basis

The inner-product method every basis inherits from the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
composite Gauss-Legendre over `panels` equal subintervals of the
interval, applied to the requested derivatives. A one-line wrapper over
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
which does the work and which
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
also calls when the period is not the width of the interval.

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- order:

  The derivative order, a single non-negative whole number, default `0`,
  or one per variable.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- panels:

  The number of equal subintervals, default `50`. For a basis of several
  variables the rule is a product and each coordinate gets
  `max(2, ceiling(panels^(1/d)))` panels, so 8 apiece at `panels = 50`
  on two variables.

- nodes:

  The number of Gauss-Legendre nodes per subinterval, default `12`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## Details

Numerical Gram Matrix of a Basis

Its accuracy is bounded by that of the derivatives it integrates, which
come from
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md):
a family with a closed-form derivative gets the accuracy of the
quadrature, and one with the finite-difference derivative carries that
error into the matrix.

All three shipped families and both wrappers register their own method,
so this one is reached only by a basis defined outside the package.

## See also

[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
which it calls;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures.
