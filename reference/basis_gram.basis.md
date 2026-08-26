# Numerical Gram Matrix of a Basis

The inner-product method every basis inherits from the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
composite Gauss-Legendre over `panels` equal subintervals of the
interval, applied to the requested derivatives. A one-line wrapper over
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
which is where the work is and which the other families also call when
their closed form does not apply.

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
  `ceiling(panels^(1/d))` panels, so 8 apiece at `panels = 50` on two
  variables.

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

Its accuracy is bounded by the derivative it integrates. On a polynomial
family the order-0 matrix agrees with the closed form to 5.2e-15, while
at order 2 the gap is 3.0e-08 relative, which is the finite-difference
error of
[`basis_deriv.basis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.basis.md)
carried through the integral, and no fault of the quadrature.

All three shipped families and both wrappers register their own method,
so this one is reached only by a basis defined outside the package.

## See also

[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
which it calls;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures.
