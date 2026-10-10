# Gram Matrix of a Transformed Basis

Returns the congruence \\T^\top G\\T\\ of the parent's Gram matrix. A
parent whose inner products are exact passes that exactness on, so the
Gram matrix of an orthonormalized B-spline is the identity up to
rounding and no quadrature is run.

## Arguments

- basis:

  A
  [TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- order:

  The derivative order, a single non-negative whole number, default `0`,
  passed to the parent unchanged.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- ...:

  Passed to the parent's method, so `panels` and `nodes` reach a parent
  whose Gram matrix is a quadrature.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## Details

The congruence preserves symmetry and positive semidefiniteness, and its
rank is at most the smaller of the rank of \\G\\ and `ncol(T)`. The
result is symmetrized as `(G + t(G))/2` before it is returned.

[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
relies on this: it chooses \\T\\ so that \\T^\top G\\T\\ is the
identity.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic;
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
which chooses \\T\\ to make this the identity.
