# Gram Matrix of a Legendre Basis

Returns the inner products of the `order`-th derivatives exactly, with
no quadrature error at any order. At `order = 0` the matrix is diagonal
and written down; above that a Gauss-Legendre rule with as many nodes as
the basis has functions integrates the product exactly, the integrand
being a polynomial of low enough degree.

## Arguments

- basis:

  A
  [PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- order:

  The derivative order, a single non-negative whole number, default `0`.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
column names `P0`, `P1`, and so on. Diagonal at `order = 0`, singular
with an `order`-dimensional null space above it, and exactly zero above
`dimension - 1`.

## Order zero

The polynomials satisfy \\\int\_{-1}^{1} P_m P_n \\\mathrm{d}t =
2\delta\_{mn}/(2n+1)\\, so on \\\[\ell, u\]\\ the matrix is diagonal
with entries \\(u - \ell)/(2n+1)\\. It is written directly and costs
nothing.

## Higher orders

Derivatives of Legendre polynomials are not orthogonal, so the matrix is
full. The integrand \\P_m^{(d)} P_n^{(d)}\\ is a polynomial of degree at
most \\2(K - 1 - d)\\, and a \\K\\-node Gauss-Legendre rule is exact to
degree \\2K - 1\\, which is larger for every \\d \ge 1\\. The quadrature
is therefore exact to rounding, and the result is symmetrized as
`(G + t(G))/2` because the two triangles of a crossproduct differ in
their last bits.

An `order` above `dimension - 1` returns the zero matrix, every
derivative having vanished.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures.
