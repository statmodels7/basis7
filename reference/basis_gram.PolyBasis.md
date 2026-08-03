# Gram Matrix of a Legendre Basis

Exact inner products: diagonal in closed form at order zero, and by an
exact Gauss-Legendre rule above it.

## Arguments

- basis:

  A
  [`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- order:

  The derivative order.

- at, weight:

  Handled by the generic before dispatch; unused here.

- ...:

  Unused.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

The Legendre polynomials satisfy \\\int\_{-1}^{1} P_m P_n \\ \mathrm{d}t
= 2\delta\_{mn}/(2n+1)\\, so on the basis interval the order-zero matrix
is diagonal with entries \\(u - \ell)/(2n+1)\\. At higher orders the
derivatives are no longer orthogonal, but the integrand stays a
polynomial of degree at most \\2(K - 1 - d)\\, so a rule with \\K\\
nodes integrates it exactly.
