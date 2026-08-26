# Gram Matrix of a Fourier Basis

Returns the inner products of the `order`-th derivatives. On a whole
period the matrix is diagonal and written down in closed form; on any
other period the orthogonality fails and the method computes a
quadrature instead.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- order:

  The derivative order, a single non-negative whole number, default `0`.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- ...:

  Passed to
  [`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
  on the non-full-period branch, where `panels` and `nodes` control the
  quadrature. Ignored on the closed branch.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
column names `const`, `sin1`, and so on. Diagonal when the interval is
one whole period, full otherwise; singular for any `order >= 1`, the
constant differentiating away.

## The closed form

Over a whole period the trigonometric functions are mutually orthogonal,
and stay orthogonal after differentiation, a derivative only shifting
the phase and rescaling. The order-\\d\\ matrix is therefore diagonal,
with

\$\$G\_{00} = \omega \\ (d = 0), \quad 0 \\ (d \ge 1), \qquad G\_{jj} =
\frac{\omega}{2}\left(\frac{2\pi j}{\omega}\right)^{2d}\$\$

for both members of pair \\j\\. At `dimension = 7` on \\\[0, 1\]\\ the
order-2 diagonal is
`0, 779.27, 779.27, 12468.4, 12468.4, 63121.1, 63121.1`, matching
\\(2\pi j)^4/2\\ exactly.

## The other period

When `basis_params$full_period` is `FALSE` the method delegates to
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
composite Gauss-Legendre over the interval, and the result is a full
matrix.

[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports `basis_gram` as `TRUE` on such a basis, the family answering
through
[`basis_numerical_route.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.FourierBasis.md)
rather than through the class the method is registered on, which is
`FourierBasis` in both branches.
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
reads the same predicate, so it holds this matrix to the tolerance a
quadrature deserves.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures;
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
for the branch taken on another period.
