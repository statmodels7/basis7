# Gram Matrix of a B-Spline Basis

Returns the inner products of the `order`-th derivatives exactly, by
integrating over one knot interval at a time with a Gauss-Legendre rule
sized so that it reproduces the integrand exactly. This is the matrix a
roughness penalty on a spline is built from.

## Arguments

- basis:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- order:

  The derivative order, a single non-negative whole number, default `0`.
  `2` is the usual roughness penalty for a cubic.

- at, weight:

  Handled in the body of
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  before dispatch, so they never arrive here. Named only because S7
  requires a method's formals to contain the generic's.

- ...:

  Unused, and accepted so that the signature matches the generic's.
  `panels` or `nodes` signals an error, the matrix being exact.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
column names `bs1`, `bs2`, and so on. Banded, entry \\(a, b)\\ being
zero whenever the two supports do not overlap.

## Exactness

On one knot interval the order-\\d\\ derivative of a spline of degree
\\m\\ is a polynomial of degree \\m - d\\, so the integrand \\B_a^{(d)}
B_b^{(d)}\\ has degree \\2(m - d)\\. A Gauss-Legendre rule with \\m -
d + 1\\ nodes is exact to degree \\2(m - d) + 1\\, which is one higher,
so the only error left is rounding.

The breaks are the boundary knots and the interior knots, so no panel
straddles a knot, and the exactness depends on that: at `order = m` the
derivative is a step function, which a rule spanning a knot does not
integrate exactly. The `weight` route of
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
uses equal panels that do not line up with the knots, and is therefore
approximate for a spline.

## Above the degree

An `order` above `degree` returns the zero matrix, every derivative
having vanished. Below it the matrix is singular with an
`order`-dimensional null space, the polynomials of lower degree
differentiating away.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic and the alternative measures;
[`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md),
which builds the composite rule.
