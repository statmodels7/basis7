# Gram Matrix of a B-Spline Basis

Returns the inner products of the `order`-th derivatives exactly, by
integrating over one knot interval at a time with a Gauss-Legendre rule
sized so that it reproduces the integrand exactly. This is the matrix a
roughness penalty on a spline is built from, so its exactness is worth
the small amount of work.

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

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
column names `bs1`, `bs2`, and so on. Banded, entry \\(a, b)\\ being
zero whenever the two supports do not overlap.

## Why it is exact

On one knot interval the order-\\d\\ derivative of a spline of degree
\\m\\ is a polynomial of degree \\m - d\\, so the integrand \\B_a^{(d)}
B_b^{(d)}\\ has degree \\2(m - d)\\. A Gauss-Legendre rule with \\m -
d + 1\\ nodes is exact to degree \\2(m - d) + 1\\, which is one higher,
so the only error left is floating point. Against the same knot-aligned
construction run at 20 nodes instead, the worst entry agrees to 3.1e-16,
2.1e-14, 8.0e-13 and 7.3e-12 at orders 0 to 3 on a cubic basis of six
functions.

The breaks are the boundary knots and the interior knots, so no panel
straddles a knot, and the exactness claim rests on that: at `order = m`
the derivative is a step function, and a rule spanning a knot would
integrate the wrong thing. Measured on the same basis, the general
`weight` route of
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
whose 50 equal panels do not line up with the knots, is out by 1.4e-3 at
order 2 and by 65 at order 3.

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
