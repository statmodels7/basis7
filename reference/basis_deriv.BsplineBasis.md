# Derivatives of a B-Spline Basis

Returns the `order`-th derivative of every B-spline, exactly, from the
derivative form of the Cox-de Boor recurrence in
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html).
An order above the degree short-circuits to the zero matrix, which is
the value of that derivative; nothing is thrown.

## Arguments

- basis:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a single non-negative whole number, default `1`.
  Above `basis@basis_params$degree` the result is exactly zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names `bs1`, `bs2`, and so on.

## Details

A spline of degree \\m\\ is piecewise polynomial of that degree, so it
has \\m\\ non-trivial derivatives and the rest vanish: a cubic gives
three, and `order = 4` is exactly zero everywhere. The short-circuit
happens here because splines2 rejects a `derivs` above the degree
instead of returning zeros.

The derivative of order \\m\\ is a step function, discontinuous at each
interior knot, and the value returned at a knot is the one the
recurrence gives there. That matters for a Gram matrix, which is why
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md)
integrates knot interval by knot interval and never across one.

## See also

[`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md),
the one call into splines2;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic.
