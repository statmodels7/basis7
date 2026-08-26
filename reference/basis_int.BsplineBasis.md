# Integral of a B-Spline Basis

Returns \\\int\_{\ell}^{x} B_j(t)\\\mathrm{d}t\\ for every function,
exactly, from
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html)
with `integral = TRUE`. That function anchors its integral at the lower
boundary knot, which is the same convention
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
states, so no correction is applied here.

## Arguments

- basis:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
exactly zero in the row at `basis@lower`.

## Details

The integral of a spline of degree \\m\\ is a spline of degree \\m + 1\\
on the same knots, and splines2 evaluates it from the recurrence, with
no quadrature anywhere.

The row at `basis@upper` holds the area under each function. Because the
basis is a partition of unity, those areas sum to the width of the
interval, which is a cheap check that the two conventions agree.

## See also

[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention;
[`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md),
the one call into splines2.
