# Integral of a B-Spline Basis

The exact integral from the lower boundary knot, from
[`bSpline`](https://wwenjie.org/splines2/reference/bSpline.html), which
follows the same convention as this package: the value at the lower
endpoint is zero.

## Arguments

- basis:

  A
  [`BsplineBasis`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.
