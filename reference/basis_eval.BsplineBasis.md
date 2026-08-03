# Evaluate a B-Spline Basis

The design matrix of the B-spline recurrence, from
[`bSpline`](https://wwenjie.org/splines2/reference/bSpline.html).

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
