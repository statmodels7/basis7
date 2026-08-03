# Numerical Integral of a Basis

The default integration method: composite Gauss-Legendre from the lower
endpoint, accumulated over the sorted evaluation points so that the
whole set costs one pass rather than one quadrature each.

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points.

- nodes:

  The number of quadrature nodes per segment.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

Numerical Integral of a Basis

The rule is placed on the segments between consecutive evaluation
points, and the results are accumulated, which makes the value at the
lower endpoint exactly zero by construction rather than by cancellation.
