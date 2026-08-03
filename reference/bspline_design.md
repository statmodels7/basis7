# Call splines2 for a B-Spline Design Matrix

The single point at which this package talks to splines2, so that the
knot arguments are assembled once and the dependency stays behind the S7
interface.

## Usage

``` r
bspline_design(basis, x, derivs = 0L, integral = FALSE)
```

## Arguments

- basis:

  A
  [`BsplineBasis`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- derivs:

  The derivative order, or zero.

- integral:

  Whether to return the integral instead.

## Value

A numeric matrix, stripped of the attributes splines2 attaches.
