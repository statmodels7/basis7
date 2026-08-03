# Map a Quadrature Rule onto Intervals

Places an `n`-point Gauss-Legendre rule on each interval given by
consecutive breakpoints, and returns the pooled nodes and weights.

## Usage

``` r
quad_rule(breaks, n)
```

## Arguments

- breaks:

  A numeric vector of at least two increasing breakpoints.

- n:

  The number of nodes per interval.

## Value

A list with components `nodes` and `weights`.
