# Plot a Basis

Draws every basis function over the interval, or its derivative or
integral.

## Arguments

- x:

  An object inheriting from class `basis`.

- order:

  What to draw: `0` for the basis functions, a positive integer for that
  derivative, `-1` for the integral from the lower endpoint.

- n:

  The number of points at which to evaluate.

- ...:

  Passed to [`matplot`](https://rdrr.io/r/graphics/matplot.html).

## Value

`x`, invisibly.

## Examples

``` r
b <- bspline_basis(dimension = 6)
plot(b)

plot(b, order = 1)

plot(b, order = -1)
```
