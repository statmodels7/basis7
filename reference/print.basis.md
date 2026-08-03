# Print a Basis

Reports the family, the number of functions, the interval, any
parameters the family carries, and which of the derived quantities are
computed numerically.

## Arguments

- x:

  An object inheriting from class `basis`.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
bspline_basis(dimension = 6)
#> Basis: bspline
#> Functions: 6   Interval: [0, 1]
#> Parameters:
#>   degree          3
#>   knots           0.3333, 0.6667
#>   boundary_knots  0, 1
#> Numerical: none
fourier_basis(dimension = 5)
#> Basis: fourier
#> Functions: 5   Interval: [0, 1]
#> Parameters:
#>   omega        1
#>   n_pairs      2
#>   full_period  TRUE
#> Numerical: none
```
