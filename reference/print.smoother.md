# Print a Smoother

Prints the family, the number of basis functions and the four decisions
the smoother carries.

## Arguments

- x:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- ...:

  Ignored.

## Value

`x`, invisibly. Called for the printing.

## Examples

``` r
bspline_smooth(k = 12)
#> <basis7::BsplineSmoother> 12 functions, degree 3
#>   penalty: derivative of order 2, lebesgue measure
#>   null space: keep     coordinates: dr
#>   interval: from the data
bspline_smooth(k = 12, null_space = "drop")
#> <basis7::BsplineSmoother> 12 functions, degree 3
#>   penalty: derivative of order 2, lebesgue measure
#>   null space: drop     coordinates: dr
#>   interval: from the data
```
