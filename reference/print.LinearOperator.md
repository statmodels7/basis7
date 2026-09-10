# How a Differential Operator Prints

Three lines: the order, the operator written out, and the functions its
null space holds.

## Arguments

- x:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
harmonic_operator(365)
#> <linear differential operator>  order 3
#>   L x = 0.0002963 Dx + D^3 x
#>   null space: 1, sin(0.0172142 t), cos(0.0172142 t)
deriv_operator(2) * harmonic_operator(365)
#> <linear differential operator>  order 5
#>   L x = 0.0002963 D^3 x + D^5 x
#>   null space: 1, t, t^2, sin(0.0172142 t), cos(0.0172142 t)
```
