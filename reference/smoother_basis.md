# The Basis a Smoother Builds On

Returns the
[basis](https://statmodels7.github.io/basis7/reference/basis.md) a
smoother expands the covariate in, before any constraint or
reparametrization. It is the one step of the construction that differs
between families, so a family is added by registering a method here
rather than by rewriting
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md).

## Usage

``` r
smoother_basis(sm, x, ...)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- x:

  The covariate, a numeric vector. A family whose interval is not fixed
  on the object reads it from these values.

- ...:

  Passed to methods.

## Value

A [basis](https://statmodels7.github.io/basis7/reference/basis.md) of
`sm@dimension` functions.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md),
which calls it.

## Examples

``` r
set.seed(1)
x <- runif(50)
b <- smoother_basis(bspline_smooth(k = 8), x)
c(dimension = b@dimension, lower = b@lower, upper = b@upper)
#>  dimension      lower      upper 
#> 8.00000000 0.01241182 0.99288461 

# The interval is the range of the data padded by a thousandth of its
# width, unless the smoother fixes it.
bf <- smoother_basis(bspline_smooth(k = 8, lower = 0, upper = 1), x)
c(lower = bf@lower, upper = bf@upper)
#> lower upper 
#>     0     1 
```
