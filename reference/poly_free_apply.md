# The Free Columns of a Polynomial Null Space at New Values

Rebuilds what
[`poly_free()`](https://statmodels7.github.io/basis7/reference/poly_free.md)
produced, at new covariate values, from the centering, the
orthogonalization coefficients and the scales recorded there. Nothing is
recomputed from `newx`.

## Usage

``` r
poly_free_apply(params, newx)
```

## Arguments

- params:

  The `params` element of a
  [`poly_free()`](https://statmodels7.github.io/basis7/reference/poly_free.md)
  result.

- newx:

  The new covariate values, a numeric vector.

## Value

A numeric matrix of `length(newx)` rows and one column per step.
