# The Free Columns of an Operator's Null Space at New Values

Rebuilds what
[`operator_span()`](https://statmodels7.github.io/basis7/reference/operator_span.md)
produced, at new covariate values, from the operator and the scales
recorded there. Only the functions themselves are evaluated at `newx`.

## Usage

``` r
operator_span_apply(params, newx)
```

## Arguments

- params:

  The `params` element of an
  [`operator_span()`](https://statmodels7.github.io/basis7/reference/operator_span.md)
  result.

- newx:

  The new covariate values, a numeric vector.

## Value

A numeric matrix of `length(newx)` rows and one column per free
direction.
