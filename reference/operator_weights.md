# The Weights of a Differential Operator

The coefficients \\w_0, \ldots, w\_{m-1}\\ of \\L x = \sum_j w_j D^j x +
D^m x\\, in increasing order of differentiation. The leading coefficient
is 1 and is not among them.

## Usage

``` r
operator_weights(op)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Value

A numeric vector of `operator_order(op)` entries, or `NA` of that length
where the operator's period has not been resolved.

## See also

[`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md),
which fills in a period.

## Examples

``` r
round(operator_weights(harmonic_operator(365)), 8)
#> [1] 0.00000000 0.00029633 0.00000000
operator_weights(deriv_operator(3))
#> [1] 0 0 0
```
