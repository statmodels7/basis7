# The Order of a Differential Operator

The highest derivative the operator takes, which is the dimension of its
null space and the length of its weight vector.

## Usage

``` r
operator_order(op)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Value

A single integer.

## See also

[`operator_weights()`](https://statmodels7.github.io/basis7/reference/operator_weights.md),
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md).

## Examples

``` r
c(deriv = operator_order(deriv_operator(2)),
  harmonic = operator_order(harmonic_operator(365)),
  two_harmonics = operator_order(harmonic_operator(365, harmonics = 2)))
#>         deriv      harmonic two_harmonics 
#>             2             3             5 
```
