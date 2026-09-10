# Fill In an Operator's Period From an Interval

Returns the operator with its period set to the width of
`[lower, upper]` where it was left `NULL`, and unchanged where it was
given. It is what a smoother calls before building its penalty, so that
`fourier_smooth(lower = 0, upper = 365)` penalizes on a cycle of 365
without the period having to be written twice.

## Usage

``` r
operator_resolve(op, lower, upper)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

- lower, upper:

  The interval, two numbers.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
with numeric weights.

## See also

[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
whose `period = NULL` this resolves.

## Examples

``` r
operator_weights(harmonic_operator())
#> [1] NA NA NA
round(operator_weights(operator_resolve(harmonic_operator(), 0, 365)), 8)
#> [1] 0.00000000 0.00029633 0.00000000

# an operator given a period keeps it, whatever the interval is
round(operator_weights(operator_resolve(harmonic_operator(12), 0, 365)), 8)
#> [1] 0.0000000 0.2741557 0.0000000
```
