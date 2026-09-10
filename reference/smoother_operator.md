# The Operator a Smoother Penalizes With, Resolved

Returns `sm@order` with its period filled in from the smoother's
interval where it was left `NULL`, which is what makes
`fourier_smooth(lower = 0, upper = 365)` penalize on a cycle of 365
without the period being written twice.

## Usage

``` r
smoother_operator(sm, x)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- x:

  The covariate, from which the interval is taken where the smoother
  does not fix it.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
with numeric weights.

## See also

[`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md),
which it calls.
