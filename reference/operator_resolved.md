# Has This Operator's Period Been Resolved?

Reports whether the weights are numbers rather than the `NA` placeholder
that
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
leave when they are given no period.

## Usage

``` r
operator_resolved(op)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Value

A single logical.
