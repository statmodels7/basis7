# The Linear Differential Operator Class

A linear differential operator with constant coefficients, \$\$L x = w_0
x + w_1 Dx + \cdots + w\_{m-1} D^{m-1} x + D^m x,\$\$ which specifies
what the penalty of a smoother measures. The leading coefficient is 1 by
construction, so an operator is determined by the weights \\w_0, \ldots,
w\_{m-1}\\ and by nothing else.

## Usage

``` r
LinearOperator(
  weights = integer(0),
  operator_name = character(0),
  operator_params = list()
)
```

## Arguments

- weights:

  The coefficients \\w_0, \ldots, w\_{m-1}\\, a numeric vector, or `NA`
  of the right length for an operator whose period has not been resolved
  yet.

- operator_name:

  A short name used when the operator prints.

- operator_params:

  A list of whatever the constructor recorded, such as the period and
  the number of harmonics.

## Value

An S7 object of class `LinearOperator`.

## Details

An operator is built with
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
or
[`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md),
and two are composed with `*`. The properties are read through the
accessors
[`operator_order()`](https://statmodels7.github.io/basis7/reference/operator_order.md),
[`operator_weights()`](https://statmodels7.github.io/basis7/reference/operator_weights.md)
and
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md).
For an operator whose period is still to be resolved, the weights are
`NA` and
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
signals an error.

## See also

[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
and
[`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md),
which build one;
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for the space a maximal penalty leaves untouched.

## Examples

``` r
op <- harmonic_operator(365)
c(class = class(op)[1], order = operator_order(op))
#>                    class                    order 
#> "basis7::LinearOperator"                      "3" 
```
