# The Null-Space Functions Evaluated

The design matrix of the functions
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
lists, evaluated at `x`: one column per function, in the order that
table gives them.

## Usage

``` r
operator_null_design(op, x)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

- x:

  The points to evaluate at, a numeric vector.

## Value

A numeric matrix of `length(x)` rows and `operator_order(op)` columns,
with the labels of
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
as column names.

## Details

These are the mathematical functions, unscaled. A smoother that restores
them as free columns scales them first and records the scale, so that a
prediction reapplies the same columns instead of recomputing them from
new data; see
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md).

A rate far from zero over a wide interval makes \\e^{at}\\ overflow, so
the design matrix of an operator with a large real root cannot be
formed. The operators built by
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
have only roots on the imaginary axis, zero included, and this
difficulty does not arise.

## See also

[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for what the columns are.

## Examples

``` r
x <- seq(0, 365, length.out = 5)
round(operator_null_design(harmonic_operator(365), x), 4)
#>      1 sin(0.0172142 t) cos(0.0172142 t)
#> [1,] 1                0                1
#> [2,] 1                1                0
#> [3,] 1                0               -1
#> [4,] 1               -1                0
#> [5,] 1                0                1

# the constant and t, the null space of the second derivative
operator_null_design(deriv_operator(2), 1:5)
#>      1 t
#> [1,] 1 1
#> [2,] 1 2
#> [3,] 1 3
#> [4,] 1 4
#> [5,] 1 5
```
