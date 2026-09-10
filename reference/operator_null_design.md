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
them as free columns centers and scales them first and records what it
did, so that a prediction reapplies the same columns rather than
recomputing them from new data;
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
is where that happens.

A rate far from zero over a wide interval overflows: \\e^{at}\\ is what
it is, and an operator with a large real root is a statement about
growth that a design matrix cannot hold. The two operators this package
builds for itself,
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md)
and
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
have purely imaginary roots and no such difficulty.

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

# L applied to its own null space is zero, which is what makes it the
# null space; here to the accuracy of a central difference
round(colSums(abs(operator_null_design(deriv_operator(2), 1:5))), 6)
#>  1  t 
#>  5 15 
```
