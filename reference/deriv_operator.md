# The Derivative Operator

The operator \\L x = D^m x\\, whose penalty \\\int (D^m x)^2\\ is the
integrated squared derivative every smoother in this package penalized
with before operators existed. Writing `order = m` on a smoother is the
shorthand for `order = deriv_operator(m)` and builds the identical
construction.

## Usage

``` r
deriv_operator(m = 2)
```

## Arguments

- m:

  The order of the derivative, a whole number of at least 1.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Details

Its null space is the polynomials of degree below \\m\\: the
characteristic polynomial is \\r^m\\, a root at zero of multiplicity
\\m\\, which contributes \\1, t, \ldots, t^{m-1}\\. A strongly penalized
fit therefore contracts to a constant at `m = 1`, to a straight line at
`m = 2` and to a parabola at `m = 3`.

## See also

[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
for the periodic one,
[`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md)
for the general form, and
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for what a maximal penalty leaves.

## Examples

``` r
deriv_operator(2)
#> <linear differential operator>  order 2
#>   L x = D^2 x
#>   null space: 1, t

# what a strongly penalized fit contracts to
operator_null(deriv_operator(3))
#>   label rate freq degree part
#> 1     1    0    0      0     
#> 2     t    0    0      1     
#> 3   t^2    0    0      2     

# `order = 2` on a smoother is this operator, and gives the same block
set.seed(1)
x <- sort(runif(100))
a <- smoother_build(bspline_smooth(k = 8, order = 2), x)
b <- smoother_build(bspline_smooth(k = 8, order = deriv_operator(2)), x)
identical(a$X, b$X) && identical(a$S, b$S)
#> [1] TRUE
```
