# Differentiate a Basis

Returns the `order`-th derivative of every basis function at the given
points, as a matrix of the same shape as
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md).

## Usage

``` r
basis_deriv(basis, x, order = 1L, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a non-negative integer.

- ...:

  Passed to methods.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

Derivative order is an argument rather than a family of generics because
it is unbounded: a Fourier basis is differentiable to any order, and a
spline of degree \\k\\ has \\k\\ non-trivial derivatives followed by
zeros. An order beyond what the family supports returns the zero matrix
rather than raising, which is the value of the derivative and not an
omission.

`order = 0` returns
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
so that a loop over orders needs no special case.

A subclass that registers no method for this generic gets the numerical
one of the
[`basis`](https://statmodels7.github.io/basis7/reference/basis.md)
class, which applies a single central difference stencil to
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md).

## See also

[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_is_numerical`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)

## Examples

``` r
b <- bspline_basis(dimension = 6, degree = 3)
round(basis_deriv(b, c(0.25, 0.5, 0.75), order = 1), 3)
#>         bs1    bs2    bs3    bs4   bs5   bs6
#> [1,] -0.562 -2.391  2.109  0.844 0.000 0.000
#> [2,]  0.000 -0.562 -1.687  1.687 0.562 0.000
#> [3,]  0.000  0.000 -0.844 -2.109 2.391 0.563

# a cubic spline has no fourth derivative
all(basis_deriv(b, 0.5, order = 4) == 0)
#> [1] TRUE
```
