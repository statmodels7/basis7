# Construct a B-Spline Basis

Returns a basis of `dimension` B-splines of the given degree on
\\\[\ell, u\]\\, with the interior knots placed at equal spacing. This
is the general-purpose choice: the functions have local support, so a
coefficient moves the fitted curve only near its own knots, and the
conditioning does not deteriorate as the dimension grows.

## Usage

``` r
bspline_basis(lower = 0, upper = 1, dimension = 5, degree = 3)
```

## Arguments

- lower, upper:

  The endpoints of the interval, each a single finite number with
  `lower < upper`. Default \\\[0, 1\]\\. They are the boundary knots,
  and evaluating outside them throws.

- dimension:

  The number of basis functions, a single whole number of at least
  `degree + 1`, default `5`. It is the number of columns, so it fixes
  the flexibility of the fit; the interior knot count follows as
  `dimension - degree - 1`.

- degree:

  The degree of the piecewise polynomials, a single non-negative whole
  number, default `3` for cubic splines. `0` gives indicators of the
  knot intervals and `1` piecewise linear functions. A spline of degree
  \\m\\ has \\m\\ non-trivial derivatives; above that
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  returns exact zeros.

## Value

An object of class
[BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md),
with `basis_name` `"bspline"`, `basis_params` holding `degree`, `knots`
and `boundary_knots`, and column names `bs1`, `bs2`, and so on.

## Dimension, degree and knots

A basis of \\K\\ functions of degree \\m\\ has \\K - m - 1\\ interior
knots, so \\K \ge m + 1\\; a smaller `dimension` throws, naming both
numbers. At equality there is no interior knot and the basis is the
polynomials of degree \\m\\ on the whole interval, which it spans
exactly. The knots are `seq(lower, upper, length.out = K - m + 1)` with
the endpoints dropped, so they are equally spaced; a quantile placement
is not offered, and a caller wanting one can build the class directly.

`degree = 0` gives indicator functions of the knot intervals, a step
basis, and `degree = 1` the piecewise linear hat functions.

## The basis is complete

All `dimension` functions are kept, so the rows of
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
sum to one and the basis spans the constant. Beside an intercept the
design is therefore rank deficient by one. Dropping a function is a
linear transformation of the basis: use
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
which keeps the object a basis, and leave the choice of constraint to
whatever owns the meaning of the term.

## Where each quantity comes from

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
and
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
call
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html),
which evaluates the recurrence and its exact derivative and integral.
The Gram matrix is integrated in this package, knot interval by knot
interval with a rule sized from the degree, and is exact to rounding.
Nothing about a B-spline basis is differenced, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports all three `FALSE`.

## References

de Boor, C. (2001). *A Practical Guide to Splines*, revised edition.
Springer.

## See also

[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
for a periodic basis and
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
for a global polynomial one;
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
to remove the constant;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
to rotate this basis into the Demmler-Reinsch form a penalized fit uses;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for its roughness penalty.

## Examples

``` r
b <- bspline_basis(dimension = 6)
b
#> Basis: bspline
#> Functions: 6   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree          3
#>   knots           0.3333, 0.6667
#>   boundary_knots  0, 1
#> Numerical: none

# Local support: each function is non-zero on a few knot intervals only.
round(basis_eval(b, c(0.1, 0.5, 0.9)), 3)
#>        bs1   bs2   bs3   bs4   bs5   bs6
#> [1,] 0.343 0.542 0.110 0.005 0.000 0.000
#> [2,] 0.000 0.031 0.469 0.469 0.031 0.000
#> [3,] 0.000 0.000 0.004 0.110 0.542 0.343

# The basis carries its own constant, so its rows sum to one.
max(abs(rowSums(basis_eval(b, seq(0, 1, length.out = 25))) - 1))
#> [1] 3.330669e-16

# The second-derivative Gram matrix, the matrix of a roughness penalty.
round(basis_gram(b, order = 2), 2)
#>        bs1     bs2     bs3     bs4     bs5    bs6
#> bs1  324.0 -445.50   94.50   27.00    0.00    0.0
#> bs2 -445.5  648.00 -182.25  -30.37   10.12    0.0
#> bs3   94.5 -182.25  121.50  -30.38  -30.37   27.0
#> bs4   27.0  -30.37  -30.38  121.50 -182.25   94.5
#> bs5    0.0   10.12  -30.37 -182.25  648.00 -445.5
#> bs6    0.0    0.00   27.00   94.50 -445.50  324.0

# At dimension = degree + 1 there is no interior knot, and the basis is
# the polynomials of that degree: a cubic is fitted exactly.
p <- bspline_basis(dimension = 4, degree = 3)
length(p@basis_params$knots)
#> [1] 0
x <- seq(0, 1, length.out = 40)
max(abs(lm.fit(basis_eval(p, x), x^3)$fitted.values - x^3))
#> [1] 1.745023e-16

# Degree 0 gives indicators of the knot intervals.
basis_eval(bspline_basis(dimension = 4, degree = 0), c(0.1, 0.3, 0.6, 0.9))
#>      bs1 bs2 bs3 bs4
#> [1,]   1   0   0   0
#> [2,]   0   1   0   0
#> [3,]   0   0   1   0
#> [4,]   0   0   0   1

# Too few functions for the degree is refused, with both numbers named.
try(bspline_basis(dimension = 3, degree = 3))
#> Error : 'dimension' (3) is too small for 'degree' (3): a B-spline basis of degree m needs at least m + 1 functions.
```
