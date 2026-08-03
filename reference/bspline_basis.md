# Construct a B-Spline Basis

A basis of B-splines of the given degree on \\\[\ell, u\]\\, with
interior knots placed at equal spacing.

## Usage

``` r
bspline_basis(lower = 0, upper = 1, dimension = 5, degree = 3)
```

## Arguments

- lower, upper:

  The endpoints of the interval.

- dimension:

  The number of basis functions, at least `degree + 1`.

- degree:

  The degree of the piecewise polynomials. Three, the default, gives
  cubic splines.

## Value

An object of class
[`BsplineBasis`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md).

## Details

The basis is complete: all `dimension` functions are kept, so the rows
of
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
sum to one. Dropping a function for identifiability is a linear
transformation of the basis and a decision for the layer that knows what
the term means.

A basis of \\K\\ functions with degree \\m\\ has \\K - m - 1\\ interior
knots, so \\K \ge m + 1\\ is required, with equality giving the
polynomials of degree \\m\\ on the whole interval and no interior knot
at all.

## References

de Boor, C. (2001). *A Practical Guide to Splines*. Springer.

## See also

[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
[`check_basis`](https://statmodels7.github.io/basis7/reference/check_basis.md)

## Examples

``` r
b <- bspline_basis(dimension = 6)
b
#> Basis: bspline
#> Functions: 6   Interval: [0, 1]
#> Parameters:
#>   degree          3
#>   knots           0.3333, 0.6667
#>   boundary_knots  0, 1
#> Numerical: none

# local support: each function is non-zero on a few knot intervals only
round(basis_eval(b, c(0.1, 0.5, 0.9)), 3)
#>        bs1   bs2   bs3   bs4   bs5   bs6
#> [1,] 0.343 0.542 0.110 0.005 0.000 0.000
#> [2,] 0.000 0.031 0.469 0.469 0.031 0.000
#> [3,] 0.000 0.000 0.004 0.110 0.542 0.343

# the second-derivative Gram matrix, which a roughness penalty integrates
round(basis_gram(b, order = 2), 2)
#>        bs1     bs2     bs3     bs4     bs5    bs6
#> bs1  324.0 -445.50   94.50   27.00    0.00    0.0
#> bs2 -445.5  648.00 -182.25  -30.37   10.12    0.0
#> bs3   94.5 -182.25  121.50  -30.38  -30.37   27.0
#> bs4   27.0  -30.37  -30.38  121.50 -182.25   94.5
#> bs5    0.0   10.12  -30.37 -182.25  648.00 -445.5
#> bs6    0.0    0.00   27.00   94.50 -445.50  324.0
```
