# The Eilers-Marx Coordinates of a Clamped B-Spline Basis

Returns the matrix \\M\\ with \\a = M c\\, where \\c\\ are the
coefficients of a function on the clamped B-spline basis `b` and \\a\\
are the coefficients of the same function on the basis of Eilers and
Marx (1996), whose knots are equally spaced with the same step beyond
the interval as inside it.

## Usage

``` r
pspline_map(b)
```

## Arguments

- b:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  built by
  [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
  with equally spaced interior knots.

## Value

A square numeric matrix of `b@dimension` rows and columns.

## Details

The two bases have the same degree and the same breakpoints inside the
interval, so on the interval they span the same splines and \\M\\ is
square and invertible. It is computed by evaluating both bases on a grid
of ten points per basis function and solving the least-squares system,
which is exact to rounding. A difference penalty \\\lVert D a\rVert^2\\
then reads \\c^\top M^\top D^\top D M c\\ on the clamped coefficients.
On the Eilers-Marx knots the Greville abscissae are equally spaced, so
coefficients that are a polynomial of degree below `diff` in their index
give exactly a polynomial of that degree.

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11(2), 89-121.

## Examples

``` r
b <- bspline_basis(0, 1, dimension = 8)
M <- basis7:::pspline_map(b)
dim(M)
#> [1] 8 8
```
