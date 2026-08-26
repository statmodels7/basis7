# Call splines2 for a B-Spline Design Matrix

The single point at which this package talks to splines2. It assembles
the knot arguments from `basis@basis_params`, calls
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html)
once, and returns a plain matrix, so the dependency stays behind the S7
interface and no caller has to know its argument names or its return
class.

## Usage

``` r
bspline_design(basis, x, derivs = 0L, integral = FALSE)
```

## Arguments

- basis:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points, already validated by the
  generic.

- derivs:

  The derivative order, default `0`. Must not exceed
  `basis@basis_params$degree`;
  [`basis_deriv.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.BsplineBasis.md)
  short-circuits above that and never calls here.

- integral:

  `TRUE` to return the integral anchored at the lower boundary knot
  instead of the functions. Default `FALSE`.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns, of
class `matrix` alone, with no dimnames and none of the attributes
splines2 attaches. Callers add the column names through
[`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md).

## Details

`intercept = TRUE` is passed always, so all `dimension` functions come
back; splines2 would otherwise drop the first.

The returned object is of class `BSpline` and carries ten attributes,
among them `x`, `knots`, `degree` and `intercept`. Rebuilding it as
`matrix(as.numeric(out), ...)` strips every one, which matters because
those attributes would survive arithmetic and reappear on a matrix that
no longer describes them. The dimensions are taken from `length(x)` and
`basis@dimension`, so a mismatch with what splines2 returned surfaces
here.

`derivs` and `integral` are mutually exclusive in practice, each caller
setting at most one.

## See also

[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html),
the function called;
[`basis_eval.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.BsplineBasis.md),
[`basis_deriv.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.BsplineBasis.md)
and
[`basis_int.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.BsplineBasis.md),
its three callers.
