# Evaluate a B-Spline Basis

Returns the B-spline design matrix at the given points, evaluated by
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html)
from the Cox-de Boor recurrence. At most `degree + 1` entries of any row
are non-zero, and the row sums are one up to rounding.

## Arguments

- basis:

  A
  [BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names `bs1`, `bs2`, and so on.

## Details

The call goes through
[`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md)
with `intercept = TRUE`, so all `dimension` functions are returned and
none is dropped for identifiability. The result is stripped of the
attributes and of the `BSpline` class that splines2 attaches, leaving a
plain matrix.

## See also

[`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md),
the one call into splines2;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
for the generic.
