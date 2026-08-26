# Compare Two Matrices Relative to Their Own Magnitude

Reports whether two matrices agree to a relative tolerance, with the
denominator taken from the values themselves and floored at a millionth
of the column's own scale. The one comparison
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
makes, so that its derivative, integral and Gram checks all read
agreement the same way.

## Usage

``` r
rel_close(a, b, tol, slack = NULL)
```

## Arguments

- a, b:

  Numeric matrices of the same shape. Neither is privileged; the
  comparison is symmetric.

- tol:

  The relative tolerance, a single positive number.

- slack:

  An optional numeric matrix of the same shape as `a`, an absolute
  allowance added entry by entry through `pmax(tol * den, slack)`.
  [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  passes
  [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)'s
  `uncertainty` here, so a point whose finite-difference reference is
  unreliable is allowed the error that reference has. `NULL`, the
  default, allows none.

## Value

A single `TRUE` or `FALSE`: `TRUE` when every informative entry agrees
within its own allowance.

## Why the denominator is not floored at one

Flooring at one would flatten a disagreement between two small numbers
into apparent agreement, and small numbers are most of what a basis
produces: a B-spline is exactly zero on most of its interval. The
denominator is `pmax(abs(a), abs(b))`.

## Why it is floored at all

A basis function's derivative crosses zero, and at the crossing the
pointwise value vanishes while the numerical reference carries its usual
rounding error. Dividing that error by nothing reports a failure of the
reference as a failure of the basis. The floor is `1e-6` times the
largest absolute value in the same column, so a proportional error stays
detectable wherever the curve is large, which is where a wrong formula
shows itself.

## Columns with nothing to say

A column whose whole scale is below `1e-8` of the largest column's is
skipped: neither side carries information there. If every column is
skipped the answer is `TRUE`.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its caller, and
[`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md),
which supplies `slack`.
