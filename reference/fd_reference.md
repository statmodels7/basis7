# A Finite-Difference Reference and Its Uncertainty

Differentiates `f` once numerically and returns both the estimate and a
bound on its error, entry by entry.
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
uses the bound as the allowance of the comparison at each point, so a
point where the reference is unreliable does not report the basis as
wrong.

## Usage

``` r
fd_reference(f, x, lower, upper)
```

## Arguments

- f:

  A function of one numeric vector returning a numeric matrix with one
  row per point.
  [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  passes a closure over
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  or
  [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md).

- x:

  A numeric vector of evaluation points.

- lower, upper:

  The endpoints of the interval, so that the stencil can be shifted to
  one side near an end instead of leaving the domain.

## Value

A list of two matrices of the same shape as `f(x)`: `value`, the first
derivative estimated at the full step, and `uncertainty`, four times the
gap between that estimate and the one at half the step. An entry where
either estimate is `NA` gets an infinite uncertainty, so the comparison
there is unconstrained.

## Details

A central difference is valid only where the function has the
derivatives the stencil assumes. At a knot the third derivative of a
cubic spline jumps, so a stencil that straddles the knot returns a
number of the order of the jump and not of the truncation error, and
comparing an exact value against it would report a failure of the
reference as one of the basis.

The estimate is therefore recomputed with the step halved. At a smooth
point the two estimates differ by about three quarters of the truncation
error of the first, so four times the gap bounds that error; at a knot
the gap is large. The bound becomes the allowance of the comparison,
point by point. The two estimates are compared with each other and not
against a denominator floored at one, since near a kink both are small
and still differ by a factor.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its only caller;
[`rel_close()`](https://statmodels7.github.io/basis7/reference/rel_close.md),
which consumes `uncertainty`;
[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md),
which computes each estimate.
