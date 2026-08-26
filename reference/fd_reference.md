# A Finite-Difference Reference, and Where It Can Be Trusted

Differentiates `f` once numerically and returns both the estimate and a
bound on its own error, entry by entry.
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
uses the second to decide how much slack the comparison at each point
deserves, so a point where the reference is unreliable does not report
the basis as wrong.

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

A central difference is only valid where the function has the
derivatives the stencil assumes. A spline does not: at a knot its third
derivative jumps, so a stencil that straddles the knot returns a number
of the order of the jump and not of the truncation error, and comparing
an exact analytical value against it reports a failure of the
*reference*.

Recomputing with the step halved says how much of the reference is
error. For a smooth point the two differ by about three quarters of the
truncation, so the gap between them bounds the reference's own
uncertainty; at a knot it is large. Nothing is discarded: the gap
becomes the slack allowed to the comparison, so each point contributes
exactly the accuracy its reference supports. This is the same device
used elsewhere in the toolkit for a parameter that is not
differentiable, and it needs the same care: the two estimates are
compared with each other, not against a denominator floored at one,
since near a kink both are small and still differ by a factor.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its only caller;
[`rel_close()`](https://statmodels7.github.io/basis7/reference/rel_close.md),
which consumes `uncertainty`;
[`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md),
which computes each estimate.
