# A Finite-Difference Reference, and Where It Can Be Trusted

Differentiates `f` numerically, and reports at which points the result
is worth comparing anything against.

## Usage

``` r
fd_reference(f, x, lower, upper)
```

## Arguments

- f:

  A function of a numeric vector returning a matrix.

- x:

  A numeric vector of evaluation points.

- lower, upper:

  The endpoints of the interval.

## Value

A list with the reference `value` and the matrix `uncertainty` bounding
its error at each entry.

## Details

A central difference is only valid where the function has the
derivatives the stencil assumes. A spline does not: at a knot its third
derivative jumps, so a stencil that straddles the knot returns a number
of the order of the jump rather than of the truncation error, and
comparing an exact analytical value against it reports a failure of the
*reference*.

Recomputing with the step halved says how much of the reference is
error. For a smooth point the two differ by about three quarters of the
truncation, so the gap between them bounds the reference's own
uncertainty; at a knot it is large, and says so. Nothing is discarded:
the gap becomes the slack the comparison is allowed, so a point tells
the check as much as it honestly can and no more. This is the same
device used elsewhere in the toolkit for a parameter that is not
differentiable, and it needs the same care: the two estimates are
compared with each other, not against a denominator floored at one,
since near a kink both are small and still differ by a factor.
