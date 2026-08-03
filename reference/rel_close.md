# Compare Two Matrices Relative to Their Own Magnitude

Whether two matrices agree to a relative tolerance, with the denominator
taken from the values themselves rather than floored at one.

## Usage

``` r
rel_close(a, b, tol)
```

## Arguments

- a, b:

  Numeric matrices of the same shape.

- tol:

  The relative tolerance.

## Value

`TRUE` or `FALSE`.

## Details

Flooring the denominator at one would flatten a disagreement between two
small numbers into apparent agreement, which is exactly the region a
basis spends most of its time in: a B-spline is zero on most of its
interval. Values whose scale is at the level of rounding error are
excluded instead of being compared, since neither side carries
information there.
