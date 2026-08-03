# Compare Two Matrices Relative to Their Own Magnitude

Whether two matrices agree to a relative tolerance, with the denominator
taken from the values themselves rather than floored at one.

## Usage

``` r
rel_close(a, b, tol, slack = NULL)
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
interval. The denominator is therefore the values themselves.

It is floored, but at a millionth of the scale of the column it belongs
to rather than at one. A basis function's derivative crosses zero, and
at the crossing the pointwise value vanishes while the numerical
reference carries its usual rounding error; dividing that error by
nothing reports a failure of the reference as a failure of the basis.
Tying the floor to the curve's own magnitude keeps a proportional error
detectable wherever the curve is large, which is where a wrong formula
shows itself.

Columns whose whole scale is at the level of rounding error are skipped,
since neither side carries information there.
