# The Periodicity Constraint of a Spline Basis

Returns the matrix whose rows say that a function of `basis` and its
first `degree - 1` derivatives take the same value at the two ends of
the interval: row \\j\\ is \\B^{(j)}(a) - B^{(j)}(b)\\ for \\j = 0,
\ldots, d - 1\\. Its null space is the periodic splines, which is what
[`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md)
builds on.

## Usage

``` r
periodic_constraint(basis, degree)
```

## Arguments

- basis:

  The parent basis, of dimension `k + degree`.

- degree:

  The degree of the spline, one row per derivative order below it.

## Value

A numeric matrix with `degree` rows and `basis@dimension` columns.

## Details

The rows are read from
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
at the two endpoints, order zero included, so the constraint is
evaluated by the same arithmetic that evaluates the basis rather than
from a formula about knots.
