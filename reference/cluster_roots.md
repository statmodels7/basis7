# Group Complex Roots by Multiplicity

Clusters the output of
[`base::polyroot()`](https://rdrr.io/r/base/polyroot.html) so that roots
within `tol` of one another are read as one root of higher multiplicity.
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
passes the roots of the polynomial scaled to unit root size, so that
`tol` is relative to the size of the roots.

## Usage

``` r
cluster_roots(rt, tol)
```

## Arguments

- rt:

  The complex roots, scaled to unit size.

- tol:

  The tolerance.

## Value

A list of `root` and `mult` pairs, ordered by the real part and then the
imaginary part.
