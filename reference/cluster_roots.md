# Group Complex Roots by Multiplicity

Clusters the output of
[`base::polyroot()`](https://rdrr.io/r/base/polyroot.html) so that roots
within `tol` of one another, relative to their own size, are read as one
root of higher multiplicity.

## Usage

``` r
cluster_roots(rt, tol)
```

## Arguments

- rt:

  The complex roots.

- tol:

  The relative tolerance.

## Value

A list of `root` and `mult` pairs, ordered by the real part and then the
imaginary part.
