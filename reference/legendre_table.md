# The Legendre Polynomials by Recurrence

A table of \\P_0, \ldots, P\_{k-1}\\ evaluated at `t` in \\\[-1, 1\]\\,
from the three-term recurrence.

## Usage

``` r
legendre_table(t, k)
```

## Arguments

- t:

  A numeric vector in \\\[-1, 1\]\\.

- k:

  The number of polynomials.

## Value

A numeric matrix with `length(t)` rows and `k` columns.
