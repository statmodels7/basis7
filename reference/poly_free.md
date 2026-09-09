# The Free Columns of a Polynomial Null Space

Builds the columns `null_space = "keep"` restores for a family whose
null space is the polynomials of degree below `order`: the powers \\x,
\ldots, x^{m}\\ made orthogonal to the constant and to one another over
the observed covariate, then standardized.

## Usage

``` r
poly_free(x, m)
```

## Arguments

- x:

  The covariate, a numeric vector.

- m:

  How many columns, `order - 1`. Zero gives a matrix of no columns.

## Value

A list of `free`, the matrix of `m` columns, and `params`, what
[`poly_free_apply()`](https://statmodels7.github.io/basis7/reference/poly_free_apply.md)
needs to rebuild it at new values.

## Details

The first column is written as `(x - mean(x)) / sd(x)` rather than as
the general regression it is a case of. The two agree to the last bit
for the quantities themselves, and the literal form is what the
construction has always computed at `order = 2`, which is every smooth
the toolkit has fitted; a change of arithmetic there would move fits
that are not being asked to move.
