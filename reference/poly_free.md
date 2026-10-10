# The Free Columns of a Polynomial Null Space

Builds the columns that `null_space = "keep"` restores for a family
whose null space is the polynomials of degree below `order`: the powers
\\x, \ldots, x^{m}\\ made orthogonal to the constant and to one another
over the observed covariate, then standardized.

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

The first column is written as `(x - mean(x)) / sd(x)` instead of as the
general regression of which it is a case. The two agree to rounding, and
the literal form keeps the arithmetic of the default `order = 2`
unchanged.
