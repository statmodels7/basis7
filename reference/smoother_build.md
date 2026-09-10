# Build a Smoother at Data

Resolves a
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md)
at a vector of covariate values and returns the design block, the
penalty matrix and what a caller needs to reapply the construction at
new values. It is the whole of a penalized smooth's arithmetic: what
comes back is the pair \\(X, S)\\ of \\(X'X + \lambda S)^{-1} X'y\\.

## Usage

``` r
smoother_build(sm, x, ...)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md),
  such as one from
  [`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md).

- x:

  The covariate, a numeric vector with no missing values.

- ...:

  Passed to methods.

## Value

A list of five elements:

- `X`:

  The design block, a numeric matrix of `length(x)` rows and one column
  per coordinate, with no dimnames.

- `S`:

  The penalty: a matrix square of `ncol(X)`, or, for a family whose
  roughness is a sum of components carrying a smoothing parameter each
  ([`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md)),
  a list of such matrices.

- `unpenalized`:

  The number of leading columns the penalty does not cover, as an
  integer.

- `blueprint`:

  What
  [`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
  needs to reapply the construction at new values.

- `names`:

  One name per column of `X`.

## Details

The block carries no column names, no `by` expansion and no label: those
belong to whichever model layer places the smooth in a formula. The
column names the family gives its own coordinates are returned
separately in `names`.

The transform is computed on `x` and recorded in `blueprint`, so a
prediction at new values reapplies it through
[`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
rather than rebuilding it. Rebuilding would give a basis over a
different interval and a rotation of a different Gram matrix, which is a
different function of the covariate.

## See also

[`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
for the block at new values,
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the smoother this builds.

## Examples

``` r
set.seed(1)
x <- sort(runif(150, 0, 1))
out <- smoother_build(bspline_smooth(k = 10), x)
vapply(out, function(e) class(e)[1], character(1))
#>           X           S unpenalized   blueprint       names 
#>    "matrix"    "matrix"   "integer"      "list" "character" 
dim(out$X)
#> [1] 150   9

# (X, S) is everything a penalized least-squares fit needs.
y <- sin(2 * pi * x) + rnorm(150, sd = 0.2)
b <- solve(crossprod(out$X) + 0.1 * out$S, crossprod(out$X, y))
round(sqrt(mean((out$X %*% b - sin(2 * pi * x))^2)), 3)
#> [1] 0.056

# The penalty leaves the first column free.
out$unpenalized
#> [1] 1
round(diag(out$S), 6)
#> [1] 0 1 1 1 1 1 1 1 1
```
