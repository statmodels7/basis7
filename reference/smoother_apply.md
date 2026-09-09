# A Smoother's Block at New Values

Evaluates the construction recorded by
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
at new covariate values. The basis, the constraint and the
reparametrization are the ones computed on the original data: nothing is
recomputed from `newx`.

## Usage

``` r
smoother_apply(sm, blueprint, newx, ...)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md),
  the one passed to
  [`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md).

- blueprint:

  The `blueprint` element of that call's result.

- newx:

  The new covariate values, a numeric vector.

- ...:

  Passed to methods.

## Value

A numeric matrix of `length(newx)` rows and as many columns as the block
the blueprint came from, with no dimnames.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md),
which records the blueprint.

## Examples

``` r
set.seed(1)
x <- sort(runif(150, 0, 1))
sm <- bspline_smooth(k = 10)
out <- smoother_build(sm, x)

# Reapplied at a subset, the rows are the rows of the original block.
i <- c(3, 40, 91)
max(abs(smoother_apply(sm, out$blueprint, x[i]) - out$X[i, ]))
#> [1] 0
```
