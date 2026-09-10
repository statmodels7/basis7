# The Roughness Matrix a Smoother Penalizes With

The Gram matrix of the derivative of order `sm@order`, integrated
against the measure the smoother carries: the length measure on the
interval for `"lebesgue"`, the empirical measure of the covariate for
`"empirical"`, and the density a function gives otherwise.

## Usage

``` r
smoother_gram(sm, b, x, ...)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- b:

  The basis
  [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md)
  returned.

- x:

  The covariate, for the empirical measure.

- ...:

  Passed to methods.

## Value

A symmetric numeric matrix of `b@dimension` rows and columns.

## Details

The measure is not decorative. Measured on a cubic B-spline of twelve
functions over \\\[-2, 2\]\\ at order 2, the correlation between the
Lebesgue Gram matrix and the one weighted by a Gaussian of standard
deviation 0.25 is 0.11: they are different penalties, and a fit under
one is not a fit under the other.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md),
which calls it, and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
the generic it reads.

## Examples

``` r
set.seed(1)
x <- sort(runif(100, -2, 2))
sm <- bspline_smooth(k = 8, lower = -2, upper = 2)
g <- smoother_gram(sm, smoother_basis(sm, x), x)
dim(g)
#> [1] 8 8

# the measure is a different penalty, not a detail
sme <- bspline_smooth(k = 8, lower = -2, upper = 2,
                      measure = "empirical")
round(cor(as.vector(g),
          as.vector(smoother_gram(sme, smoother_basis(sme, x), x))), 3)
#> [1] 0.995
```
