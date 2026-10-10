# The Roughness Matrix of a Smoother

The Gram matrix of the derivative of order `sm@order`, integrated
against the measure of the smoother: the length measure on the interval
for `"lebesgue"`, the empirical measure of the covariate for
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

  The basis returned by
  [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md).

- x:

  The covariate, for the empirical measure.

- ...:

  Passed to methods.

## Value

A symmetric numeric matrix of `b@dimension` rows and columns, or a list
of them where the family's roughness is a sum of components with a
smoothing parameter each, as
[`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md)'s
is.

## Details

The measure changes the penalty: the Gram matrices under two measures
are different matrices, and a fit under one measure differs from a fit
under another.

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

# another measure gives a different penalty
sme <- bspline_smooth(k = 8, lower = -2, upper = 2,
                      measure = "empirical")
round(cor(as.vector(g),
          as.vector(smoother_gram(sme, smoother_basis(sme, x), x))), 3)
#> [1] 0.995
```
