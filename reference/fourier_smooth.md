# A Fourier Smoother

The periodic smoother: `k` Fourier functions over an interval, penalized
by the integrated squared derivative of order `order`, rotated to the
Demmler-Reinsch coordinates. Every column of the block is periodic, so a
fit built on it takes the same value at the two ends of the interval.

## Usage

``` r
fourier_smooth(
  k = 9,
  order = 2,
  measure = "lebesgue",
  null_space = "keep",
  reparam = "dr",
  penalty = NULL,
  omega = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- k:

  The number of basis functions, an **odd** whole number of at least 3:
  a Fourier basis holds a constant plus complete sine-cosine pairs.

- order:

  The order of derivative the penalty integrates.

- measure:

  The measure the roughness is integrated against.

- null_space:

  What becomes of the null space. It is the constant alone here, and the
  constant is removed whatever this says, so the argument is accepted
  for symmetry with the other families and changes nothing.

- reparam:

  The coordinates the coefficients live in.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory.

- omega:

  The period, `NULL` for the width of the interval.

- lower, upper:

  The interval, which for a periodic basis is the period. See the
  section above: give both.

## Value

An S7 object of class
[FourierSmoother](https://statmodels7.github.io/basis7/reference/FourierSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## What it does not have, and why

`degree` is a B-spline's, and a Fourier basis has none.

`constrain` in the form "the polynomials up to degree c" has no reading
here: a periodic basis contains no linear function, so its penalty's
null space is the **constant at every order**, not a space growing with
the order. Measured on the Gram matrix of a nine-function Fourier basis
at orders 1, 2 and 3, the null function has a standard deviation of
exactly zero, which is to say it is constant, and the null space is
one-dimensional in all three.

The constant is removed, so a model carrying an intercept spans the
level and the smooth carries the shape. Nothing is restored as a free
column, which is why `null_space` has no effect here and `k` functions
give `k - 1` columns.

## The interval is the period

A periodic basis needs its interval fixed by the modeller, not read from
the data: the period is a fact about the covariate, and the range of one
sample is not it. Give `lower` and `upper`, as in
`fourier_smooth(k = 9, lower = 0, upper = 365)` for a day of the year.
Without them the interval is the padded range of the data, which makes
the basis periodic on an interval nothing else knows about.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
for what it produces at data,
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the non-periodic family,
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
for the basis alone.

## Examples

``` r
# A periodic covariate, and a truth that is periodic on [0, 1].
set.seed(3)
x <- sort(runif(300))
f <- function(u) sin(2 * pi * u) + 0.4 * cos(4 * pi * u)
y <- f(x) + rnorm(300, sd = 0.2)

sm <- fourier_smooth(k = 9, lower = 0, upper = 1)
out <- smoother_build(sm, x)
dim(out$X)
#> [1] 300   8

# THE FIT IS PERIODIC, which is the property the basis was chosen for.
b <- solve(crossprod(out$X) + 0.01 * out$S, crossprod(out$X, y))
ends <- smoother_apply(sm, out$blueprint, c(0, 1)) %*% b
format(diff(as.vector(ends)), digits = 3)
#> [1] "-3.33e-16"
```
