# A Fourier Smoother

The periodic smoother: `k` Fourier functions over an interval, penalized
by a linear differential operator and rotated to the Demmler-Reinsch
coordinates. Every column of the block is periodic, so a fit built on it
takes the same value at the two ends of the interval.

## Usage

``` r
fourier_smooth(
  k = 9,
  order = harmonic_operator(),
  measure = "lebesgue",
  null_space = NULL,
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

  What the penalty measures: a
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  or a whole number `m` as the shorthand for `deriv_operator(m)`. The
  default is the harmonic acceleration operator at the period of the
  interval.

- measure:

  The measure the roughness is integrated against.

- null_space:

  What becomes of the directions the penalty does not see, other than
  the constant, which is always removed. `NULL` takes `"shrink"`, or
  `"keep"` where a `penalty` factory is given, the two being refused
  together. See the section above for what the choice costs.

- reparam:

  The coordinates the coefficients live in.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- omega:

  The period **of the basis**, `NULL` for the width of the interval. It
  is not the operator's, which
  [`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
  carries.

- lower, upper:

  The interval, which for a periodic basis is the period. See the
  section above: give both.

## Value

An S7 object of class
[FourierSmoother](https://statmodels7.github.io/basis7/reference/FourierSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## The penalty is the harmonic acceleration operator

The default `order` is
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
not a derivative, and the period it uses is the interval's. The reason
is what a strongly penalized fit contracts **to**: the derivative
penalty asks for a straight line, which is not periodic and is not what
a cyclic phenomenon simplifies to, while the harmonic operator leaves
the level and the fundamental cycle alone and removes everything above
them.

Measured at `k = 21` over 200 observations with the smoothing parameter
chosen by generalized cross-validation on eight samples: on a truth
dominated by its fundamental the harmonic penalty reaches a root mean
square error of 0.0596 against the derivative penalty's 0.0654 and wins
on seven samples of eight, at 9.28 effective degrees of freedom against
10.72. On a truth with no periodic signal at all the two are level,
0.0110 against 0.0121 at 1.16 degrees of freedom against 1.20, because
`null_space = "shrink"` lets the fundamental leave the model.

⚠️ That second row is why the default null space is `"shrink"` here and
`"keep"` elsewhere. With `"keep"` the fundamental is free, the fit
contracts to a sinusoid exactly as the operator promises, and the two
columns are spent whether or not the covariate carries a cycle: on that
same flat truth `"keep"` reads 0.0281 at 3.00 degrees of freedom and
wins on none of the eight. Use `"keep"` where the cycle is known to be
there and is the thing being estimated, and `"shrink"` where it is a
hypothesis.

`order = 2` restores the integrated squared second derivative, in one
argument, and gives the construction this family had before operators
existed.

## What it does not have, and why

`degree` is a B-spline's, and a Fourier basis has none.

`constrain` in the form "the polynomials up to degree c" has no reading
here: a periodic basis contains no linear function, so a derivative
penalty's null space is the **constant at every order**, not a space
growing with the order. Measured on the Gram matrix of a nine-function
Fourier basis at orders 1, 2 and 3, the null function has a standard
deviation of exactly zero, which is to say it is constant, and the null
space is one-dimensional in all three.

## Which operators this family accepts

The constant is always removed, so a model carrying an intercept spans
the level and the smooth carries the shape. Every **other** function of
the operator's null space is restored as a free column, and it must be
periodic on the basis's own period, which means a pure sine or cosine at
a whole multiple of the fundamental frequency. An operator whose null
space holds \\t\\, or \\e^{at}\\, or a frequency that is not a multiple,
is rejected: restoring such a column is what makes a fit on a periodic
basis lose the property the basis was chosen for.

So
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
are accepted at any number of harmonics the basis is wide enough to
hold, `deriv_operator(m)` is accepted and restores nothing, and
`deriv_operator(2) * oscillator_operator()` is rejected with the reason,
being right for a B-spline and wrong here.

## The interval is the period

A periodic basis needs its interval fixed by the modeller, not read from
the data: the period is a fact about the covariate, and the range of one
sample is not it. Give `lower` and `upper`, as in
`fourier_smooth(k = 9, lower = 0, upper = 365)` for a day of the year.
Without them the interval is the padded range of the data, which makes
the basis periodic on an interval nothing else knows about.

## See also

[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
for the default penalty,
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
sm
#> <basis7::FourierSmoother> 9 functions
#>   penalty: harmonic acceleration operator of order 3, lebesgue measure
#>   null space: shrink     coordinates: dr
#>   interval: [0, 1]
out <- smoother_build(sm, x)
dim(out$X)
#> [1] 300   8

# THE FIT IS PERIODIC, which is the property the basis was chosen for.
b <- solve(crossprod(out$X) + 0.01 * out$S, crossprod(out$X, y))
ends <- smoother_apply(sm, out$blueprint, c(0, 1)) %*% b
format(diff(as.vector(ends)), digits = 3)
#> [1] "-2.78e-16"

# THE FUNDAMENTAL IS WHAT A STRONG PENALTY LEAVES. With the null space
# kept it is free, and the heavily penalized fit is a pure sinusoid.
sk <- fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep")
ok <- smoother_build(sk, x)
ok$unpenalized
#> [1] 2
bk <- solve(crossprod(ok$X) + 1e8 * ok$S, crossprod(ok$X, y - mean(y)))
fv <- as.vector(ok$X %*% bk)
round(sqrt(mean(resid(lm(fv ~ sin(2 * pi * x) + cos(2 * pi * x)))^2)), 8)
#> [1] 3.7e-07

# The old construction, in one argument.
fourier_smooth(k = 9, order = 2, lower = 0, upper = 1)
#> <basis7::FourierSmoother> 9 functions
#>   penalty: derivative of order 2, lebesgue measure
#>   null space: shrink     coordinates: dr
#>   interval: [0, 1]

# Two harmonics left alone instead of one.
fourier_smooth(k = 13, lower = 0, upper = 365,
               order = harmonic_operator(harmonics = 2))
#> <basis7::FourierSmoother> 13 functions
#>   penalty: harmonic acceleration operator of order 5, lebesgue measure
#>   null space: shrink     coordinates: dr
#>   interval: [0, 365]

# An operator whose null space a periodic basis cannot carry.
try(smoother_build(
  fourier_smooth(k = 9, lower = 0, upper = 1,
                 order = deriv_operator(2) * oscillator_operator(1)), x))
#> Error : a periodic basis cannot carry every function of this operator's null space.
#>   It would restore t as free columns, and a column that is not periodic on
#>   the basis's period costs the fit the property the basis was chosen for.
#>   Use harmonic_operator() or oscillator_operator() at this period, a whole
#>   'order', or a non-periodic family such as bspline_smooth().
```
