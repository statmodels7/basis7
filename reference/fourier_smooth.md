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

  What happens to the directions that the penalty does not see, other
  than the constant, which is always removed. `NULL` takes `"shrink"`,
  or `"keep"` where a `penalty` factory is given, the two being
  incompatible. See the section above.

- reparam:

  The coordinates in which the coefficients are expressed.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- omega:

  The period **of the basis**, `NULL` for the width of the interval. The
  period of the operator is carried by
  [`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md).

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
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and not a derivative, and the period that it uses is the width of the
interval. The two differ in what a strongly penalized fit contracts
**to**. On a periodic basis a derivative penalty has the constant as its
only null direction, so the fit contracts to a constant; the harmonic
operator leaves the level and the fundamental cycle unpenalized and
penalizes every higher harmonic, so the fit contracts to a sinusoid.

The default null space is `"shrink"` for this family and `"keep"` for
the others. With `"keep"` the sine and cosine of the fundamental are
free columns, so they are spent whether or not the covariate carries a
cycle; with `"shrink"` they are penalized lightly and can leave the
model. The first suits a cycle known to be present and estimated, the
second a cycle whose presence is a hypothesis.

`order = 2` gives the integrated squared second derivative.

## Arguments the family does not take

`degree` is a B-spline's, and a Fourier basis has none.

`constrain` in the form "the polynomials up to degree c" has no reading
here: a periodic basis contains no linear function, so a derivative
penalty's null space is the **constant at every order**, not a space
growing with the order: the null space of the Gram matrix of a Fourier
basis is one-dimensional at every order, spanned by the constant.

## Accepted operators

The constant is always removed, including when the operator penalizes
it, as
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
does, so a model carrying an intercept spans the level and the smooth
carries the shape. Every **other** function of the operator's null space
is restored as a free column, and it must be periodic on the period of
the basis, which means a pure sine or cosine at a whole multiple of the
fundamental frequency. An operator whose null space holds \\t\\, or
\\e^{at}\\, or a frequency that is not a multiple, is rejected, because
such a column would make the fit non-periodic.

So
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
are accepted at any number of harmonics that leaves at least one
function to penalize, which `fourier_smooth()` checks at construction;
`deriv_operator(m)` is accepted and restores nothing; and
`deriv_operator(2) * oscillator_operator(p)`, with a period `p`, is
rejected because its null space holds \\t\\, which suits a B-spline and
not a periodic basis.

## The interval is the period

The interval of a periodic basis is fixed by the modeler and not read
from the data, because the period is a property of the covariate and the
range of one sample does not determine it. `lower` and `upper` give it,
as in `fourier_smooth(k = 9, lower = 0, upper = 365)` for a day of the
year. Without them the interval is the range of the data padded by a
thousandth of its width, and the basis is periodic on that interval.

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
