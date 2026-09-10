# A Cyclic Smoother

The local periodic smoother: `k` B-spline functions over one period,
constrained so that the fit and its first `degree - 1` derivatives take
the same value at the two ends of the interval, penalized by the
integrated squared derivative of order `order` and rotated to the
Demmler-Reinsch coordinates. It is to
[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
what
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
is to
[`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md):
a local basis where the other is global, so a feature at one point of
the cycle leaves the rest of it alone.

## Usage

``` r
cyclic_smooth(
  k = 10,
  degree = 3,
  order = 2,
  measure = "lebesgue",
  null_space = "keep",
  reparam = "dr",
  penalty = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- k:

  The number of periodic functions, a whole number of at least 3. The
  block carries `k - 1` columns, the constant being removed.

- degree:

  The degree of the underlying B-spline, `3` for a cubic. The fit
  matches at the two ends in its value and its first `degree - 1`
  derivatives.

- order:

  What the penalty measures: a
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md)
  from
  [`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
  [`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
  [`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
  or
  [`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md),
  or a whole number `m` as the shorthand for `deriv_operator(m)`. It
  says what a strongly penalized fit contracts toward, which for `m` is
  a constant at 1, a straight line at 2 and a parabola at 3, and for any
  operator is
  [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md).
  Its order is at most `degree`.

- measure:

  The measure the roughness is integrated against.

- null_space:

  What becomes of the directions the penalty does not see. A periodic
  basis restores none, so the settings differ only in what they refuse.

- reparam:

  The coordinates the coefficients live in.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- lower, upper:

  The ends of one period, `NULL` to read them from the data. See the
  section above.

## Value

An S7 object of class
[CyclicSmoother](https://statmodels7.github.io/basis7/reference/CyclicSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## The construction

A spline of degree \\d\\ on \\\[a, b\]\\ is periodic exactly when
\$\$f^{(j)}(a) = f^{(j)}(b), \qquad j = 0, 1, \ldots, d - 1,\$\$ which
is \\d\\ linear conditions on its coefficients. The periodic splines are
therefore a subspace of an ordinary spline space, and the whole
construction is
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
applied to a B-spline of dimension `k + degree`, whose null space has
dimension `k`.

Building it this way rather than by folding a widened knot sequence is
what makes every quantity exact. The parent basis lives on \\\[a, b\]\\,
the period itself, so its Gram matrix integrates over one period and the
penalty needs no numerical fallback: the roughness matrix is \\T^\top G
T\\ with \\G\\ the parent's own exact Gram. A folded construction puts
its parent on a widened interval, and its Gram then integrates over more
than one period.

## What it buys over a Fourier basis, and where it buys nothing

Both families are periodic, so the choice between them is locality and
nothing else, and what it is worth depends entirely on the truth. The
measurement sweeps the concentration of a seasonal peak, \\\exp(\kappa
\cos t)\\, at 400 observations with both bases ten columns wide and both
smoothing parameters chosen by \\\mathrm{REML}\\, three samples per
setting, reporting the root mean square error against the truth:

|            |        |         |          |
|------------|--------|---------|----------|
| \\\kappa\\ | cyclic | fourier | ratio    |
| 1          | 0.0600 | 0.0612  | 1.02     |
| 3          | 0.0652 | 0.0655  | 1.00     |
| 8          | 0.0717 | 0.1783  | **2.49** |
| 20         | 0.3962 | 0.6156  | 1.55     |
| 50         | 0.9282 | 1.1168  | 1.20     |

At a low concentration the two are the same fit: \\\exp(\kappa \cos t)\\
is a von Mises density, whose Fourier coefficients are modified Bessel
functions and decay geometrically, so it is the function a global basis
is best at. The gain appears where the peak becomes narrow against what
ten columns can carry, and falls back again at \\\kappa = 50\\, where
neither basis of that width represents the spike at all and the
comparison stops being about locality. So the family is worth choosing
for a localized seasonal feature and is worth nothing for a smooth one,
which is the same reading
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
and
[`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md)
have against each other away from the circle.

Measured on a cubic with `k = 9` over \\\[0, 1\]\\: the basis functions
agree at the two ends to 5.6e-17 in value, 7.1e-15 in the first
derivative and 8.5e-14 in the second, where an ordinary B-spline of the
same `k` disagrees by 4.12. The roughness matrix agrees with a fine
trapezoid of the second derivatives over one period, the gap falling by
exactly 4.00 at each halving of the step, which is the reference's own
order of convergence rather than an error of the matrix. Integrating
over 99 per cent of the period instead moves that matrix by 167.6
against a size of 2820.4.

## What it does not have, and why

`constrain` in the form "the polynomials up to degree c" has no reading
here, for the reason it has none on
[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md):
a non-constant periodic function is never a polynomial, so the null
space of the roughness matrix is the constant at every order rather than
a space growing with the order. Measured on the pair at orders 1, 2 and
3, the null space is one-dimensional in all three and its function is
constant to 2.7e-15.

The constant is removed, so a model carrying an intercept spans the
level and the smooth carries the shape, and nothing is restored as a
free column: a linear column is not periodic, and restoring one is what
makes a fit lose the property the basis was chosen for. A block
therefore carries `k - 1` columns.

## The interval is the period

`lower` and `upper` are the ends of one cycle, and they are a property
of the problem rather than of the sample: day-of-year data run over
\\\[0, 365\]\\ whether or not an observation falls on the first day of
the year. Left `NULL` they are read from the data, which makes the
fitted period the observed range and is almost never what a periodic
model means.

## See also

[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
for the global periodic family,
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the local non-periodic one.

## Examples

``` r
set.seed(3)
doy <- sort(runif(200, 0, 365))
sm <- cyclic_smooth(k = 10, lower = 0, upper = 365)
out <- smoother_build(sm, doy)
dim(out$X)
#> [1] 200   9

# the block takes the same value at the two ends of the period
ends <- smoother_apply(sm, out$blueprint, c(0, 365))
max(abs(ends[1, ] - ends[2, ]))
#> [1] 1.110223e-16

# 'order' may not exceed the degree.
try(cyclic_smooth(k = 10, degree = 3, order = 4))
#> Error : 'order' (4) exceeds 'degree' (3): the derivative of that order of a
#>   spline of degree m is zero, so the penalty would be the zero matrix.
```
