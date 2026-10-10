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
  or a whole number `m` as the shorthand for `deriv_operator(m)`. The
  block is constrained against the constant only, so with a derivative
  penalty a strongly penalized fit contracts to a constant at every
  order. Its order is at most `degree`.

- measure:

  The measure against which the roughness is integrated.

- null_space:

  What happens to the directions that the penalty does not see. A
  periodic basis restores no free column, so the three settings give the
  same block, and `"shrink"` differs only in that it cannot be combined
  with a `penalty` factory.

- reparam:

  The coordinates in which the coefficients are expressed.

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

Building it this way instead of by folding a widened knot sequence makes
every quantity exact. The parent basis lives on \\\[a, b\]\\, the period
itself, so its Gram matrix integrates over one period and the penalty
needs no numerical fallback: the roughness matrix is \\T^\top G T\\ with
\\G\\ the parent's own exact Gram. A folded construction puts its parent
on a widened interval, and its Gram then integrates over more than one
period.

## Comparison with the Fourier basis

Both families are periodic and differ in locality: a Fourier basis is
global and a periodic B-spline basis is local. For a smooth periodic
function, such as \\\exp(\kappa \cos t)\\ at a low concentration
\\\kappa\\, whose Fourier coefficients are modified Bessel functions and
decay rapidly, the two give the same fit. A narrow seasonal feature, at
a higher concentration, is represented better by the local basis at the
same number of columns.

At the two ends of the period the basis functions agree in value and in
their first `degree - 1` derivatives, up to rounding.

## Arguments the family does not take

`constrain` in the form "the polynomials up to degree c" has no meaning
here, for the same reason as on
[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md):
a non-constant periodic function is never a polynomial, so the null
space of the roughness matrix of a derivative penalty is the constant at
every order, instead of a space growing with the order.

The constant is removed, so a model carrying an intercept spans the
level and the smooth carries the shape, and no free column is restored:
a linear column is not periodic. A block therefore carries `k - 1`
columns.

## The interval is the period

`lower` and `upper` are the ends of one cycle, and they are a property
of the problem rather than of the sample: day-of-year data run over
\\\[0, 365\]\\ whether or not an observation falls on the first day of
the year. Left `NULL` they are read from the data, which makes the
period the observed range.

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
#> [1] 0

# 'order' may not exceed the degree.
try(cyclic_smooth(k = 10, degree = 3, order = 4))
#> Error : 'order' (4) exceeds 'degree' (3): the derivative of that order of a
#>   spline of that degree is zero, so the penalty would be the zero matrix.
```
