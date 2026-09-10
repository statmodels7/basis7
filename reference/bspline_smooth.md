# A B-Spline Smoother

The B-spline smoother: `k` B-spline functions of degree `degree` over an
interval, penalized by the integrated squared derivative of order
`order`, rotated to the Demmler-Reinsch coordinates and carrying the
unpenalized direction as a free column. It is the construction
`modelterms7::s()` has always used, named as an object.

## Usage

``` r
bspline_smooth(
  k = 10,
  degree = 3,
  order = 2,
  measure = "lebesgue",
  constrain = NULL,
  null_space = "keep",
  reparam = "dr",
  penalty = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- k:

  The number of basis functions, a whole number of at least 3. Two
  directions are removed by the constraint, so a smaller `k` leaves
  nothing to smooth.

- degree:

  The degree of the B-spline pieces, a whole number of at least 1. `3`,
  the default, is the cubic spline.

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
  A spline of degree `d` has no derivative above `d`, so an operator of
  order above `degree` is rejected.

- measure:

  The measure the roughness is integrated against.

- constrain:

  The directions the smooth is made orthogonal to. `NULL`, the default,
  is the null space of the basis and the penalty together.

- null_space:

  What becomes of the directions the penalty does not see: `"keep"`
  leaves them as free columns, `"drop"` removes them, `"shrink"`
  penalizes them under the same smoothing parameter.

- reparam:

  The coordinates the coefficients live in. `"dr"`, the default, is the
  Demmler-Reinsch rotation.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- lower, upper:

  The interval. `NULL`, the default for each, reads it from the data at
  build; give both to fix it, which is what a prediction outside the
  observed range needs.

## Value

An S7 object of class
[BsplineSmoother](https://statmodels7.github.io/basis7/reference/BsplineSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).
It is a recipe: pass it to
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
with the covariate to obtain the block and its penalty matrix.

## The construction

At a vector of covariate values
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
performs, in order:

1.  a
    [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
    of `k` functions over the interval, which is `lower` and `upper`
    when both are given and otherwise the range of the data padded by a
    thousandth of its width;

2.  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
    which restricts the basis to the orthogonal complement of the
    constant and the linear function over the observed values, then
    diagonalizes the pencil of the empirical Gram matrix against the
    roughness matrix. The result is a basis whose columns are orthogonal
    over the data and ordered from the smoothest to the most
    oscillatory, and whose penalty is the identity;

3.  with `null_space = "keep"`, the standardized covariate prepended as
    a free column, so the penalty is `diag(0, 1, ..., 1)` and a strongly
    penalized fit contracts to a straight line rather than to a
    constant.

A basis of `k` functions therefore gives `k - 1` columns when the null
space is kept and `k - 2` when it is dropped: the constraint against the
constant and the linear function removes two directions, and the free
column adds one back.

## What `order` means

`order` is the order of derivative the penalty integrates, and it fixes
what a fit contracts toward as the smoothing parameter grows: a constant
at `order = 1`, a straight line at `order = 2`, a parabola at
`order = 3`. The null space of the penalty is the polynomials of degree
below `order`, so `order` may not exceed `degree`; above it the
roughness matrix is identically zero and penalizes nothing.

## References

Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
smoothing. *Numerische Mathematik*, 24, 375–382.

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science*, 11, 89–121.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
for what it produces at data,
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md)
for the four decisions it carries,
[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
for the basis alone.

## Examples

``` r
sm <- bspline_smooth(k = 10)
sm
#> <basis7::BsplineSmoother> 10 functions, degree 3
#>   penalty: derivative of order 2, lebesgue measure
#>   null space: keep     coordinates: dr
#>   interval: from the data

# k = 10 gives nine columns: the free linear column and eight deviations.
set.seed(1)
x <- sort(runif(200, -2, 3))
out <- smoother_build(sm, x)
dim(out$X)
#> [1] 200   9
out$names
#> [1] "lin" "z1"  "z2"  "z3"  "z4"  "z5"  "z6"  "z7"  "z8" 

# The penalty is diag(0, 1, ..., 1), and one column is unpenalized.
round(diag(out$S), 6)
#> [1] 0 1 1 1 1 1 1 1 1
out$unpenalized
#> [1] 1

# The free column is the standardized covariate, orthogonal to the rest
# over the observed values.
round(cor(out$X[, 1], x), 12)
#> [1] 1
max(abs(crossprod(out$X[, 1], out$X[, -1])))
#> [1] 2.198242e-14

# Dropping the null space removes it.
dim(smoother_build(bspline_smooth(k = 10, null_space = "drop"), x)$X)
#> [1] 200   8

# A penalty factory is STORED AND NEVER CALLED here: one that raises
# still builds, because it is the model layer that calls it.
sm2 <- bspline_smooth(k = 10, penalty = function(n_coef) stop("not here"))
out2 <- smoother_build(sm2, x)
identical(out2$S, out$S)
#> [1] TRUE

# It must be a function of the count, and it cannot be combined with a
# shrunk null space, which is a weight inside the matrix it replaces.
try(bspline_smooth(k = 10, penalty = 3))
#> Error : 'penalty' must be NULL, or a function of the number of coefficients giving a
#>   penalty. A penalties7 constructor passes bare -- penalty = penalties7::lasso_penalty
#>   -- and anything else is written out, as function(n_coef) my_penalty(n_coef).
try(bspline_smooth(k = 10, penalty = function(n) n, null_space = "shrink"))
#> Error : 'penalty' and null_space = "shrink" cannot both be given: the shrinkage is
#>   a weight inside the roughness matrix, and a penalty factory replaces that
#>   matrix. Use null_space = "keep" to leave the unpenalized directions free,
#>   or "drop" to remove them.

# 'k' must leave something after the constraint.
try(bspline_smooth(k = 2))
#> Error : 'k' (2) is too small for 'degree' (3): a B-spline basis of degree m
#>   needs at least m + 1 functions.
```
