# A Legendre Smoother

The orthogonal polynomial smoother: `k` shifted Legendre polynomials
over an interval, penalized by the integrated squared derivative of
order `order`, rotated to the Demmler-Reinsch coordinates. It is a
global basis where a B-spline is local: every function is supported on
the whole interval, so a feature at one end moves the fit at the other.

## Usage

``` r
legendre_smooth(
  k = 8,
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

  The number of polynomials, a whole number of at least 2. The highest
  degree is `k - 1`.

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
  determines the functions toward which a strongly penalized fit
  contracts: a constant for `m = 1`, a straight line for 2, a parabola
  for 3, and the functions of
  [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
  for any operator. Its order is at most `k - 1`, the highest degree the
  basis carries.

- measure:

  The measure against which the roughness is integrated.

- constrain:

  The directions to which the smooth is made orthogonal, `NULL` for the
  null space of the penalty.

- null_space:

  What happens to the directions that the penalty does not see.

- reparam:

  The coordinates in which the coefficients are expressed.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- lower, upper:

  The interval, `NULL` to read it from the data.

## Value

An S7 object of class
[LegendreSmoother](https://statmodels7.github.io/basis7/reference/LegendreSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## Details

The null space of the roughness matrix is the polynomials of degree
below `order`, exactly as for a B-spline, so `order` and `constrain`
mean what they mean there and `null_space = "keep"` restores `order - 1`
free columns.

## See also

[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the local family,
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
for the basis alone.

## Examples

``` r
set.seed(3)
x <- sort(runif(200, -1, 1))
out <- smoother_build(legendre_smooth(k = 8), x)
dim(out$X)
#> [1] 200   7
out$names
#> [1] "lin" "z1"  "z2"  "z3"  "z4"  "z5"  "z6" 

# 'order' may not exceed what the basis can differentiate away.
try(legendre_smooth(k = 3, order = 4))
#> Error : 'order' (4) exceeds the highest degree the basis carries (2): the
#>   derivative of that order is zero, so the penalty would be the zero matrix.
```
