# The Legendre Smoother Class

The class
[`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md)
returns: an orthogonal polynomial basis with a roughness penalty and the
Demmler-Reinsch reparametrization. Its null space is the polynomials of
degree below `order`, as a B-spline's is, so it inherits the default
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md).

## Usage

``` r
LegendreSmoother(
  smoother_name = character(0),
  dimension = integer(0),
  order = NULL,
  measure = NULL,
  constrain = NULL,
  null_space = character(0),
  reparam = character(0),
  penalty = NULL,
  lower = NULL,
  upper = NULL,
  smoother_params = list()
)
```

## Arguments

- smoother_name:

  A single string naming the family, printed by
  [`print.smoother()`](https://statmodels7.github.io/basis7/reference/print.smoother.md).
  Not read by any computation.

- dimension:

  The number of basis functions before any constraint or
  reparametrization, a single positive integer. Must be of storage mode
  integer.

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

- measure:

  The measure the roughness is integrated against. `"lebesgue"` is the
  length measure on the interval.

- constrain:

  The directions the smooth is made orthogonal to, or `NULL` for the
  null space of the basis and the penalty together. The form it takes
  belongs to the family, a periodic basis having no reading for
  "polynomials up to degree c".

- null_space:

  What becomes of the directions the penalty does not see: `"keep"`
  leaves them as free columns, `"drop"` removes them, `"shrink"`
  penalizes them under the same smoothing parameter.

- reparam:

  The coordinates the coefficients live in, a single string. `"dr"` is
  the Demmler-Reinsch rotation.

- penalty:

  `NULL` for the quadratic roughness matrix, or a factory building a
  penalties7 penalty from a coefficient count.

- lower, upper:

  The endpoints of the interval, each a single finite number or `NULL`
  to read it from the data at build. Given, both must be given, with
  `lower < upper`.

- smoother_params:

  A named list of whatever else the family needs.

## Value

An S7 object of class `LegendreSmoother`, inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).
Construct one with
[`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md).

## See also

[`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md),
which is the way to build one.

## Examples

``` r
sm <- legendre_smooth(k = 8)
c(class = class(sm)[1], dimension = sm@dimension)
#>                      class                  dimension 
#> "basis7::LegendreSmoother"                        "8" 
```
