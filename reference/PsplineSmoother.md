# The P-spline Smoother Class

The class
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
returns: a B-spline basis with a difference penalty on its coefficients.
It carries `degree` and `diff` beside the properties every
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md)
has.

## Usage

``` r
PsplineSmoother(
  smoother_name = character(0),
  dimension = integer(0),
  order = integer(0),
  measure = NULL,
  constrain = NULL,
  null_space = character(0),
  reparam = character(0),
  penalty = NULL,
  lower = NULL,
  upper = NULL,
  smoother_params = list(),
  degree = integer(0),
  diff = integer(0)
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

  The order of derivative the penalty integrates, a single positive
  integer. It says what a strongly penalized fit contracts toward: a
  constant at 1, a straight line at 2, a parabola at 3.

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

- degree:

  The degree of the B-spline.

- diff:

  The order of difference the penalty takes.

## Value

An S7 object of class `PsplineSmoother`, inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).
Construct one with
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md).

## Examples

``` r
sm <- pspline_smooth(k = 20)
c(class = class(sm)[1], dimension = sm@dimension)
#>                     class                 dimension 
#> "basis7::PsplineSmoother"                      "20" 
```
