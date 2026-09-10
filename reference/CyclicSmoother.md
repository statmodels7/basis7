# The Cyclic Smoother Class

The class
[`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md)
returns: a periodic B-spline basis with a roughness penalty and the
Demmler-Reinsch reparametrization. It carries `degree` beside the
properties every
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md)
has.

## Usage

``` r
CyclicSmoother(
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
  degree = integer(0)
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

  The degree of the underlying B-spline, `3` for a cubic.

## Value

An S7 object of class `CyclicSmoother`, inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).
Construct one with
[`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md).

## Examples

``` r
sm <- cyclic_smooth(k = 10)
c(class = class(sm)[1], dimension = sm@dimension)
#>                    class                dimension 
#> "basis7::CyclicSmoother"                     "10" 
```
