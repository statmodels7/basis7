# A Smoother: the Four Decisions of a Penalized Smooth

A smoother is a recipe for a penalized smooth. It carries the four
choices that make up a smooth, which a single basis does not determine:

1.  which functions span the space, the **basis**;

2.  what counts as roughness, the **penalty**;

3.  which directions the penalty leaves unpenalized and what happens to
    them, the **null space**;

4.  in which coordinates the coefficients are expressed, the
    **reparametrization**.

It is a recipe and not a built object because the third and fourth
choices need the data: the Demmler-Reinsch rotation diagonalizes the
pencil of the empirical Gram matrix against the penalty, and the
interval of the default basis is read from the covariate.
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
resolves a smoother at a vector of covariate values and returns the
block and its penalty matrix;
[`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
reapplies the transform recorded there at new values.

## Usage

``` r
smoother(
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

  A single string naming the family. It is stored and validated;
  [`print.smoother()`](https://statmodels7.github.io/basis7/reference/print.smoother.md)
  prints the class name instead.

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
  determines the functions toward which a strongly penalized fit
  contracts: a constant for `m = 1`, a straight line for 2, a parabola
  for 3, and the functions of
  [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
  for any operator.

- measure:

  The measure against which the roughness is integrated. `"lebesgue"` is
  the length measure on the interval.

- constrain:

  The directions to which the smooth is made orthogonal, or `NULL` for
  the null space of the penalty. The form it takes belongs to the
  family, a periodic basis having no reading for "polynomials up to
  degree c".

- null_space:

  What happens to the directions that the penalty does not see: `"keep"`
  leaves them as free columns, `"drop"` removes them, `"shrink"`
  penalizes them under the same smoothing parameter with a weight of one
  tenth.

- reparam:

  The coordinates in which the coefficients are expressed, a single
  string. `"dr"` is the Demmler-Reinsch rotation.

- penalty:

  `NULL` for the quadratic roughness matrix, or a factory building a
  penalties7 penalty from a coefficient count.

- lower, upper:

  The endpoints of the interval, each a single finite number, or `NULL`
  to read that endpoint from the data at build. When both are given,
  `lower` must be less than `upper`.

- smoother_params:

  A named list of whatever else the family needs.

## Value

An S7 class object. `smoother` itself is abstract and cannot be
constructed; a family constructor such as
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
returns an object inheriting from it.

## One object for the basis and the penalty

The basis and the penalty cannot be chosen independently. The null space
is a property of the **pair**: the second-derivative Gram matrix of a
cubic B-spline has a two-dimensional null space, that of a Fourier basis
is one-dimensional at every order because the basis contains no linear
function, and the Gram matrix of a B-spline of degree \\d\\ at an order
above \\d\\ is identically zero, which penalizes nothing. With a basis
and a penalty passed separately, that compatibility would have to be
checked at every call site. A smoother is one object, so each
constructor validates its own arguments.

## A user-supplied penalty

`penalty` replaces the roughness matrix with a penalty built by a
factory of the coefficient count:
`bspline_smooth(penalty = penalties7::lasso_penalty)`. It is a factory
and not a built penalty because how many coefficients a smooth has is
settled by the data, the constraint and the null space moving with them.

A smoother **stores the function and never calls it**. Of the toolkit,
this package imports numericals7 only, so it cannot name penalties7 and
cannot test whether the function returns a penalty;
[`check_penalty()`](https://statmodels7.github.io/basis7/reference/check_penalty.md)
checks only that it is a function with at least one argument. The layer
that builds the term calls it, at the coefficient count that only the
data settle, and checks the result there.

The construction is unaffected.
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
returns the roughness matrix in `S` whether or not a factory is given,
because the reparametrization reads that matrix: it orders the
coordinates from the smoothest to the most wiggly and makes the penalty
on them the identity, and that ordering is what gives a penalty of
another shape its meaning. `unpenalized` counts the leading columns that
the roughness leaves free, which a caller building a separable penalty
needs, since such a penalty has no zero row with which to leave a column
unpenalized.

`smoother` is abstract: construct one through a family, of which
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
is the first.

## See also

[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the B-spline family,
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
for what a smoother produces at data,
[basis](https://statmodels7.github.io/basis7/reference/basis.md) for the
object a smoother builds on.

## Examples

``` r
sm <- bspline_smooth(k = 10)
S7::S7_inherits(sm, smoother)
#> [1] TRUE
sm@order
#> <linear differential operator>  order 2
#>   L x = D^2 x
#>   null space: 1, t

# Abstract: there is no direct constructor.
try(smoother())
#> Error in S7::new_object(S7::S7_object(), smoother_name = smoother_name,  : 
#>   Can't construct an object from abstract class <smoother>
```
