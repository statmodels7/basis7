# A Smoother: the Four Decisions of a Penalized Smooth

A smoother is a recipe for a penalized smooth. It carries the four
decisions a smooth is made of, which are independent of one another and
which a single basis does not determine:

1.  which functions span the space, the **basis**;

2.  what counts as roughness, the **penalty**;

3.  which directions the penalty leaves alone and what becomes of them,
    the **null space**;

4.  which coordinates the coefficients live in, the
    **reparametrization**.

It is a recipe rather than a built object because the third and fourth
decisions need the data: the Demmler-Reinsch rotation diagonalizes the
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
  order = integer(0),
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

## Value

An S7 class object. `smoother` itself is abstract and cannot be
constructed; a family constructor such as
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
returns an object inheriting from it.

## Why one object rather than a basis and a penalty

The two are not independently choosable. The null space is a property of
the **pair**: the second-derivative Gram matrix of a cubic B-spline has
a two-dimensional null space, that of a Fourier basis is one-dimensional
at every order because the basis contains no linear function, and the
Gram matrix of a B-spline of degree \\d\\ at an order above \\d\\ is
identically zero, which penalizes nothing. Passing a basis and a penalty
separately would leave the caller to satisfy that compatibility at every
call site. A smoother is one object, so each constructor validates its
own arguments where the caller wrote them.

## A penalty of your own

`penalty` replaces the roughness matrix with a penalty built by a
factory of the coefficient count:
`bspline_smooth(penalty = penalties7::lasso_penalty)`. It is a factory
and not a built penalty because how many coefficients a smooth has is
settled by the data, the constraint and the null space moving with them.

A smoother **stores the function and never calls it**. This package sits
at the bottom of the dependency graph and imports numericals7 alone, so
it cannot name penalties7 and cannot ask whether what the function
returns is a penalty;
[`check_penalty()`](https://statmodels7.github.io/basis7/reference/check_penalty.md)
asks only that it be a function of one argument. Whichever layer builds
the term calls it, at the count only the data settle, and checks the
result there.

The construction is unaffected.
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
returns the roughness matrix in `S` whether or not a factory is given,
because the reparametrization reads that matrix: it is what orders the
coordinates from the smoothest to the most wiggly and makes the penalty
on them the identity, and that ordering is the reason a penalty of
another shape is worth reaching for. `unpenalized` counts the columns
the roughness leaves free, which a caller building a separable penalty
needs, since such a penalty has no zero row with which to leave a column
alone.

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
#> [1] 2

# Abstract: there is no direct constructor.
try(smoother())
#> Error in S7::new_object(S7::S7_object(), smoother_name = smoother_name,  : 
#>   Can't construct an object from abstract class <smoother>
```
