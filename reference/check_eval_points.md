# Validate Evaluation Points Against a Basis

Checks that `x` is numeric and lies inside the basis interval, reshapes
it for a basis of several variables, and returns it with near-endpoint
values clamped exactly onto the endpoints. Every method that evaluates a
basis calls it first, so the three families and the fallbacks agree on
what a legal evaluation point is. Missing values pass through and become
a missing row of the result.

## Usage

``` r
check_eval_points(basis, x)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  A numeric vector of evaluation points for a basis of one variable, or
  a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns for a basis of several. A one-column matrix is accepted for a
  basis of one variable and flattened; a vector is accepted for a basis
  of several and taken by row. `NA` is allowed anywhere and is neither
  range-checked nor clamped.

## Value

`x` with near-endpoint values clamped onto the endpoints: a numeric
vector for a basis of one variable, a numeric matrix of
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
columns otherwise.

## Details

A basis is defined on its interval and nowhere else, so a point outside
is rejected with an error naming how many points were outside and what
the interval is. No extrapolation rule is applied: which one a model
wants, if any, belongs to the layer that owns the meaning of the
covariate, and a silent answer here would take that decision away from
it.

The comparison carries a tolerance of `1e-8` times the width of the
interval, so a point that is an endpoint up to rounding is accepted and
then set to the endpoint exactly. On \\\[0, 1\]\\ that admits `1 + 5e-9`
and refuses `1 + 5e-8`; on \\\[0, 1000\]\\ it admits `1000 + 5e-6`. The
clamp matters because a B-spline evaluated a rounding step past its last
knot is zero in every column.

For a basis of several variables `x` is a matrix of one column per
variable, each column checked against its own endpoints. A plain vector
is accepted there and reshaped **by row**, so on two variables
`c(0.1, 0.2, 0.5, 0.6)` is the two points `(0.1, 0.2)` and `(0.5, 0.6)`.
A matrix of the wrong number of columns throws, naming the number
wanted.

## See also

[`clamp_to_range()`](https://statmodels7.github.io/basis7/reference/clamp_to_range.md),
which does the per-variable work;
[`check_basis_args()`](https://statmodels7.github.io/basis7/reference/check_basis_args.md),
the same role for a constructor's arguments.
