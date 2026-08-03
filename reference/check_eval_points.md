# Validate Evaluation Points Against a Basis

Checks that `x` is numeric and lies inside the basis interval, and
returns it unchanged. Missing values are allowed and travel through to a
missing row.

## Usage

``` r
check_eval_points(basis, x)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points.

## Value

`x`, with near-endpoint values clamped onto the endpoints.

## Details

A basis is defined on its interval and nowhere else, so a point outside
it is refused rather than extrapolated. Which extrapolation rule a model
wants, if any, is a decision for the layer that knows what the covariate
means; a silent answer here would take that decision away from it.

The comparison uses a tolerance relative to the width of the interval,
so that a point which is an endpoint up to rounding is accepted and then
clamped onto the endpoint exactly.
