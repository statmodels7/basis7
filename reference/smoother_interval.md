# The Interval a Smoother Expands Over

Returns the interval a smoother's basis is built on: the endpoints
stored on the object when both are given, and otherwise the range of `x`
padded by a thousandth of its width. The padding keeps the observed
values strictly inside the interval, which is what a basis whose
validator requires an open interval needs at its endpoints.

## Usage

``` r
smoother_interval(sm, x)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- x:

  The covariate, a numeric vector.

## Value

A list of two numbers, the lower and upper endpoints.
