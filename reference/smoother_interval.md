# The Interval of a Smoother

Returns the interval on which the basis of a smoother is built. An
endpoint stored on the object is used as it stands, and an endpoint
stored as `NULL` is the corresponding end of the range of `x`, moved
outward by a thousandth of the width of that range. The padding places
the observed values strictly inside the interval.

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

A list of two numbers, the lower and upper endpoints. A covariate that
takes a single value, with neither endpoint given, signals an error.
