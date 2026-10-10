# Check a Smoother's Measure

Validates the `measure` argument of a smoother constructor: one of the
names that the package integrates against, or a function of the points.
For a function only the presence of an argument is checked here; its
weights are checked when the Gram matrix is computed.

## Usage

``` r
check_measure(measure)
```

## Arguments

- measure:

  The value given.

## Value

`measure`, unchanged.
