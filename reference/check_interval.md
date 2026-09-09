# Check a Smoother's Interval Arguments

Validates the `lower` and `upper` arguments of a smoother constructor.
Each may be `NULL`, meaning the endpoint is read from the data at build.

## Usage

``` r
check_interval(lower, upper)
```

## Arguments

- lower, upper:

  The endpoints.

## Value

`NULL`, invisibly. Called for the error it signals.
