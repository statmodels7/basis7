# Validate the Arguments Every Basis Constructor Takes

Checks the interval and the number of functions, and returns the
dimension as an integer. Called by every constructor before anything
family-specific.

## Usage

``` r
check_basis_args(lower, upper, dimension)
```

## Arguments

- lower, upper:

  The endpoints of the interval.

- dimension:

  The number of basis functions.

## Value

`dimension`, as an integer.
