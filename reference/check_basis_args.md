# Validate the Arguments Every Basis Constructor Takes

Checks the three arguments common to
[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
and
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md),
and returns the dimension coerced to integer, which is the storage mode
the class validator requires. Every constructor calls it first, before
anything family-specific, so the three report the same errors in the
same words.

## Usage

``` r
check_basis_args(lower, upper, dimension)
```

## Arguments

- lower, upper:

  The endpoints of the interval, each a single finite number with
  `lower < upper`.

- dimension:

  The number of basis functions, a single whole number of at least 1,
  integer or double.

## Value

`dimension` as an integer of length one.

## Details

Three conditions, each throwing with `call. = FALSE`:

- `lower` and `upper` must each be a single finite number, else
  `'lower' and 'upper' must be single finite numbers.`

- `lower` must be strictly below `upper`, else
  `'lower' must be strictly less than 'upper'.`

- `dimension` must be a single finite number, at least 1 and whole, else
  `'dimension' must be a single positive integer.`

A whole-valued double passes and is converted, so a caller writing `6`
rather than `6L` is served. The multivariate case is not reachable here:
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
builds its endpoints from margins already validated.

## See also

[`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md),
the same role for the points a basis is evaluated at.
