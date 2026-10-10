# Print the Outcome of check_basis

Prints the six-row table
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
shows when `verbose` is `TRUE`: a header naming the family and its
dimension, one line per check with a description and a verdict, and a
closing line listing the quantities that came from the numerical
fallback.

## Usage

``` r
print_basis_checks(basis, res, num)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md). Read
  for its `@basis_name` and `@dimension` only.

- res:

  The named logical vector of results, in the order `shape`, `deriv`,
  `integral`, `partition`, `gram`, `missing`.

- num:

  The named logical vector
  [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
  returns, printed as the closing line when any entry is `TRUE`.

## Value

`NULL`, invisibly. Called for the printed output.

## Details

A verdict is `[PASSED]`, `[FAILED]`, or one of two labels for an `NA`. A
check that was not run because the quantity is numerical prints
`[numerical]`; the partition-of-unity check on a family that is not a
partition of unity prints `[not claimed]`. The first marks a value that
was not verified, the second a property that does not apply.

Only `deriv`, `integral` and `partition` can be `NA`. The rows `shape`,
`gram` and `missing` are assigned on every path, so their `[numerical]`
label is never printed.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its only caller.
