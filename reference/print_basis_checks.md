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

A verdict is `[PASSED]`, `[FAILED]`, or one of two things for an `NA`. A
check not run because the quantity is a finite difference prints
`[numerical]`; the partition-of-unity check on a family that is not one
prints `[not claimed]`. Saying which matters: a value that was not
verified is a gap, where a property never promised is not.

Only `deriv`, `integral` and `partition` can be `NA` today. The other
three rows carry a `[numerical]` label in the table that nothing
currently reaches, `shape`, `gram` and `missing` all being assigned on
every path.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its only caller.
