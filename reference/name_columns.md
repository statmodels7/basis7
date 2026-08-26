# Name the Columns of a Basis Matrix

Sets `colnames(m)` to
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
of the basis and returns the matrix. Every method returning an `n` by
`@dimension` matrix ends with a call to it, so the evaluation, the
derivatives of every order, the integral and the Gram matrix of one
basis all carry the same column names in the same order.

## Usage

``` r
name_columns(m, basis)
```

## Arguments

- m:

  A numeric matrix with exactly `basis@dimension` columns.

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

## Value

`m`, with its column names set and its contents untouched.

## Details

The names come from the basis, so a matrix with the wrong number of
columns produces R's own recycling error instead of a silently
mislabeled result. Nothing else is checked.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md),
which supplies the names.
