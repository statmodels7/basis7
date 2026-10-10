# Name the Columns of a Basis Matrix

Sets `colnames(m)` to
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
of the basis and returns the matrix. The methods that return an `n` by
`@dimension` matrix end with a call to it, so the evaluation, the
derivatives of every order and the integral of one basis carry the same
column names in the same order; the Gram methods set the same names on
both margins.

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

A matrix with the wrong number of columns signals R's error that the
length of the dimnames does not equal the extent of the array. The
function checks nothing else.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md),
which supplies the names.
