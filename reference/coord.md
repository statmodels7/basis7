# One Coordinate of the Evaluation Points

A matched pair for reading and writing one variable of a set of
evaluation points. `coord(x, j)` returns the `j`th variable;
`replace_coord(x, j, z)` returns the points with that variable replaced
by `z` and the others left alone.
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
uses them to sweep one coordinate at a time when the basis takes
several.

## Usage

``` r
coord(x, j)

replace_coord(x, j, z)
```

## Arguments

- x:

  A numeric vector of evaluation points, or a matrix of one column per
  variable.

- j:

  The coordinate to read or write, a column index. Ignored when `x` is a
  vector.

- z:

  The replacement values, for `replace_coord()` only: a numeric vector
  as long as `x` has rows.

## Value

`coord()` returns a numeric vector: column `j` of `x`, or `x` itself
when it is not a matrix. `replace_coord()` returns points of the same
shape as `x`: the matrix with column `j` overwritten, or `z` itself when
`x` is not a matrix.

## Details

Both are the identity for a basis of one variable, where `x` is a plain
vector with no coordinate to pick: `coord()` returns `x` and
`replace_coord()` returns `z`, so the caller writes one loop over
`seq_len(basis_nvar(basis))` with no branch for the univariate case.

Neither validates anything. `j` outside the columns of `x` gives R's own
subscript error, and a `z` of the wrong length is recycled by R's usual
rules.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
their only caller, and
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
for the count they loop over.
