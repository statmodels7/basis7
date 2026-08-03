# Row-Wise Kronecker Product of Two Matrices

For matrices with the same number of rows, the matrix whose \\i\\th row
is the Kronecker product of their \\i\\th rows, with the columns of the
second varying fastest.

## Usage

``` r
khatri_rao(a, b)
```

## Arguments

- a, b:

  Numeric matrices with the same number of rows.

## Value

A numeric matrix with `ncol(a) * ncol(b)` columns.
