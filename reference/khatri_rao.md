# Row-Wise Kronecker Product of Two Matrices

For two matrices with the same number of rows, returns the matrix whose
\\i\\th row is the Kronecker product of their \\i\\th rows. Reducing the
margins' design matrices with it gives the product's, one row of
\\\prod_j K_j\\ entries per observation.

## Usage

``` r
khatri_rao(a, b)
```

## Arguments

- a, b:

  Numeric matrices with the same number of rows. The row counts are not
  checked; a mismatch gives R's own recycling behavior.

## Value

A numeric matrix with `nrow(a)` rows and `ncol(a) * ncol(b)` columns, no
dimnames.

## Details

The columns come out with `b` varying fastest, matching
[`base::kronecker()`](https://rdrr.io/r/base/kronecker.html) and
[`basis_colnames.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.TensorBasis.md).
That is the whole column-order convention of a tensor basis, and the
reason
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
reverses an array's dimensions before flattening it.

The result has `ncol(a) * ncol(b)` columns, so reducing across several
margins grows geometrically; nothing here bounds that, and
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
is what avoids paying it.

## See also

[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md),
its only caller;
[`basis_colnames.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.TensorBasis.md),
which names its columns.
