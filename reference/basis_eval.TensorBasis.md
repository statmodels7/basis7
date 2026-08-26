# Evaluate a Tensor Product Basis

Evaluates each margin at its own column of the points and takes the
row-wise Kronecker product, so row \\i\\ of the result is \\B_1(x\_{i1})
\otimes \cdots \otimes B_D(x\_{iD})\\. The matrix has \\\prod_j K_j\\
columns, so at several variables it is worth avoiding, and
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
avoids it.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable, or a vector taken by
  row. Each column is checked against its own margin's interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns, with
column names pasting the margins' with dots.

## Details

Cost is one marginal evaluation per variable plus the products, and the
result is `nrow(x)` by `basis@dimension`, which at four margins of eight
functions and 20000 points is 625 MB. Nothing is cached.

## See also

[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md)
and
[`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md),
which do the work;
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md),
for the value of a fit without this matrix.
