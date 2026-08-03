# Contract a Tensor Product Basis Against Factor Matrices

The canonical polyadic contraction: with \\\Gamma_j\\ of size \\K_j
\times F\\, the value is \\\sum_f \prod_j B_j(x_j)^\top \gamma\_{j,f}\\.

## Usage

``` r
contract_cp(basis, x, coef)
```

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- coef:

  A list of factor matrices, one per marginal, with a common number of
  columns.

## Value

A numeric vector with one value per row of `x`.

## Details

Each marginal is evaluated once and multiplied by its own factor matrix,
giving \\D\\ matrices of size \\n \times F\\; their elementwise product,
summed across columns, is the answer. Neither the design matrix nor the
coefficient array appears anywhere, so the cost is linear in the number
of variables where the array is exponential in it.
