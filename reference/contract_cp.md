# Contract a Tensor Product Basis Against Factor Matrices

The canonical polyadic contraction. With \\\Gamma_j\\ of size \\K_j
\times F\\, the coefficient array is a sum of \\F\\ outer products and
the value is \\\sum_f \prod_j B_j(x_j)^\top \gamma\_{j,f}\\. Neither the
design matrix nor the coefficient array appears anywhere.

## Usage

``` r
contract_cp(basis, x, coef)
```

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- coef:

  A list of `basis_nvar(basis)` numeric matrices, the \\j\\th with
  `basis@marginals[[j]]@dimension` rows, all with the same number of
  columns \\F\\.

## Value

A numeric vector with one value per row of `x`.

## Details

Each margin is evaluated once and multiplied by its own factor matrix,
giving \\D\\ matrices of size \\n \times F\\; their elementwise product,
summed across columns, is the answer. The cost is \\O(nF\sum_j K_j)\\ in
time and memory, linear in the number of variables where the array is
exponential in it: measured at four margins of eight functions, 20000
points and \\F = 3\\, 0.03 s against the blocked array route's 0.83 s.

Against the array it stands for, the two routes agree to 1.7e-16.

## References

Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
*Proceedings of AISTATS*.

## See also

[`basis_contract.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_contract.TensorBasis.md),
its only caller.
