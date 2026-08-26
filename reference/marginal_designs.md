# Evaluate Every Marginal at Its Own Column

Returns a list of the margins' design matrices, each evaluated at its
own column of the points, and each differentiated or integrated as
asked. The first step of
[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md).

## Usage

``` r
marginal_designs(basis, x, order = NULL, integral = FALSE)
```

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- order:

  An integer vector with one entry per variable, or `NULL` for the
  functions themselves.

- integral:

  `TRUE` to ask each margin for its anchored integral.

## Value

A list of `basis_nvar(basis)` numeric matrices, the \\j\\th with
`nrow(x)` rows and `basis@marginals[[j]]@dimension` columns.

## See also

[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md),
its only caller.
