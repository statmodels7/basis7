# The Row-Wise Kronecker Product of the Marginal Designs

Evaluates each margin, differentiated or integrated as `order` and
`integral` specify, and reduces the results with
[`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md)
into the product's own design matrix. The one body behind
[`basis_eval.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.TensorBasis.md),
[`basis_deriv.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.TensorBasis.md)
and
[`basis_int.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.TensorBasis.md),
which differ only in the quantity computed from the margins.

## Usage

``` r
tensor_design(basis, x, order = NULL, integral = FALSE)
```

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable, already validated.

- order:

  An integer vector with one entry per variable, or `NULL` for the
  functions themselves.

- integral:

  `TRUE` to use the anchored integral of each margin instead of its
  evaluation. Not combined with `order`; each caller sets at most one.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns. Its
column names, inherited from the first margin, are replaced by the
callers through
[`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md).

## See also

[`marginal_designs()`](https://statmodels7.github.io/basis7/reference/marginal_designs.md)
and
[`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md),
its two steps.
