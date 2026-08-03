# The Row-Wise Kronecker Product of the Marginal Designs

Evaluates each marginal, differentiated or integrated as asked, and
multiplies the results row by row.

## Usage

``` r
tensor_design(basis, x, order = NULL, integral = FALSE)
```

## Arguments

- basis:

  A
  [`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- order:

  An integer vector of derivative orders, or `NULL`.

- integral:

  Whether to integrate instead.

## Value

A numeric matrix with `nrow(x)` rows.
