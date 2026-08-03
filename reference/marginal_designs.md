# Evaluate Every Marginal at Its Own Column

Evaluate Every Marginal at Its Own Column

## Usage

``` r
marginal_designs(basis, x, order = NULL, integral = FALSE)
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

A list of numeric matrices, one per marginal.
