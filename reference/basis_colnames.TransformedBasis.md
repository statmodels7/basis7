# Column Names of a Transformed Basis

The transform's prefix numbered from one. The parent's names cannot be
kept: each new function is a combination of all of them, and a transform
may produce fewer functions than it consumed.

## Arguments

- basis:

  A
  [`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- ...:

  Unused.

## Value

A character vector of length `basis@dimension`.
