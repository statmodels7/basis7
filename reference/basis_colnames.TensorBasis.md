# Column Names of a Tensor Product Basis

Pastes the margins' column names with dots, in the order the columns
come out: the last marginal varies fastest, following
[`base::kronecker()`](https://rdrr.io/r/base/kronecker.html). A column
named `bs2.cos1` is the product of the parent B-spline's second function
with the Fourier margin's first cosine, so a coefficient's name says
which marginal function it belongs to in each variable.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`.

## Details

The names are built with the same recycling
[`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md)
uses, the second factor repeated within each element of the first.
[`outer()`](https://rdrr.io/r/base/outer.html) would give the transpose
of this, which reads plausibly and labels every column but the first and
last wrongly.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic;
[`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md),
whose ordering this matches.
