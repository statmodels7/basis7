# Column Names of a Transformed Basis

Numbers the columns under the two-letter prefix that the transformation
recorded: `on1 ... onk` after
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
`cn1 ...` after
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
`dr1 ...` after
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md).

## Arguments

- basis:

  A
  [TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`.

## Details

The parent's names cannot be carried over. Each new function is a
combination of all the parent's, so the parent's names do not fit it,
and a constraint returns fewer functions than it consumed. The prefix
records which transformation produced the column. The default
[`basis_colnames.basis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.basis.md)
takes the first two characters of `@basis_name`, which here begins with
the name of the transformation, and would give `or1` for both an
orthonormalization and its re-orthonormalization.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic.
