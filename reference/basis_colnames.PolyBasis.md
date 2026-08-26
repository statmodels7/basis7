# Column Names of a Legendre Basis

Names the columns `P0`, `P1`, ..., `P(dimension - 1)`, so that a
column's name is the degree of the polynomial in it. This replaces the
default `le1`, `le2` of
[`basis_colnames.basis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.basis.md),
where the number would be the degree plus one and a reader would have to
remember the offset.

## Arguments

- basis:

  A
  [PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`, `"P0"` first.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic.
