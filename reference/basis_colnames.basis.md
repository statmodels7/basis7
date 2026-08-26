# Default Column Names

Numbers the basis functions after the family, taking the first two
characters of `@basis_name` and appending `1` to `@dimension`:
`bs1 ... bs6` for a B-spline, `on1 ... on5` for an orthonormalized
basis. The method every class inherits unless it registers one of its
own, as
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
and
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
do.

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`.

## Details

Two characters is enough to tell the families apart at a glance in a
coefficient table without making the names long. Nothing depends on the
names being distinct across bases, and a model combining two B-spline
blocks will see `bs1` twice unless whatever assembles the design
disambiguates them.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic.
