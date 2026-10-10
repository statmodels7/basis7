# Default Column Names

Numbers the basis functions after the family, taking the first two
characters of `@basis_name` and appending the numbers from 1 to
`@dimension`: `bs1 ... bs6` for a B-spline. Every class inherits this
method unless it registers its own, as the Fourier, Legendre, tensor and
transformed bases do.

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`.

## Details

Two characters tell the families apart in a coefficient table without
making the names long. The names are not distinct across bases: a model
that combines two B-spline blocks has `bs1` twice, unless the code that
assembles the design distinguishes them.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic.
