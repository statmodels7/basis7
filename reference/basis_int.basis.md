# Numerical Integral of a Basis

The integration method every basis inherits from the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
composite Gauss-Legendre from the lower endpoint, accumulated over the
sorted evaluation points, so the whole set costs one pass where a
quadrature each would cost as many.

## Arguments

- basis:

  A basis object of one variable, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md). More
  than one variable throws.

- x:

  A numeric vector of evaluation points inside the basis interval, the
  upper limits of the integrals. Need not be sorted or unique. `NA`
  gives an `NA` row.

- nodes:

  The number of Gauss-Legendre nodes per segment, default `12`, so a
  polynomial integrand of degree up to 23 is integrated exactly.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md),
exactly zero in the row at `basis@lower`.

## Details

Numerical Integral of a Basis

## How it is accumulated

The points are sorted and made unique, the rule is placed on each
segment between consecutive ones, and the segment integrals are
cumulated. A point equal to `basis@lower` gives an empty first segment,
so its row is exactly zero and the anchoring convention of
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
holds by construction with no cancellation behind it. Duplicated points
are computed once and matched back.

## Accuracy

The `nodes`-point rule integrates a polynomial of degree up to
`2 * nodes - 1` exactly on each segment, so on a polynomial family the
result is exact to rounding: measured against the closed-form Legendre
integrals at 21 points, 3.3e-16. On a family that is not polynomial the
error is the rule's own on each segment, which shrinks with the spacing
of the evaluation points; a single distant point is integrated by one
rule over the whole span.

## One variable only

A basis of several variables throws. The integral there is over a box,
one iterated integral per variable, and a family that wants it supplies
its own, as
[`basis_int.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.TensorBasis.md)
does from its margins.

## See also

[`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md),
which places the rule;
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention.
