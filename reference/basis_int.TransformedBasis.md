# Integral of a Transformed Basis

Returns the parent's anchored integral multiplied by the transform. The
anchoring survives without a correction: every column of the parent's
integral is zero at the lower endpoint, and a linear combination of
zeros is zero.

## Arguments

- basis:

  A
  [TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
exactly zero in the row at `basis@lower`.

## Details

Integration is linear and \\T\\ does not depend on `x`, so the identity
is exactly the one
[`basis_deriv.TransformedBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.TransformedBasis.md)
uses in the other direction. The accuracy is again the parent's.

## See also

[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention.
