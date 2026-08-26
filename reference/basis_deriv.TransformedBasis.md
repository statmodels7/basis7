# Derivatives of a Transformed Basis

Returns the parent's `order`-th derivative multiplied by the transform,
\\\tilde{B}^{(d)}(x) = B^{(d)}(x)\\T\\. Differentiation is linear and
\\T\\ does not depend on `x`, so the same matrix serves at every order
and no derivative of the transformation itself enters.

## Arguments

- basis:

  A
  [TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a single non-negative whole number, default `1`,
  passed to the parent unchanged.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

The accuracy is the parent's. Where the parent differentiates exactly so
does this; where the parent falls back to a stencil the result is that
stencil's answer arranged by \\T\\, which is why
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports the parent's flags for a transformed basis in place of its own
registered methods.

## See also

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
for whose accuracy this inherits.
