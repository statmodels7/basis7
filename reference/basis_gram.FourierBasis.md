# Gram Matrix of a Fourier Basis

The inner products of the basis functions, diagonal in closed form when
the interval is a whole period.

## Arguments

- basis:

  A
  [`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- order:

  The derivative order.

- ...:

  Passed to the numerical method when the period is not the interval
  width.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

Over a whole period the trigonometric functions are mutually orthogonal,
and remain so after differentiation, since a derivative only shifts the
phase and rescales. The order-\\d\\ Gram matrix is therefore diagonal,
with entries \\\omega\\ for the constant at order zero, zero for it at
every higher order, and \\(\omega/2)(2\pi j/\omega)^{2d}\\ for both
members of pair \\j\\.

When `omega` is not the width of the interval the orthogonality fails,
and the numerical method of the
[`basis`](https://statmodels7.github.io/basis7/reference/basis.md) class
is used instead.
