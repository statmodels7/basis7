# Derivatives of a Fourier Basis

Returns the `order`-th derivative of every column, exactly and at any
order, from the phase-shift identity of
[FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md):
differentiating \\\sin(jz)\\ shifts its phase by \\k\pi/2\\ and
multiplies it by \\(2\pi j/\omega)^{k}\\. The constant column is zero at
every order above 0.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a single non-negative whole number, default `1`.
  Any order is available.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
first column zero, with column names `const`, `sin1`, and so on.

## Details

Reaching order \\k\\ costs the same as reaching order 1: the shift and
the scale are both computed directly from \\k\\, with no recursion over
the orders below it. A Fourier basis therefore has no order at which its
derivatives stop being available, in contrast with a spline, whose
derivatives run out at its degree.

The scale grows as \\j^{k}\\, so a high frequency differentiated many
times is a large number: at `dimension = 21` and `order = 4` the largest
entry is \\(20\pi)^4\\, about 1.6e+06 on the unit interval. That is the
value, not a loss of accuracy.

## See also

[`fourier_trig()`](https://statmodels7.github.io/basis7/reference/fourier_trig.md),
which applies the identity;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic.
