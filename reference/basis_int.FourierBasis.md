# Integral of a Fourier Basis

Returns \\\int\_{\ell}^{x} \varphi_j(t)\\\mathrm{d}t\\ for every column
in closed form, with no quadrature. The constant integrates to \\x -
\ell\\ and the sinusoids come from the phase-shift identity of
[FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
taken at order \\-1\\.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
exactly zero in the row at `basis@lower`.

## Details

The identity at \\k = -1\\ gives *an* antiderivative, which is not the
one
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
promises. The two differ by a constant that is not the same in every
column: at the lower endpoint the sine columns of the raw antiderivative
are \\-\omega/(2\pi j)\\ while the cosine columns are already zero.
Subtracting the row at the lower endpoint corrects every column at once
and makes the anchoring exact.

Over a full period every sinusoid integrates to zero, so the row at
`basis@upper` is \\(\omega, 0, 0, \ldots)\\, which is the area under a
fitted curve being the constant's coefficient times the period.

## See also

[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention;
[`fourier_trig()`](https://statmodels7.github.io/basis7/reference/fourier_trig.md),
which supplies the raw antiderivative.
