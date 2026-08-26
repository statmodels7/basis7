# Integral of a Legendre Basis

Returns \\\int\_{\ell}^{x} P_n(t)\\\mathrm{d}t\\ for every polynomial,
in closed form and with no quadrature. The zero at the lower endpoint
that
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
promises falls out of the identity used and needs no correction term.

## Arguments

- basis:

  A
  [PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
exactly zero in the row at `basis@lower`.

## Details

The identity is \\\int P_n = (P\_{n+1} - P\_{n-1})/(2n+1)\\ for \\n \ge
1\\, with \\\int P_0 = t\\, so one table of \\K + 1\\ polynomials
supplies all \\K\\ integrals: the integral of the last one reaches one
degree beyond the basis.

The anchoring is automatic. At \\t = -1\\, which is the lower endpoint,
\\P\_{n+1} - P\_{n-1}\\ is \\(-1)^{n+1} - (-1)^{n-1} = 0\\, and the \\n
= 0\\ column is written as \\(t + 1)/2 \cdot (u - \ell)\\, which also
vanishes there. Every column is therefore exactly zero at `basis@lower`,
with no cancellation of large numbers behind it.

## See also

[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention;
[`legendre_table()`](https://statmodels7.github.io/basis7/reference/legendre_table.md),
which supplies the polynomials.
