# Integral of a Legendre Basis

The exact integral from the lower endpoint.

## Arguments

- basis:

  A
  [`PolyBasis`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

From \\\int P_n = (P\_{n+1} - P\_{n-1})/(2n+1)\\ for \\n \ge 1\\ and
\\\int P_0 = t\\. The convention is satisfied without a correction: at
\\t = -1\\ the combination \\P\_{n+1} - P\_{n-1}\\ is \\(-1)^{n+1} -
(-1)^{n-1} = 0\\, and the constant term is anchored by construction.
