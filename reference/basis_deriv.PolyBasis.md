# Derivatives of a Legendre Basis

Returns the `order`-th derivative of every Legendre polynomial, exactly
and at any order, from the recurrence \\P\_{n+1}^{(d)} =
P\_{n-1}^{(d)} + (2n+1) P_n^{(d-1)}\\. No difference is taken anywhere,
so the fourth derivative is as accurate as the first, and an order above
`basis@basis_params$degree` returns exact zeros.

## Arguments

- basis:

  A
  [PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a single non-negative whole number, default `1`.
  Above `dimension - 1` the result is exactly zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names `P0`, `P1`, and so on.

## Details

Each order is built from the one below it, so reaching order \\d\\ costs
\\d\\ passes over a table of \\K\\ columns. The shift onto \\\[-1, 1\]\\
contributes a factor \\(2/(u - \ell))^d\\ by the chain rule, applied
once at the end.

## See also

[`legendre_deriv()`](https://statmodels7.github.io/basis7/reference/legendre_deriv.md),
which holds the recurrence;
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the generic.
