# Evaluate a Legendre Basis

Evaluates \\P_0, \ldots, P\_{K-1}\\ at the given points by the
three-term recurrence \$\$(n+1) P\_{n+1}(t) = (2n+1)\\ t\\ P_n(t) - n\\
P\_{n-1}(t),\$\$ after shifting the points from \\\[\ell, u\]\\ onto
\\\[-1, 1\]\\ by \\t = 2(x - \ell)/(u - \ell) - 1\\. The recurrence is
stable on that interval, every \\\|P_n\| \le 1\\ there, so no rescaling
is needed at any dimension.

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
with column names `P0`, `P1`, and so on.

## Details

One pass builds the whole table, column \\n+1\\ from columns \\n\\ and
\\n-1\\, at a cost of \\O(nK)\\ for `n` points. Missing points give
missing rows.

The values at the ends are exact: \\P_n(u) = 1\\ for every \\n\\, and
\\P_n(\ell) = (-1)^n\\.

## See also

[`legendre_table()`](https://statmodels7.github.io/basis7/reference/legendre_table.md),
which holds the recurrence;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
for the generic.
