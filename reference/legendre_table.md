# The Legendre Polynomials by Recurrence

Builds a table of \\P_0, \ldots, P\_{k-1}\\ evaluated at `t` on \\\[-1,
1\]\\, one column per polynomial, from the three-term recurrence
\$\$(n+1) P\_{n+1}(t) = (2n+1)\\ t\\ P_n(t) - n\\ P\_{n-1}(t),\$\$
seeded with \\P_0 = 1\\ and \\P_1 = t\\. Used by
[`basis_eval.PolyBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.PolyBasis.md)
and, at one extra column, by
[`basis_int.PolyBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.PolyBasis.md).

## Usage

``` r
legendre_table(t, k)
```

## Arguments

- t:

  A numeric vector on \\\[-1, 1\]\\, the evaluation points already
  shifted.

- k:

  The number of polynomials, a whole number of at least 1.

## Value

A numeric matrix with `length(t)` rows and `k` columns, no dimnames.
Column 1 is the constant 1.

## Details

The recurrence is evaluated in the standard variable, never in the basis
interval's own, because \\\|P_n(t)\| \le 1\\ on \\\[-1, 1\]\\ and the
recurrence is stable there. Cost is \\O(nk)\\ for `n` points, one column
per step.

`t` is not range-checked here; the caller has already shifted validated
points onto \\\[-1, 1\]\\. Outside that interval the recurrence still
runs and the polynomials grow without bound.

## See also

[`legendre_deriv()`](https://statmodels7.github.io/basis7/reference/legendre_deriv.md),
which starts from this table.
