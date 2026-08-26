# Derivatives of the Legendre Polynomials

Computes the `order`-th derivative of \\P_0, \ldots, P\_{K-1}\\ with
respect to `x` on the basis interval. The one routine behind both
[`basis_eval.PolyBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.PolyBasis.md),
which calls it at order 0, and
[`basis_deriv.PolyBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.PolyBasis.md).

## Usage

``` r
legendre_deriv(basis, x, order)
```

## Arguments

- basis:

  A
  [PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- order:

  The derivative order, a single non-negative whole number. Order 0
  returns the polynomials themselves.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns, no
dimnames, `NA` in the rows where `x` is `NA`.

## Details

Differentiating \\P\_{n+1}' - P\_{n-1}' = (2n+1) P_n\\ repeatedly gives
\$\$P\_{n+1}^{(d)} = P\_{n-1}^{(d)} + (2n+1)\\ P_n^{(d-1)},\$\$ which
builds order \\d\\ from order \\d-1\\ in closed form, with no difference
taken anywhere. The seeds are \\P_0^{(d)} = 0\\ for \\d \ge 1\\ and
\\P_1' = 1\\, \\P_1^{(d)} = 0\\ above.

The loop runs `order` times over a `length(x)` by `K` table, so the cost
is \\O(d n K)\\, and it produces exact zeros once every polynomial has
been differentiated away. The shift onto \\\[-1, 1\]\\ contributes
\\(2/(u - \ell))^d\\, applied once after the loop.

## See also

[`legendre_table()`](https://statmodels7.github.io/basis7/reference/legendre_table.md),
which supplies order 0.
