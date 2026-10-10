# Gram Matrix by Composite Quadrature

Computes the inner products of the order-`order` derivatives by
Gauss-Legendre on equal subintervals. Shared by
[`basis_gram.basis()`](https://statmodels7.github.io/basis7/reference/basis_gram.basis.md),
by
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
when the period is not the interval width, and by
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
as the independent reference it compares a family's own Gram matrix
against.

## Usage

``` r
numerical_gram(basis, order = 0L, panels = 50L, nodes = 12L)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- order:

  The derivative order, default `0`. Passed through
  [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md),
  so a single non-zero order on a basis of several variables signals an
  error.

- panels:

  The number of equal subintervals, default `50`, divided among the
  coordinates as above for a basis of several variables.

- nodes:

  The number of Gauss-Legendre nodes per subinterval, default `12`,
  exact for a polynomial integrand of degree up to 23 on each panel.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## One variable

The interval is cut into `panels` equal pieces with an `nodes`-point
rule on each, and the matrix is formed as a crossproduct of
\\\sqrt{w_i}\\B^{(d)}(t_i)\\, which keeps it positive semidefinite
whatever the integrand does. It is then symmetrized as `(G + t(G))/2`.

## Several variables

The rule is a product over the box: the nodes are the product grid of
the marginal rules and the weights their products. Each coordinate gets
`max(2, ceiling(panels^(1/d)))` panels, so the total node count stays
near `panels * nodes^d` and does not grow as `panels^d`. At the defaults
on two variables that is 8 panels of 12 nodes per coordinate, 9216
points.

## Accuracy

Equal panels are not aligned with the knots of a spline, so a family
whose derivative has kinks is integrated less accurately than a smooth
one; on a polynomial family the result is exact to rounding. A larger
`panels` reduces the error where the breaks cannot be aligned.

## See also

[`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md),
which builds the composite rule;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generic.
