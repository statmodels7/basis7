# Gram Matrix Against a Weighted Lebesgue Measure

Computes \\\int_a^b B^{(d)}(t)\\ B^{(d)}(t)^\top w(t)\\\mathrm{d}t\\ by
composite Gauss-Legendre. A weight is an arbitrary function, so no
family has a closed form for it and the quadrature is always run. Called
from the body of
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
when `weight` is supplied.

## Usage

``` r
weighted_gram(basis, order, weight, panels = 50L, nodes = 12L, ...)
```

## Arguments

- basis:

  A basis object of one variable, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md). More
  than one variable throws.

- order:

  The derivative order, already validated by
  [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md).

- weight:

  A function of one numeric vector returning one non-negative value per
  point. A wrong length, an `NA` or a negative value throws.

- panels:

  The number of equal subintervals, default `50`.

- nodes:

  The number of Gauss-Legendre nodes per subinterval, default `12`.

- ...:

  Unused.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## Details

The interval is cut into `panels` equal pieces and an `nodes`-point
Gauss-Legendre rule is placed on each, which integrates a polynomial of
degree up to `2 * nodes - 1` exactly on every panel. The weight is
folded into the quadrature weights and the matrix formed as a
crossproduct of \\\sqrt{w_i}\\B^{(d)}(t_i)\\, which keeps the result
positive semidefinite whatever the weight does.

Measured against the closed forms at \\w \equiv 1\\ with the defaults,
worst absolute entry: 5e-15 at order 0 and 2e-11 at order 2 for Fourier
and Legendre; 5e-12 at order 0 and 1.4e-3 at order 2 for a cubic
B-spline over eight knots, 2e-6 of the matrix's own scale, the second
derivative there having kinks the panel breaks do not line up with.
Raise `panels` when a family's derivative is not smooth.

A basis of several variables is refused: the rule above is
one-dimensional.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
its only caller, and
[`empirical_gram()`](https://statmodels7.github.io/basis7/reference/empirical_gram.md)
for the other alternative measure.
