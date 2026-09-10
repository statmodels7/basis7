# The Roughness Matrix of an Operator on a B-Spline Basis

Exact for the length measure, integrating knot interval by knot interval
as
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md)
does, and with one guard in front of it: an operator of order above the
degree of the spline is rejected rather than integrated.

## Arguments

- basis:

  A [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- at:

  The covariate values, for the empirical measure, or `NULL`.

- weight:

  A density to integrate against, or `NULL`.

- ...:

  Passed on to the quadrature (`panels`, `nodes`).

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns.

## Why it is exact

\\Lb\\ is a linear combination of derivatives of a spline of degree
\\p\\, so it is piecewise polynomial of degree at most \\p\\ with the
same breakpoints, and \\(Lb)(Lb)^\top\\ is piecewise of degree at most
\\2p\\. Gauss-Legendre with \\p + 1\\ nodes on each knot interval is
exact to degree \\2p + 1\\, so nothing is approximated. The panels of
the base method's rule do not line up with the knots, where a derivative
of the spline jumps, and it is that misalignment rather than the node
count that costs the accuracy: measured at `dimension = 10`,
`degree = 3`, the base rule at 50 panels differs from the derivative
route by 5e-8 relative, and this one by 3e-15.

## The guard

The \\m\\-th derivative of a spline of degree \\p\\ is identically zero
for \\m \> p\\, so the leading term of \\Lb\\ vanishes and what would be
integrated is the operator with its highest derivative deleted. That is
a different penalty with a different null space, and nothing about the
result would say so.
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
returns an exact zero matrix in the same situation for a plain
derivative, where the answer is at least unmistakable; here it is not,
so the refusal names the degree to raise.

A spline does not contain the null space of a periodic operator exactly,
which is admissible and is not this guard's business: see the section on
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)'s
page.

## See also

[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
for the generic.

## Examples

``` r
b <- bspline_basis(lower = 0, upper = 1, dimension = 10, degree = 3)
dim(basis_gram(b, order = harmonic_operator(1)))
#> [1] 10 10

# an operator of order 5 on a cubic spline is refused, not truncated
try(basis_gram(b, order = harmonic_operator(1, harmonics = 2)))
#> Error : the operator has order 5 and the spline has degree 3, so D^5 of every
#>   basis function is zero and the penalty would be the operator with its leading
#>   term deleted. Raise 'degree' to at least 5.
```
