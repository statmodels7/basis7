# The Roughness Matrix of an Operator on a B-Spline Basis

Exact for the length measure, integrating knot interval by knot interval
as
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md)
does. An operator of order above the degree of the spline is rejected.

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

  Passed on to the quadrature (`panels`, `nodes`) when `weight` is
  given. Without `at` or `weight` the matrix is exact, and `panels` or
  `nodes` signals an error.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns.

## Exactness

\\Lb\\ is a linear combination of derivatives of a spline of degree
\\p\\, so it is piecewise polynomial of degree at most \\p\\ with the
same breakpoints, and \\(Lb)(Lb)^\top\\ is piecewise of degree at most
\\2p\\. Gauss-Legendre with \\p + 1\\ nodes on each knot interval is
exact to degree \\2p + 1\\, so nothing is approximated. The panels of
the rule of the base method are equally spaced and do not line up with
the knots, where a derivative of the spline jumps, so the base method is
approximate for a spline.

## Operators above the degree

The \\m\\-th derivative of a spline of degree \\p\\ is identically zero
for \\m \> p\\, so the leading term of \\Lb\\ vanishes and what would be
integrated is the operator with its highest derivative deleted. That is
a different penalty with a different null space, and the result would
give no sign of it. For a plain derivative
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
returns the zero matrix in the same situation, which cannot be mistaken
for a penalty. For an operator the method signals an error that names
the degree to raise.

A spline contains the null space of a periodic operator only
approximately. That is admissible, and this method does not test for it;
see
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md).

## See also

[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
for the generic.

## Examples

``` r
b <- bspline_basis(lower = 0, upper = 1, dimension = 10, degree = 3)
dim(basis_gram(b, order = harmonic_operator(1)))
#> [1] 10 10

# an operator of order 5 on a cubic spline signals an error
try(basis_gram(b, order = harmonic_operator(1, harmonics = 2)))
#> Error : the operator has order 5 and the spline has degree 3, so D^5 of every
#>   basis function is zero and the penalty would be the operator with its leading
#>   term deleted. Raise 'degree' to at least 5.
```
