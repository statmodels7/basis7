# B-Spline Basis

The S7 class of B-spline bases, the objects
[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
returns. It adds no property to
[basis](https://statmodels7.github.io/basis7/reference/basis.md) and
exists as the class the spline methods dispatch on. A B-spline basis is
piecewise polynomial, each function supported on a few knot intervals,
and its functions sum to one.

## Usage

``` r
BsplineBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list()
)
```

## Arguments

- basis_name:

  A single string naming the family, printed by
  [`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  and used by wrappers to build their own name. Not read by any
  computation.

- dimension:

  The number of functions in the basis, a single integer of at least 1.
  Must be of storage mode integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several. Both finite, of the same length, and
  `lower[j] < upper[j]` for every `j`. Their length is what
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  reports.

- basis_params:

  A named list of whatever else the subclass needs: the knots and degree
  of a B-spline, the frequency of a Fourier basis, the marginal
  dimensions of a product.
  [`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  shows it, abbreviating any numeric entry of more than four values.
  Defaults to an empty list.

## Value

An object of class `BsplineBasis`, inheriting from
[basis](https://statmodels7.github.io/basis7/reference/basis.md), with
the same five properties and `basis_params` holding `degree`, `knots`
and `boundary_knots`. Call
[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
instead of this class directly; it places the knots and checks the
dimension against the degree.

## The Cox-de Boor recurrence

On a knot sequence \\t_1 \le \cdots \le t\_{d+m+1}\\ the functions are
defined from the indicators upward:

\$\$B\_{j,0}(x) = \mathbf{1}\\t_j \le x \< t\_{j+1}\\,\$\$

\$\$B\_{j,q}(x) = \frac{x - t_j}{t\_{j+q} - t_j} B\_{j,q-1}(x) +
\frac{t\_{j+q+1} - x}{t\_{j+q+1} - t\_{j+1}} B\_{j+1,q-1}(x), \qquad q =
1, \dots, m,\$\$

a term with a zero denominator being taken as zero.

## The two properties that follow

\\B\_{j,m}\\ vanishes outside \\\[t_j, t\_{j+m+1}\]\\, so at most \\m +
1\\ columns are non-zero in any row and the design matrix is banded: at
`degree = 3` exactly four of them, whatever the dimension. And \\\sum_j
B\_{j,m}(x) = 1\\ on the interval, measured to 2.2e-16, so the basis
carries its own constant and is collinear with an intercept in the same
design.

## Where the numbers come from

Evaluation, derivatives and integrals come from
[`splines2::bSpline()`](https://wwenjie.org/splines2/reference/bSpline.html),
through the single wrapper
[`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md),
which computes all three from the recurrence and takes no differences.
The Gram matrix is integrated here instead, knot interval by knot
interval, and is exact: see
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md).

## Its `basis_params`

Three entries: `degree`, the interior `knots` as a numeric vector, and
`boundary_knots`, which is `c(lower, upper)`.

## References

de Boor, C. (2001). *A Practical Guide to Splines*, revised edition.
Springer.

## See also

[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
the constructor;
[`basis_eval.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.BsplineBasis.md),
[`basis_deriv.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.BsplineBasis.md),
[`basis_int.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.BsplineBasis.md)
and
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md)
for the methods registered on it.

## Examples

``` r
b <- bspline_basis(dimension = 5)
S7::S7_inherits(b, BsplineBasis)
#> [1] TRUE
b@basis_params
#> $degree
#> [1] 3
#> 
#> $knots
#> [1] 0.5
#> 
#> $boundary_knots
#> [1] 0 1
#> 

# Local support: at most degree + 1 columns are non-zero in a row.
E <- basis_eval(bspline_basis(dimension = 12), seq(0.05, 0.95, by = 0.1))
rowSums(E != 0)
#>  [1] 4 4 4 4 4 4 4 4 4 4

# And the rows sum to one.
max(abs(rowSums(E) - 1))
#> [1] 2.220446e-16
```
