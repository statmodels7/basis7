# Linearly Transformed Basis

The S7 class of bases obtained from another by a fixed linear map of its
functions, \\\tilde{B}(x) = B(x)\\T\\. It carries the parent and the
matrix, and it is itself a basis, so a constrained or rotated basis goes
on answering every generic and nothing downstream has to know a
transformation happened. Constructed by
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
or
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md).

## Usage

``` r
TransformedBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list(),
  parent_basis = basis(),
  transform = integer(0)
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

- parent_basis:

  The basis being transformed, any object inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md). Kept
  whole, so it can still be evaluated and asked what it is.

- transform:

  The matrix \\T\\, numeric, with one row per parent function and one
  column per function of the result.

## Value

An object of class `TransformedBasis`, inheriting from
[basis](https://statmodels7.github.io/basis7/reference/basis.md), with
the five properties of a basis plus `parent_basis` and `transform`.

## Three operations, one class

Orthonormalizing a basis, restricting it to satisfy a linear constraint,
and rebuilding it so that it diagonalizes an inner product are the same
operation with different matrices. Each generic follows from linearity:

\$\$\tilde{B}^{(d)}(x) = B^{(d)}(x)\\T, \qquad \int_a^x \tilde{B} =
\left(\int_a^x B\right) T, \qquad \tilde{G} = T^\top G\\ T.\$\$

Differentiation and integration are linear and \\T\\ does not depend on
\\x\\, so both transform by the same matrix; the Gram matrix transforms
by congruence. A parent with exact derivatives and an exact Gram matrix
therefore passes its exactness on, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports the parent's answer in place of this class's own methods.

The anchored integral survives too: a linear combination of columns that
are all zero at the lower endpoint is zero there.

## Fewer columns than rows

\\T\\ may be \\K \times m\\ with \\m \< K\\, which is how
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
reduces the dimension. `@dimension` is `ncol(transform)` and the
parent's is `nrow(transform)`.

## Transforms compose by multiplication

Transforming a `TransformedBasis` again produces one object holding the
product of the two matrices, with the original as its parent. A chain of
transformations therefore costs one matrix multiplication per evaluation
however long it is. The `@basis_name` still nests, so an orthonormalized
orthonormal basis reads `orthonorm(orthonorm(bspline))` while
`@parent_basis` is the B-spline.

## See also

[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
the three constructors;
[`new_transformed()`](https://statmodels7.github.io/basis7/reference/new_transformed.md),
which they share.

## Examples

``` r
o <- orthonorm_basis(bspline_basis(dimension = 6))
S7::S7_inherits(o, TransformedBasis)
#> [1] TRUE
dim(o@transform)
#> [1] 6 6

# A constraint takes a column away.
x <- seq(0, 1, length.out = 200)
b <- bspline_basis(dimension = 6)
dim(constrain_basis(b, colSums(basis_eval(b, x)))@transform)
#> [1] 6 5

# Transforming twice keeps one matrix and the original parent.
oo <- orthonorm_basis(o)
dim(oo@transform)
#> [1] 6 6
class(oo@parent_basis)
#> [1] "basis7::BsplineBasis" "basis7::basis"        "S7_object"           
```
