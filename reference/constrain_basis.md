# Restrict a Basis to a Linear Constraint

Returns a basis whose functions are exactly those of the parent
satisfying \\C\beta = 0\\, with the dimension reduced by the rank of
\\C\\. This is how a basis that carries its own constant is made to sit
beside an intercept, and how a smooth term is made identifiable.

## Usage

``` r
constrain_basis(basis, constraint, tol = 1e-10)
```

## Arguments

- basis:

  The basis to restrict, any object inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- constraint:

  A numeric matrix with one column per basis function, or a plain vector
  for a single constraint, which is taken as one row. Rows need not be
  independent; the rank is computed.

- tol:

  The relative tolerance below which a singular value counts as zero
  when the rank is determined, default `1e-10`, applied as
  `tol * max(d)`.

## Value

An object of class
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
of dimension `basis@dimension - rank`, with `basis_name`
`constrained(<parent>)`, column names `cn1`, `cn2`, and so on, and
`basis_params$constraint_rank`.

## The construction

The transform is an orthonormal basis of the null space of \\C\\, taken
from the right singular vectors of its singular value decomposition
belonging to the zero singular values. The constrained basis therefore
spans precisely the admissible functions, and rescaling or reordering
the rows of \\C\\ changes nothing about the space produced.

The rank is counted as the number of singular values above
`tol * max(d)`, and the result has `basis@dimension - rank` columns,
with the rank recorded in `basis_params$constraint_rank`.

## Writing the constraint

`constraint` has one column per basis function, so a condition on the
fitted curve is expressed by first evaluating the basis. The usual
sum-to-zero identifiability constraint over a grid or over the observed
covariate is `colSums(basis_eval(b, x))`, which makes the fitted values
sum to zero there; several constraints are the rows of a matrix.

What the package supplies is the mechanics. Which constraint a model
term should carry, whether a sum-to-zero condition for identifiability
or orthogonality to a linear part, needs to know what the term means and
belongs to the layer that does.

## Errors

A `constraint` that is not numeric, or whose column count is not
`basis@dimension`, throws with both numbers named; a missing value
throws; and a constraint of full rank leaves no functions and throws, in
place of returning a basis of zero columns.

## See also

[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
for a transformation that keeps the dimension;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
which applies a constraint and a rotation together.

## Examples

``` r
b <- bspline_basis(dimension = 6)
x <- seq(0, 1, length.out = 200)

# Sum to zero over a grid: the usual identifiability constraint, which
# removes the constant a B-spline basis carries.
cs <- constrain_basis(b, colSums(basis_eval(b, x)))
cs@dimension
#> [1] 5
cs@basis_params$constraint_rank
#> [1] 1
max(abs(colSums(basis_eval(cs, x))))
#> [1] 2.034527e-14

# Two constraints take two columns: sum to zero and orthogonal to x.
C <- rbind(colSums(basis_eval(b, x)), colSums(basis_eval(b, x) * x))
constrain_basis(b, C)@dimension
#> [1] 4

# A constraint of full rank leaves nothing, and is refused.
try(constrain_basis(b, diag(6)))
#> Error : The constraint leaves no functions: its rank equals the dimension of the basis.
```
