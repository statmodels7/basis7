# Gram Matrix of a Basis

Returns the matrix of inner products of the `order`-th derivatives of
the basis functions, \$\$G\_{ab} = \int\_{a}^{b}
\varphi_a^{(d)}(t)\\\varphi_b^{(d)}(t)\\ \mathrm{d}t,\$\$ which is the
matrix of a roughness penalty: \\\beta^\top G_2 \beta\\ is exactly
\\\int (f'')^2\\ for the function \\f\\ the coefficients describe. The
three shipped families compute it in closed form, so no quadrature error
enters a penalized fit.

## Usage

``` r
basis_gram(basis, order = 0L, at = NULL, weight = NULL, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- order:

  What the inner products are taken of. A single non-negative whole
  number, default `0`, giving the inner products of the basis functions
  themselves; `2` is the usual roughness penalty; one entry per variable
  for a basis of several. A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md)
  instead gives \\\int (Lb)(Lb)^\top\\, the roughness matrix of that
  operator, and routes to
  [`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md).

- at:

  An optional numeric vector of points, or a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns. When given, the inner products are taken against the
  empirical measure of those points and divided by their number. Missing
  values are dropped before evaluation; no usable point left throws.

- weight:

  An optional function of one numeric vector returning one non-negative
  value per point, a density to weight the integral by. Refused for a
  basis of several variables. A function returning the wrong length, an
  `NA` or a negative value throws.

- ...:

  Passed to methods, and to
  [`weighted_gram()`](https://statmodels7.github.io/basis7/reference/weighted_gram.md)
  when `weight` is given, where `panels` and `nodes` control the
  quadrature.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins. Positive semidefinite, and singular for any
`order >= 1`.

## What it is and what it is not

The Gram matrix is an inner product of basis functions, so it depends on
the basis and the interval and on nothing else. It says which
combinations of coefficients are wiggly; how hard to shrink them is a
separate decision belonging to whatever fits the model.

It is symmetric and positive semidefinite by construction. At
`order = 0` it is positive definite for a basis of linearly independent
functions. At any `order >= 1` it is singular, the constant
differentiating to zero, and its null space has dimension `order` for a
family that contains the polynomials of that degree.

## Three measures

The default is Lebesgue measure on the basis interval, computed from a
closed form by every shipped family.

`at` replaces it with the empirical measure of the points given,
\\B^{(d)\top} B^{(d)} / n\\. That is the matrix a design matrix
produces, and the one to diagonalize against when the construction
should depend on where the data lie;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
uses it.

`weight` takes a weighted Lebesgue measure \\\int B B^\top w\\, by
composite Gauss-Legendre over 50 panels of 12 nodes. A weight is an
arbitrary function, so no family has a closed form and the quadrature is
always run. Measured against the closed forms with \\w \equiv 1\\: 5e-15
for Fourier and Legendre at order 0, and 2e-11 at order 2; for a cubic
B-spline 5e-12 at order 0 and 2e-6 relative at order 2, where the second
derivative has kinks at knots that the panel breaks do not line up with.
Give at most one of `at` and `weight`; both together throws.

## Where the arguments are handled

`at` and `weight` are dealt with in the body of the generic, before
dispatch, so a method never sees either and returns the plain Lebesgue
matrix alone. A method must still carry both names in its signature, S7
requiring a method's formals to contain the generic's.

## See also

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
whose columns it takes the inner products of;
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
which makes the `order = 0` matrix the identity;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
which diagonalizes one against the other.

## Examples

``` r
b <- fourier_basis(dimension = 5)

# Diagonal at order 0, the Fourier functions being orthogonal.
round(basis_gram(b), 6)
#>       const sin1 cos1 sin2 cos2
#> const     1  0.0  0.0  0.0  0.0
#> sin1      0  0.5  0.0  0.0  0.0
#> cos1      0  0.0  0.5  0.0  0.0
#> sin2      0  0.0  0.0  0.5  0.0
#> cos2      0  0.0  0.0  0.0  0.5
round(basis_gram(b, order = 1), 4)
#>       const    sin1    cos1    sin2    cos2
#> const     0  0.0000  0.0000  0.0000  0.0000
#> sin1      0 19.7392  0.0000  0.0000  0.0000
#> cos1      0  0.0000 19.7392  0.0000  0.0000
#> sin2      0  0.0000  0.0000 78.9568  0.0000
#> cos2      0  0.0000  0.0000  0.0000 78.9568

# The order-2 matrix is the penalty matrix: beta' G beta is the integrated
# squared second derivative, which a fine trapezoid rule confirms to 6e-11.
p <- poly_basis(dimension = 5)
beta <- c(0.2, 1.1, -0.4, 0.8, 0.1)
tt <- seq(0, 1, length.out = 200001)
fpp <- drop(basis_deriv(p, tt, order = 2) %*% beta)^2
c(quadratic_form = drop(t(beta) %*% basis_gram(p, order = 2) %*% beta),
  quadrature = sum(head(fpp, -1) + tail(fpp, -1)) / 2 * (tt[2] - tt[1]))
#> quadratic_form     quadrature 
#>         807.84         807.84 

# Singular from order 1 up: order d leaves exactly d zero eigenvalues,
# the polynomials of degree below d differentiating away.
bs <- bspline_basis(dimension = 8)
vapply(0:3, function(d) {
  ev <- eigen(basis_gram(bs, order = d), only.values = TRUE)$values
  sum(abs(ev) < 1e-8 * max(abs(ev)))
}, integer(1))
#> [1] 0 1 2 3

# Against the empirical measure of a sample instead of the interval.
set.seed(1)
round(basis_gram(b, at = runif(2000)), 3)
#>        const   sin1   cos1   sin2   cos2
#> const  1.000  0.021  0.019 -0.021  0.011
#> sin1   0.021  0.494 -0.010  0.010  0.010
#> cos1   0.019 -0.010  0.506  0.031  0.009
#> sin2  -0.021  0.010  0.031  0.503 -0.010
#> cos2   0.011  0.010  0.009 -0.010  0.497

# And against a weighted Lebesgue measure.
round(basis_gram(b, weight = function(x) dbeta(x, 2, 5)), 4)
#>         const    sin1    cos1   sin2    cos2
#> const  1.0000  0.5999 -0.0670 0.1095 -0.1467
#> sin1   0.5999  0.5733  0.0547 0.0045 -0.2828
#> cos1  -0.0670  0.0547  0.4267 0.3171 -0.0714
#> sin2   0.1095  0.0045  0.3171 0.5224  0.0074
#> cos2  -0.1467 -0.2828 -0.0714 0.0074  0.4776
```
