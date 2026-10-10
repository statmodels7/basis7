# Print a Basis

Prints a summary of a basis object: the family it comes from, how many
functions it holds and over how many variables, the interval each
variable runs over, the parameters the family carries, and which of the
three derived quantities are computed numerically instead of from a
formula. The summary has four lines, plus a header and one line per
parameter when the basis carries parameters.

## Arguments

- x:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- ...:

  Unused, and accepted so that the signature matches the
  [`print()`](https://rdrr.io/r/base/print.html) generic's.

## Value

`x`, invisibly.

## The lines

`Basis:` is `@basis_name`, which records how the object was built. A
wrapper or a product names its parents, so an orthonormalized B-spline
reads `orthonorm(bspline)` and a two-way product reads
`tensor(bspline, fourier)`.

`Functions:` is `@dimension`, the number of columns
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
returns, and `Variables:` is
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md),
the number of columns of `x` it expects. Every family here has one
variable except
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md).

`Domain:` is `@lower` and `@upper` as a closed interval per variable,
separated by `x` for a product. It is the interval every generic checks
its evaluation points against. A point outside signals an error that
states how many of the points were outside and what the interval is.

`Parameters:` appears only when `@basis_params` is non-empty, and prints
one line per entry. A numeric entry of more than four values is
abbreviated to its length, so a B-spline over many knots reads
`<8 values>` and leaves the knots to `b@basis_params$knots`.

`Numerical:` reads
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
and names those of `basis_deriv`, `basis_int` and `basis_gram` that are
computed by the numerical route: finite differences for the derivatives,
quadrature for the integral and the Gram matrix. The three shipped
families report `none`, except a Fourier basis whose interval is not a
whole period, whose Gram matrix is computed by quadrature. A basis
defined from its evaluation alone reports all three.

## See also

[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
for the last line as a named logical vector,
[`plot.basis()`](https://statmodels7.github.io/basis7/reference/plot.basis.md)
to see the functions themselves, and
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
for a verification of the components this summary only names.

## Examples

``` r
# The parameter block differs by family: the degree and the knots for a
# B-spline, the period, the number of pairs and a full-period flag for
# Fourier, the degree alone for Legendre.
bspline_basis(dimension = 6)
#> Basis: bspline
#> Functions: 6   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree          3
#>   knots           0.3333, 0.6667
#>   boundary_knots  0, 1
#> Numerical: none
fourier_basis(dimension = 5)
#> Basis: fourier
#> Functions: 5   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   omega        1
#>   n_pairs      2
#>   full_period  TRUE
#> Numerical: none
poly_basis(dimension = 4)
#> Basis: legendre
#> Functions: 4   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree  3
#> Numerical: none

# A product names its margins and reports one interval per variable.
tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#> Basis: tensor(bspline, fourier)
#> Functions: 12   Variables: 2
#> Domain: [0, 1] x [0, 1]
#> Parameters:
#>   marginal_dimensions  4, 3
#> Numerical: none

# More than four values in a parameter are abbreviated to a count. The
# values themselves stay reachable on the object.
b <- bspline_basis(dimension = 12)
b
#> Basis: bspline
#> Functions: 12   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree          3
#>   knots           <8 values>
#>   boundary_knots  0, 1
#> Numerical: none
b@basis_params$knots
#> [1] 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667 0.7777778
#> [8] 0.8888889

# A basis defined from its evaluation alone reports all three quantities
# as numerical, where every shipped family reports none.
Bumps <- S7::new_class("Bumps", parent = basis)
S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
  out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
                          "-")^2 / 0.12^2)
  colnames(out) <- basis_colnames(basis)
  out
}
#> Overwriting method basis_eval(<Bumps>)
Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)
#> Basis: bumps
#> Functions: 4   Variables: 1
#> Domain: [0, 1]
#> Numerical: basis_deriv, basis_int, basis_gram
```
