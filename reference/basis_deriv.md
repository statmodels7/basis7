# Differentiate a Basis

Returns the `order`-th derivative of every basis function at the given
points, as a matrix of the same shape
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
returns. Since an expansion is linear in its coefficients, this matrix
times \\\beta\\ is the `order`-th derivative of the fitted function, so
a derivative of a fit needs no refitting.

## Usage

``` r
basis_deriv(basis, x, order = 1L, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  Evaluation points inside the basis interval: a numeric vector for a
  basis of one variable, or a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns for a basis of several.

- order:

  The derivative order. A single non-negative whole number for a basis
  of one variable, default `1`; a vector of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  such numbers for a basis of several, or a single `0`. A negative,
  fractional or missing order throws.

- ...:

  Passed to methods. No shipped family and no fallback reads anything
  from it.

## Value

A numeric matrix of `n` rows and `basis@dimension` columns, with column
names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md).
All zero above the smoothness of the family.

## Order is an argument

Derivative order is an argument, not a family of generics, because it is
unbounded: a Fourier basis is differentiable to any order, and a spline
of degree \\k\\ has \\k\\ non-trivial derivatives and zeros above that.
An order beyond what the family carries returns the zero matrix, which
is the value of the derivative; nothing is thrown, and the zeros are
exact.

`order = 0` short-circuits to
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
in the generic body, so a loop over orders needs no special case and
pays nothing for the zero.

## A basis of several variables takes a multi-index

For a product basis `order` has one entry per variable and names a mixed
partial: `c(2, 0)` is \\\partial^2/\partial x_1^2\\ and `c(1, 1)` is
\\\partial^2/\partial x_1 \partial x_2\\. A single non-zero number is
refused, having two readings; `0` alone is accepted, meaning no
derivative under either. See
[`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md).

## The numerical fallback

A subclass registering no method gets the one on the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class,
which applies **one** stencil of the order asked for to
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
never a composition of lower-order differences. The offsets, weights and
step come from
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html),
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
and
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
and the stencil is shifted to one side near an endpoint so that no node
leaves the interval.
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
says whether this is the route in use, and
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
measures the agreement.

For a basis of several variables the fallback differentiates one
coordinate at a time, so a mixed partial such as `c(1, 1)` throws there:
a stencil in the plane has the product of two errors, and the one family
that needs mixed partials,
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
computes them exactly from its margins.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
which is `order = 0`;
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the opposite direction;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
to learn whether a family answers this from a formula or a stencil.

## Examples

``` r
b <- bspline_basis(dimension = 6, degree = 3)
round(basis_deriv(b, c(0.25, 0.5, 0.75), order = 1), 3)
#>         bs1    bs2    bs3    bs4   bs5   bs6
#> [1,] -0.562 -2.391  2.109  0.844 0.000 0.000
#> [2,]  0.000 -0.562 -1.687  1.687 0.562 0.000
#> [3,]  0.000  0.000 -0.844 -2.109 2.391 0.563

# A cubic spline has three derivatives and then exact zeros.
all(basis_deriv(b, 0.5, order = 4) == 0)
#> [1] TRUE

# Order 0 is the evaluation itself.
identical(basis_deriv(b, 0.5, order = 0), basis_eval(b, 0.5))
#> [1] TRUE

# The derivative of a fit is the derivative of the basis times the same
# coefficients, so no refitting is involved.
beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
h <- 1e-5
fp <- (drop(basis_eval(b, 0.4 + h) %*% beta) -
       drop(basis_eval(b, 0.4 - h) %*% beta)) / (2 * h)
c(exact = drop(basis_deriv(b, 0.4, order = 1) %*% beta), difference = fp)
#>      exact difference 
#>      0.153      0.153 

# On a product basis the order is a multi-index, one entry per variable.
tb <- tensor_basis(bspline_basis(dimension = 4), poly_basis(dimension = 3))
round(basis_deriv(tb, cbind(0.5, 0.5), order = c(1, 1)), 4)
#>      bs1.P0 bs1.P1 bs1.P2 bs2.P0 bs2.P1 bs2.P2 bs3.P0 bs3.P1 bs3.P2 bs4.P0
#> [1,]      0   -1.5      0      0   -1.5      0      0    1.5      0      0
#>      bs4.P1 bs4.P2
#> [1,]    1.5      0
try(basis_deriv(tb, cbind(0.5, 0.5), order = 1))
#> Error : 'order' must have one entry per variable (2), or be 0. A single non-zero order is ambiguous for a basis of several variables: it could mean that order in each coordinate, or that total order.
```
