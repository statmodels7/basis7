# Basis Expansion

The abstract parent of every basis in the package. A basis is a finite
collection of functions on an interval, and the object carries what is
needed to evaluate that collection, differentiate it, integrate it, and
compute its inner products.

## Usage

``` r
basis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list()
)
```

## Arguments

- basis_name:

  A short name for the family, used when printing.

- dimension:

  The number of functions in the basis, a positive integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several.

- basis_params:

  A named list of whatever else the subclass needs.

## Value

An object inheriting from class `basis`.

## Details

A basis of dimension \\d\\ on \\\[a, b\]\\ is the collection
\\\varphi_1, \dots, \varphi_d\\, and the object exists so that a
function may be written as a linear combination of them,

\$\$f(x) = \sum\_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\beta, \qquad
B(x)\_{ij} = \varphi_j(x_i),\$\$

with \\B(x)\\ the \\n \times d\\ design matrix
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
returns. Fitting \\f\\ is then a linear problem in \\\beta\\ whatever
the family, which is what makes derivatives, the anchored integral and
the Gram matrix properties of the basis rather than of the fit.

Concrete bases are subclasses. Each implements at least
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md);
[`basis_deriv`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
fall back to numerical methods registered on this class, so a subclass
that implements only its evaluation is immediately complete. Every
closed form registered later takes over through dispatch, with no change
to calling code.

Bases are complete: a B-spline basis carries all its functions and its
rows sum to one. Restricting a basis, whether for identifiability or to
separate a linear from a nonlinear part, is a linear transformation of
it and belongs to the layer that owns that decision.

A basis lives on an interval, or, when it is a product of several, on a
box: `lower` and `upper` then have one entry per variable and
[`basis_nvar`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
reports how many. Everything else is unchanged, and a univariate basis
is the case of one variable rather than a separate kind of object.

## See also

[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`check_basis`](https://statmodels7.github.io/basis7/reference/check_basis.md),
[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)

## Examples

``` r
# `basis` is abstract; construct a concrete subclass
b <- bspline_basis(dimension = 6)
b@dimension
#> [1] 6
c(b@lower, b@upper)
#> [1] 0 1
```
