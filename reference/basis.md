# Basis Expansion

The abstract parent of every basis in the package, and the class to
inherit from when writing one of your own. A basis is a finite
collection of functions on an interval; the object carries the interval,
how many functions there are, and whatever else the family needs, and
the generics evaluate that collection, differentiate it, integrate it
and take its inner products. The class is abstract, so `basis(...)`
throws; construct one of the concrete families or a subclass.

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

This class is abstract and cannot be constructed. A subclass constructed
from it is an S7 object with properties `basis_name` (character),
`dimension` (integer), `lower` and `upper` (numeric, one entry per
variable) and `basis_params` (list).

## What a basis is

A basis of dimension \\d\\ on \\\[a, b\]\\ is a collection \\\varphi_1,
\dots, \varphi_d\\, and the object exists so that a function may be
written as a linear combination of them,

\$\$f(x) = \sum\_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\beta, \qquad
B(x)\_{ij} = \varphi_j(x_i),\$\$

with \\B(x)\\ the \\n \times d\\ design matrix
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
returns. Fitting \\f\\ is a linear problem in \\\beta\\ whatever the
family, so the derivative, the anchored integral and the Gram matrix are
all properties of the basis alone and can be computed once, before any
data arrive.

## Writing a subclass

A concrete basis is a subclass of this one. It must implement
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md);
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
have numerical methods registered on this class, so a subclass supplying
its evaluation alone answers all four generics immediately. A closed
form registered later takes over through dispatch, with no change to
calling code, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports which of the three are still on the fallback. The vignette
[`vignette("defining-a-basis")`](https://statmodels7.github.io/basis7/articles/defining-a-basis.md)
works one through.

## Bases here are complete

A B-spline basis carries all its functions, so its rows sum to one and
it spans the constant. Restricting a basis, for identifiability or to
separate a linear part from a nonlinear one, is a linear transformation
of it:
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
each return a
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
that is itself a basis, so nothing downstream has to know a restriction
happened.

## One variable or several

A basis lives on an interval, or, when it is a product of several, on a
box. `@lower` and `@upper` then hold one endpoint per variable and
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
reports how many;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
takes a matrix of that many columns instead of a vector. Nothing else
changes, and a basis of one variable is the case \\d = 1\\ of the same
object.

## What the validator enforces

`@basis_name` must be one string, `@dimension` one integer of at least
1, and `@lower` and `@upper` must be finite, of equal length and
strictly ordered within each variable. A double `@dimension` is rejected
by S7's own property check before the validator runs, so pass `6L` or an
[`as.integer()`](https://rdrr.io/r/base/integer.html); the constructors
take a plain `6` and convert it.

## See also

[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
and
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
for the three concrete families;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the generics every basis answers;
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
to verify a subclass of your own.

## Examples

``` r
# The class is abstract, so it is a subclass that gets constructed.
b <- bspline_basis(dimension = 6)
b@dimension
#> [1] 6
c(b@lower, b@upper)
#> [1] 0 1
b@basis_params$degree
#> [1] 3

# Its five properties are the same five on every basis in the package.
S7::prop_names(b)
#> [1] "basis_name"   "dimension"    "lower"        "upper"        "basis_params"
S7::prop_names(fourier_basis(dimension = 5))
#> [1] "basis_name"   "dimension"    "lower"        "upper"        "basis_params"

# Inheriting from it and writing one method gives a working basis.
Bumps <- S7::new_class("Bumps", parent = basis)
S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
  out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
                          "-")^2 / 0.12^2)
  colnames(out) <- basis_colnames(basis)
  out
}
bump <- Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)

# The three unwritten generics answer anyway, from finite differences.
basis_is_numerical(bump)
#> basis_deriv   basis_int  basis_gram 
#>        TRUE        TRUE        TRUE 
round(basis_deriv(bump, c(0.25, 0.75), order = 1), 4)
#>         bu1     bu2     bu3   bu4
#> [1,] -1.982  4.5471  0.0697 0.000
#> [2,]  0.000 -0.0697 -4.5471 1.982

# An interval the wrong way round is caught at construction.
try(Bumps(basis_name = "b", dimension = 4L, lower = 1, upper = 0))
#> Error : <Bumps> object is invalid:
#> - every @lower must be strictly less than its @upper
```
