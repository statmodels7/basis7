# Tensor Product Basis

The S7 class of tensor product bases: all products of the functions of
several bases, one basis per variable, \$\$B(x_1, \ldots, x_D) =
B_1(x_1) \otimes \cdots \otimes B_D(x_D).\$\$ It carries the marginals
whole, and every quantity it answers with is built from theirs, so a
product of exactly integrated marginals is itself exact at any number of
variables. Constructed by
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md).

## Usage

``` r
TensorBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list(),
  marginals = list()
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

- marginals:

  The list of bases being multiplied, one per variable, each of one
  variable itself. Kept whole, so each can still be evaluated and asked
  what it is.

## Value

An object of class `TensorBasis`, inheriting from
[basis](https://statmodels7.github.io/basis7/reference/basis.md), with
the five properties of a basis plus `marginals`, and `basis_params`
holding `marginal_dimensions`.

## Everything follows from the marginals

The product separates, so each generic reduces to its marginals':

- a partial derivative differentiates one marginal to its own order and
  leaves the others alone, `order` being a multi-index;

- the integral over the box from the lower corner is the product of the
  marginal integrals, and the anchor survives, a product whose every
  factor is zero at the corner being zero there;

- the Gram matrix is the Kronecker product of the marginal Gram
  matrices, which is a product of one-dimensional integrals and never a
  quadrature over the box.

That last point is what keeps the construction affordable. A quadrature
over a box of \\D\\ variables costs a node count exponential in \\D\\
and carries an error to match; a Kronecker product of exact marginal
matrices is exact, and costs one marginal Gram matrix per variable.

## Column order

Columns follow
[`base::kronecker()`](https://rdrr.io/r/base/kronecker.html): the
**last** marginal varies fastest. At two variables of 3 and 2 functions
the names are `bs1.P0`, `bs1.P1`, `bs2.P0`, `bs2.P1`, `bs3.P0`,
`bs3.P1`.

That order matters when coefficients are supplied as an array, R storing
an array with its *first* index fastest.
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
reverses the dimensions for you; a hand-written `as.numeric(coef)` does
not, and pairs every coefficient with the wrong function.

## What the validator enforces

`@marginals` must be non-empty and hold only bases, and `@dimension`
must equal the product of their dimensions.

## See also

[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
the constructor;
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md),
which evaluates a fit without forming the product;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the Kronecker identity.

## Examples

``` r
t2 <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
c(basis_nvar(t2), t2@dimension)
#> [1]  2 12

# The last marginal varies fastest, and the names record both margins.
basis_colnames(tensor_basis(bspline_basis(dimension = 3, degree = 2),
                            poly_basis(dimension = 2)))
#> [1] "bs1.P0" "bs1.P1" "bs2.P0" "bs2.P1" "bs3.P0" "bs3.P1"

# A product is exact wherever its margins are.
basis_is_numerical(t2)
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 
```
