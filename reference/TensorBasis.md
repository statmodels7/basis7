# Tensor Product Basis

The basis of products of the functions of several bases, one per
variable: \$\$B(x_1, \ldots, x_D) = B_1(x_1) \otimes \cdots \otimes
B_D(x_D).\$\$ Constructed by
[`tensor_basis`](https://statmodels7.github.io/basis7/reference/tensor_basis.md).

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

  A short name for the family, used when printing.

- dimension:

  The number of functions in the basis, a positive integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several.

- basis_params:

  A named list of whatever else the subclass needs.

- marginals:

  The list of bases being multiplied, one per variable.

## Value

An object of class `TensorBasis`.

## Details

Everything a tensor product needs follows from the marginals, because
the product separates. A partial derivative differentiates one marginal
and leaves the others alone; the integral over the box from its lower
corner is the product of the marginal integrals; and the Gram matrix is
the Kronecker product of the marginal Gram matrices, so a tensor of
exactly integrated marginals is exactly integrated too, however many
variables it has.

Columns follow the convention of
[`kronecker`](https://rdrr.io/r/base/kronecker.html): the last marginal
varies fastest.

## See also

[`tensor_basis`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
[`basis_contract`](https://statmodels7.github.io/basis7/reference/basis_contract.md)

## Examples

``` r
t2 <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
c(basis_nvar(t2), t2@dimension)
#> [1]  2 12
```
