# How Many Variables a Basis Takes

One for an ordinary basis, and the number of marginals for a product of
several.

## Usage

``` r
basis_nvar(basis)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

## Value

A positive integer.

## Details

The number is read off the endpoints, which carry one entry per
variable, so a basis declares its input dimension by construction rather
than by saying so separately and possibly disagreeing.

## See also

[`basis_colnames`](https://statmodels7.github.io/basis7/reference/basis_colnames.md),
[`basis_is_numerical`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)

## Examples

``` r
basis_nvar(bspline_basis(dimension = 5))
#> [1] 1
basis_nvar(tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3)))
#> [1] 2
```
