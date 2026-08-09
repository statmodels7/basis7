# Which of a Basis's Methods Are Numerical

Reports, for each of the three derived generics, whether the basis
supplies its own method or falls back to the numerical one.

## Usage

``` r
basis_is_numerical(basis)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

## Value

A named logical vector with elements `basis_deriv`, `basis_int` and
`basis_gram`, `TRUE` where the numerical fallback is in force.

## Details

A value computed by finite differences or by quadrature cannot be
checked against a numerical reference: the comparison would be the same
arithmetic twice, agreeing however wrong the basis is.
[`check_basis`](https://statmodels7.github.io/basis7/reference/check_basis.md)
uses this to report such an order as not checked rather than as passed,
which is the difference between a validator and a formality.

## See also

[`basis_colnames`](https://statmodels7.github.io/basis7/reference/basis_colnames.md),
[`basis_nvar`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)

## Examples

``` r
basis_is_numerical(bspline_basis(dimension = 5))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 
```
