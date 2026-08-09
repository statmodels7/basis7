# Column Names of a Basis Matrix

The names every matrix the basis produces carries. The default numbers
the functions after the family name; a subclass whose functions have
their own identities overrides it.

## Usage

``` r
basis_colnames(basis, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- ...:

  Passed to methods.

## Value

A character vector of length `basis@dimension`.

## See also

[`basis_nvar`](https://statmodels7.github.io/basis7/reference/basis_nvar.md),
[`basis_is_numerical`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)

## Examples

``` r
basis_colnames(bspline_basis(dimension = 4))
#> [1] "bs1" "bs2" "bs3" "bs4"
basis_colnames(fourier_basis(dimension = 5))
#> [1] "const" "sin1"  "cos1"  "sin2"  "cos2" 
```
