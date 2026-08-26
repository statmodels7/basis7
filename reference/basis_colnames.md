# Column Names of a Basis Matrix

Returns the names every matrix the basis produces carries, so that the
evaluation, the derivatives of every order, the anchored integral and
the Gram matrix of one basis agree on their columns. The default numbers
the functions after the family; a family whose functions have identities
of their own overrides it.

## Usage

``` r
basis_colnames(basis, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- ...:

  Passed to methods. No shipped family reads anything from it.

## Value

A character vector of length `basis@dimension`.

## Details

The default method takes the first two characters of `@basis_name` and
appends `1` to `@dimension`, giving `bs1 ... bs6` for a B-spline.
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
and
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
override it with names carrying meaning: `const`, `sin1`, `cos1`, `sin2`
for the first, `P0`, `P1`, `P2` for the second.

A wrapper numbers its own columns under a short prefix, so an
orthonormalized basis reads `on1 ... on5`. A
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
pastes its margins' names, one term per pair, as `bs1.const`,
`bs1.sin1`, `bs2.const`, so a coefficient's name says which marginal
function it belongs to in each variable.

Overriding it is how a subclass gives its columns meaning. The method
must return exactly `basis@dimension` strings;
[`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md)
sets them without checking, so a shorter vector is recycled by R and
silently mislabels.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
whose columns these name, and
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
for the object's summary.

## Examples

``` r
# The default numbers the functions after the family name.
basis_colnames(bspline_basis(dimension = 4))
#> [1] "bs1" "bs2" "bs3" "bs4"

# Fourier and Legendre name theirs instead.
basis_colnames(fourier_basis(dimension = 5))
#> [1] "const" "sin1"  "cos1"  "sin2"  "cos2" 
basis_colnames(poly_basis(dimension = 4))
#> [1] "P0" "P1" "P2" "P3"

# Every matrix the basis produces carries the same names.
b <- bspline_basis(dimension = 4)
identical(colnames(basis_eval(b, 0.5)), basis_colnames(b))
#> [1] TRUE
identical(colnames(basis_gram(b, order = 2)), basis_colnames(b))
#> [1] TRUE
```
