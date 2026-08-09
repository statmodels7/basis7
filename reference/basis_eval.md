# Evaluate a Basis

Returns the design matrix of the basis at the given points: one row per
evaluation point, one column per basis function.

## Usage

``` r
basis_eval(basis, x, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Passed to methods.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with the column names the basis declares.

## Details

Entry \\(i, j)\\ of the result is \\\phi_j(x_i)\\, so an expansion with
coefficients \\\beta\\ is evaluated as

\$\$f(x) = \sum\_{j=1}^{d} \beta_j \phi_j(x) = B(x)\\\beta,\$\$

and the matrix is the design block a regression on the basis uses.

This is the only generic a basis must implement. Everything else in the
package has a numerical method registered on the
[`basis`](https://statmodels7.github.io/basis7/reference/basis.md) class
and is therefore available from this one alone.

The generic validates the evaluation points before dispatching, so every
method, including one written outside the package, rejects a point
outside the basis interval and receives points that are endpoints up to
rounding already clamped onto them.

## See also

[`basis_deriv`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int`](https://statmodels7.github.io/basis7/reference/basis_int.md),
[`basis_gram`](https://statmodels7.github.io/basis7/reference/basis_gram.md)

## Examples

``` r
b <- bspline_basis(dimension = 6)
round(basis_eval(b, c(0, 0.5, 1)), 4)
#>      bs1    bs2    bs3    bs4    bs5 bs6
#> [1,]   1 0.0000 0.0000 0.0000 0.0000   0
#> [2,]   0 0.0312 0.4687 0.4687 0.0312   0
#> [3,]   0 0.0000 0.0000 0.0000 0.0000   1

# a B-spline basis is a partition of unity
rowSums(basis_eval(b, seq(0, 1, length.out = 5)))
#> [1] 1 1 1 1 1
```
