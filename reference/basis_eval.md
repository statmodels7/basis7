# Evaluate a Basis

Returns the design matrix of a basis at the given points: one row per
evaluation point, one column per basis function, entry \\(i, j)\\ equal
to \\\varphi_j(x_i)\\. This is the matrix a regression on the basis is
fitted against, and the one every other generic in the package is
defined in terms of.

## Usage

``` r
basis_eval(basis, x, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  Evaluation points inside the basis interval: a numeric vector for a
  basis of one variable, or a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns for a basis of several, where a plain vector is taken by row.
  `NA` is allowed and produces a missing row. A point outside the
  interval throws.

- ...:

  Passed to methods. No shipped family reads anything from it.

## Value

A numeric matrix of `n` rows and `basis@dimension` columns, `n` being
`length(x)` for one variable and `nrow(x)` for several, with column
names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md).

## The design matrix

An expansion with coefficients \\\beta\\ is evaluated as

\$\$f(x) = \sum\_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\\\beta,\$\$

so `basis_eval(b, x) %*% beta` is the fitted function at `x`, and the
matrix is the design block a linear model on the basis uses. Its columns
carry the names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
declares, and every other matrix the basis produces carries the same
ones in the same order.

## The one generic a basis must implement

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
all have numerical methods registered on the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class,
computed from this one, so a subclass supplying its evaluation alone
answers all four. Registering a closed form for any of the three later
takes over through dispatch;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports which are still on the fallback.

## What the generic does before dispatching

The generic body validates `x` through
[`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md)
and passes the validated version on, so every method, including one
written outside the package, gets the same guarantees without writing
them: a point outside the interval throws, a point that is an endpoint
up to a relative `1e-8` arrives clamped exactly onto that endpoint, and
a basis of several variables receives a matrix of
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
columns whatever shape the caller passed.

`NA` is neither checked nor clamped and flows through arithmetic, so a
missing evaluation point gives a row of `NA`.

## See also

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for its derivatives,
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for its anchored integral,
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for its inner products, and
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
to get \\B(x)\beta\\ without forming \\B(x)\\ for a product basis.

## Examples

``` r
b <- bspline_basis(dimension = 6)
round(basis_eval(b, c(0, 0.5, 1)), 4)
#>      bs1    bs2    bs3    bs4    bs5 bs6
#> [1,]   1 0.0000 0.0000 0.0000 0.0000   0
#> [2,]   0 0.0312 0.4687 0.4687 0.0312   0
#> [3,]   0 0.0000 0.0000 0.0000 0.0000   1

# A B-spline basis is a partition of unity: every row sums to one.
rowSums(basis_eval(b, seq(0, 1, length.out = 5)))
#> [1] 1 1 1 1 1

# Evaluating an expansion is one matrix product.
beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
drop(basis_eval(b, c(0.3, 0.7)) %*% beta)
#> [1] 0.242525 0.421275

# Outside the interval it throws; an endpoint up to rounding is clamped.
try(basis_eval(b, 1.5))
#> Error : 1 of 1 evaluation points fall outside the basis interval [0, 1].
all.equal(basis_eval(b, 1 + 5e-9), basis_eval(b, 1))
#> [1] TRUE

# A basis of several variables takes one column per variable.
tb <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
dim(basis_eval(tb, cbind(c(0.1, 0.5), c(0.2, 0.6))))
#> [1]  2 12
```
