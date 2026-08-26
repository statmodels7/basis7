# Integrate a Basis

Returns the definite integral of every basis function from the lower
endpoint of the basis interval up to each evaluation point, as a matrix
of the same shape
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
returns. Times a coefficient vector it gives the antiderivative of the
fitted function that vanishes at the lower endpoint.

## Usage

``` r
basis_int(basis, x, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  Evaluation points inside the basis interval, the upper limits of the
  integrals: a numeric vector for a basis of one variable, or a matrix
  of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns for a basis of several.

- ...:

  Passed to methods. The numerical fallback reads `nodes` from it, the
  number of Gauss-Legendre nodes per panel, default `12`.

## Value

A numeric matrix of `n` rows and `basis@dimension` columns, with column
names
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md).
The row at `basis@lower` is exactly zero.

## The anchored integral

Column \\j\\ of the result is

\$\$\int\_{a}^{x} \varphi_j(t)\\\mathrm{d}t,\$\$

with \\a\\ the lower endpoint of the basis interval, so the integral of
an expansion is the expansion against the same coefficients: \\\int_a^x
\sum_j \beta_j \varphi_j = \sum_j \beta_j I_j(x)\\. The integral over
the whole interval is the row at `basis@upper`, and the integral over
\\\[u, v\]\\ is the difference of two rows.

## Why the anchor is fixed

The value at `basis@lower` is exactly zero, for every basis and every
column, and that is part of the contract every implementation owes. Any
antiderivative satisfies the differentiation check, so with the constant
of integration left free two bases could disagree while both being
right, and a sum of them would be wrong with nothing to report it.

## A basis of several variables

The integral is taken over the box from the lower corner to the point,
one iterated integral per variable, so on two variables column \\j\\ is
\\\int\_{a_1}^{x_1}\int\_{a_2}^{x_2} \varphi_j\\.

## The numerical fallback

A subclass registering no method gets the one on the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class:
composite Gauss-Legendre from the lower endpoint to each point,
`nodes = 12` per panel by default. Exact for a polynomial integrand of
degree up to `2 * nodes - 1`, and accurate to the panel width elsewhere.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
and
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
for the other direction, and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for integrals of products of basis functions.

## Examples

``` r
b <- fourier_basis(dimension = 5)

# Zero at the lower endpoint, by the convention above.
basis_int(b, b@lower)
#>      const sin1 cos1 sin2 cos2
#> [1,]     0    0    0    0    0
round(basis_int(b, c(0.25, 0.5, 1)), 4)
#>      const   sin1   cos1   sin2 cos2
#> [1,]  0.25 0.1592 0.1592 0.1592    0
#> [2,]  0.50 0.3183 0.0000 0.0000    0
#> [3,]  1.00 0.0000 0.0000 0.0000    0

# Differentiating the integral returns the basis.
bs <- bspline_basis(dimension = 6)
x <- c(0.2, 0.55, 0.9)
h <- 1e-6
fd <- (basis_int(bs, x + h) - basis_int(bs, x - h)) / (2 * h)
max(abs(fd - basis_eval(bs, x)))
#> [1] 1.532541e-11

# The integral of a fitted curve, and the integral over a subinterval.
beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
drop(basis_int(bs, 1) %*% beta)
#> [1] 0.2666667
drop((basis_int(bs, 0.75) - basis_int(bs, 0.25)) %*% beta)
#> [1] 0.1399089
```
