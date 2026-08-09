# Integrate a Basis

Returns the definite integral of every basis function from the lower
endpoint of the basis interval up to each evaluation point.

## Usage

``` r
basis_int(basis, x, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Passed to methods.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

Column \\j\\ of the result is

\$\$\int\_{a}^{x} \phi_j(t)\\\mathrm{d}t,\$\$

with \\a\\ the basis interval's lower endpoint, so that the integral of
an expansion is the expansion of the integral: \\\int_a^x \sum_j \beta_j
\phi_j = \sum_j \beta_j\\ times column \\j\\.

The convention is fixed and is part of the contract: the value at
`basis@lower` is exactly zero, for every basis and every column. Any
antiderivative would satisfy the differentiation check, so without a
fixed constant of integration two bases could disagree while both being
right, and a sum of them would be wrong in a way nothing would report.

## See also

[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_deriv`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)

## Examples

``` r
b <- fourier_basis(dimension = 5)
basis_int(b, b@lower) # exactly zero, by the convention above
#>      const sin1 cos1 sin2 cos2
#> [1,]     0    0    0    0    0
round(basis_int(b, c(0.25, 0.5, 1)), 4)
#>      const   sin1   cos1   sin2 cos2
#> [1,]  0.25 0.1592 0.1592 0.1592    0
#> [2,]  0.50 0.3183 0.0000 0.0000    0
#> [3,]  1.00 0.0000 0.0000 0.0000    0
```
