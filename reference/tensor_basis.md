# Construct a Tensor Product Basis

Multiplies bases, one per variable, into the basis of all products of
their functions.

## Usage

``` r
tensor_basis(...)
```

## Arguments

- ...:

  The bases to multiply, or a single list of them.

## Value

An object of class
[`TensorBasis`](https://statmodels7.github.io/basis7/reference/TensorBasis.md).

## Details

The result has \\\prod_j K_j\\ functions and takes \\D\\ variables, so
the evaluation points become a matrix with one column per variable. That
growth is the reason
[`basis_contract`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
exists: it computes what a fit needs from the marginal evaluations
alone, without ever forming the product.

Each marginal must take one variable. A product of products is flattened
rather than nested, so the marginals of the result are always the
original bases.

## References

Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics* 62, 1025-1036.

## See also

[`basis_contract`](https://statmodels7.github.io/basis7/reference/basis_contract.md),
[`basis_gram`](https://statmodels7.github.io/basis7/reference/basis_gram.md)

## Examples

``` r
b <- tensor_basis(bspline_basis(dimension = 4), bspline_basis(dimension = 3, degree = 2))
b
#> Basis: tensor(bspline, bspline)
#> Functions: 12   Variables: 2
#> Domain: [0, 1] x [0, 1]
#> Parameters:
#>   marginal_dimensions  4, 3
#> Numerical: none

x <- cbind(c(0.2, 0.5), c(0.7, 0.1))
round(basis_eval(b, x), 4)
#>      bs1.bs1 bs1.bs2 bs1.bs3 bs2.bs1 bs2.bs2 bs2.bs3 bs3.bs1 bs3.bs2 bs3.bs3
#> [1,]  0.0461  0.2150  0.2509  0.0346  0.1613  0.1882  0.0086  0.0403  0.0470
#> [2,]  0.1013  0.0225  0.0013  0.3038  0.0675  0.0038  0.3038  0.0675  0.0038
#>      bs4.bs1 bs4.bs2 bs4.bs3
#> [1,]  0.0007  0.0034  0.0039
#> [2,]  0.1013  0.0225  0.0013

# the Gram matrix is the Kronecker product of the marginal ones
g <- basis_gram(b)
max(abs(g - kronecker(
  basis_gram(bspline_basis(dimension = 4)),
  basis_gram(bspline_basis(dimension = 3, degree = 2))
)))
#> [1] 0
```
