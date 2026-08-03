# Evaluate a Basis Against Coefficients

Returns the values of the function the coefficients describe, without
necessarily forming the design matrix.

## Usage

``` r
basis_contract(basis, x, coef, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  Evaluation points: a numeric vector, or a matrix with one column per
  variable.

- coef:

  The coefficients: a vector, an array with one dimension per marginal,
  or a list of factor matrices.

- ...:

  Passed to methods.

## Value

A numeric vector with one value per evaluation point, or a matrix with
one column per column of `coef` when several sets are given.

## Details

For an ordinary basis this is `basis_eval(basis, x) %*% coef` and there
is nothing to save. For a tensor product there is: the design matrix has
\\\prod_j K_j\\ columns, so forming it is what makes a model with
several variables expensive, while the value it is used to compute needs
only the marginal evaluations.

Coefficients come in two shapes.

- An **array** of dimension \\(K_1, \ldots, K_D)\\, which is the general
  case. The rows are processed in blocks, so the peak memory is bounded
  by the block size rather than by the number of observations, however
  large the product.

- A **list of factor matrices** \\\Gamma_j\\ of size \\K_j \times F\\,
  the canonical polyadic form, in which the coefficient array is a sum
  of \\F\\ outer products. Here the value is \\\sum_f \prod_j
  B_j(x_j)^\top \gamma\_{j,f}\\, which costs \\O(nF\sum_j K_j)\\ in both
  time and memory: neither the design matrix nor the coefficient array
  is ever formed.

The second shape is what makes a model with high-order interactions
affordable, and what the factorised tensor product spline models of
Ruegamer (2024) estimate. Choosing the factors is a modelling decision
and belongs to the layer that owns the parameters; evaluating them is
basis arithmetic and belongs here.

## References

Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
*Proceedings of AISTATS*.

## See also

[`tensor_basis`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
[`basis_eval`](https://statmodels7.github.io/basis7/reference/basis_eval.md)

## Examples

``` r
b <- tensor_basis(bspline_basis(dimension = 5), bspline_basis(dimension = 4))
x <- cbind(runif(10), runif(10))

# a full array of coefficients
set.seed(1)
cf <- array(rnorm(20), dim = c(5, 4))
max(abs(basis_contract(b, x, cf) - basis_eval(b, x) %*% as.numeric(cf)))
#> [1] 0.9842723

# or a rank-two factorisation, which never forms either matrix
g <- list(matrix(rnorm(10), 5, 2), matrix(rnorm(8), 4, 2))
head(basis_contract(b, x, g))
#> [1]  0.08658967 -0.37591495  0.27853745  0.52864812  0.35983032  0.70173222
```
