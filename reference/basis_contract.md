# Evaluate a Basis Against Coefficients

Returns the values of the function a coefficient vector describes,
\\B(x)\beta\\, computed without forming \\B(x)\\ in full where that is
possible. For an ordinary basis it is the design matrix times the
coefficients; for a
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
the design matrix has \\\prod_j K_j\\ columns, and avoiding it keeps a
model of several variables affordable.

## Usage

``` r
basis_contract(basis, x, coef, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  Evaluation points: a numeric vector for a basis of one variable, or a
  matrix with one column per variable.

- coef:

  The coefficients. For an ordinary basis, a numeric vector of length
  `basis@dimension`, or a matrix of several such columns. For a
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md),
  an array of dimension `(K_1, ..., K_D)`, a vector of `basis@dimension`
  values in the design's own column order, or a list of `D` factor
  matrices each with `F` columns. A wrong length or a wrong set of
  dimensions signals an error that states what is expected.

- ...:

  Passed to methods. The
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  method reads `block` from it, the number of rows processed at once,
  default `1024`.

## Value

A numeric vector with one value per evaluation point, or a matrix with
one column per column of `coef` when several sets are given to an
ordinary basis.

## Two shapes of coefficient

An **array** of dimension \\(K_1, \ldots, K_D)\\ is the general case.
The rows are processed in blocks of `block`, so the peak memory is that
of one block's design matrix instead of the whole one, and the values
are those of the design-matrix product.

A **list of factor matrices** \\\Gamma_j\\ of size \\K_j \times F\\ is
the canonical polyadic form, in which the coefficient array is a sum of
\\F\\ outer products. The value is then \\\sum_f \prod_j B_j(x_j)^\top
\gamma\_{j,f}\\, costing \\O(nF\sum_j K_j)\\ in both time and memory:
neither the design matrix nor the coefficient array is formed.

The second shape makes a model with high-order interactions affordable,
and is what the factorized tensor product spline models of Ruegamer
(2024) estimate. The choice of the factors is a modeling decision made
in the layer that owns the parameters; this function evaluates the
factors it receives.

## An array and a vector are read differently

A tensor design's columns run with the **last** margin fastest,
following [`base::kronecker()`](https://rdrr.io/r/base/kronecker.html),
and an R array is stored with its **first** index fastest. An array is
therefore transposed with [`aperm()`](https://rdrr.io/r/base/aperm.html)
before it is flattened, so that `coef[m1, ..., mD]` is the coefficient
of the column named for those marginal functions.

A plain **vector** has no dimensions to reverse and is taken in the
design's own column order. So a vector and an array holding the same
numbers in the same storage order describe different functions, and the
identity to check against is
`basis_eval(b, x) %*% as.numeric(aperm(coef))`. Flattening without the
[`aperm()`](https://rdrr.io/r/base/aperm.html) pairs most coefficients
with the wrong function and returns a finite value with no error, as the
example below shows.

## References

Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
*Proceedings of AISTATS*.

## See also

[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
the case this exists for;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
the design matrix it avoids forming.

## Examples

``` r
b <- tensor_basis(bspline_basis(dimension = 5), bspline_basis(dimension = 4))
set.seed(1)
x <- cbind(runif(10), runif(10))

# A full array of coefficients. The identity holds against the flattening
# that matches the column order, which reverses the array's dimensions.
cf <- array(rnorm(20), dim = c(5, 4))
max(abs(basis_contract(b, x, cf) -
        basis_eval(b, x) %*% as.numeric(aperm(cf))))
#> [1] 0

# Flattening without that gives a wrong value and no error.
max(abs(basis_contract(b, x, cf) - basis_eval(b, x) %*% as.numeric(cf)))
#> [1] 1.162171

# A rank-two factorization, which forms neither the design nor the array.
g <- list(matrix(rnorm(10), 5, 2), matrix(rnorm(8), 4, 2))
head(basis_contract(b, x, g))
#> [1]  0.12417711  0.04423014  0.14065631 -0.33696805  0.02249578 -0.25496420

# And agrees with the array it stands for.
full <- outer(g[[1]][, 1], g[[2]][, 1]) + outer(g[[1]][, 2], g[[2]][, 2])
max(abs(basis_contract(b, x, g) - basis_contract(b, x, full)))
#> [1] 1.665335e-16

# On an ordinary basis it is the design matrix times the coefficients.
p <- poly_basis(dimension = 4)
max(abs(basis_contract(p, c(0.2, 0.8), 1:4) -
        basis_eval(p, c(0.2, 0.8)) %*% 1:4))
#> [1] 0
```
