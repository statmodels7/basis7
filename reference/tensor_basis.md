# Construct a Tensor Product Basis

Multiplies bases, one per variable, into the basis of all products of
their functions. This is how a smooth surface of several covariates is
built from one-dimensional pieces: the result is a basis like any other,
evaluated at a matrix of points instead of a vector, and every generic
answers exactly where the margins do.

## Usage

``` r
tensor_basis(...)
```

## Arguments

- ...:

  The bases to multiply, two or more, or a single list of them. Each
  must take one variable; a
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  among them is flattened into its own margins.

## Value

An object of class
[TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
with `basis_name` `tensor(<names>)`, `basis_params` holding
`marginal_dimensions`, and column names pasting the margins' with dots.

## Size

The result has \\\prod_j K_j\\ functions and takes \\D\\ variables, so
the evaluation points become a matrix of \\D\\ columns. The dimension
grows geometrically: four cubic B-splines of eight functions each give
4096 columns, and the design matrix at 20000 observations would be 625
MB.

[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
exists for that reason. It computes what a fit needs from the marginal
evaluations, in blocks for a full coefficient array and without forming
anything at all for a factorized one.

## Flattening

Each marginal must take one variable. A product of products is flattened
and never nested, so `@marginals` always holds the original bases and
[`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
counts them all: multiplying a two-way product by a third basis gives
one three-way product with three margins.

## The interval

`@lower` and `@upper` hold one entry per variable, taken from the
margins, so the domain is the box they span and a point outside any
margin's interval throws.

## References

Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
generalized additive mixed models. *Biometrics* **62**, 1025-1036.

## See also

[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md),
for evaluating a fit without the design matrix;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
which is a Kronecker product here;
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
and
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
which apply to a product as to any basis.

## Examples

``` r
b <- tensor_basis(bspline_basis(dimension = 4),
                  bspline_basis(dimension = 3, degree = 2))
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

# The Gram matrix is the Kronecker product of the marginal ones, exactly.
max(abs(basis_gram(b) - kronecker(
  basis_gram(bspline_basis(dimension = 4)),
  basis_gram(bspline_basis(dimension = 3, degree = 2))
)))
#> [1] 0

# The integral is anchored at the lower corner of the box.
basis_int(b, cbind(0, 0))
#>      bs1.bs1 bs1.bs2 bs1.bs3 bs2.bs1 bs2.bs2 bs2.bs3 bs3.bs1 bs3.bs2 bs3.bs3
#> [1,]       0       0       0       0       0       0       0       0       0
#>      bs4.bs1 bs4.bs2 bs4.bs3
#> [1,]       0       0       0

# A product of products is flattened, so the margins stay one-dimensional.
t3 <- tensor_basis(b, fourier_basis(dimension = 3))
c(basis_nvar(t3), length(t3@marginals), t3@dimension)
#> [1]  3  3 36

# The dimension grows geometrically. basis_contract() exists for that.
vapply(2:5, function(D) {
  tensor_basis(rep(list(bspline_basis(dimension = 8)), D))@dimension
}, numeric(1))
#> [1]    64   512  4096 32768
```
