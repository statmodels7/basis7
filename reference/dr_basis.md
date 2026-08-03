# Demmler-Reinsch Basis

Returns the basis that simultaneously diagonalises the empirical inner
product at the given points and the penalty, and is empirically
orthogonal to the constraint functions.

## Usage

``` r
dr_basis(basis, x, penalty = NULL, constraints = NULL, scale = TRUE)
```

## Arguments

- basis:

  The basis to transform.

- x:

  The points the empirical inner product is taken at, normally the
  observed values of the covariate.

- penalty:

  A square penalty matrix with one row per basis function. Defaults to
  the second-derivative Gram matrix, the integrated squared second
  derivative. A discrete difference penalty is passed explicitly, for
  instance `crossprod(diff(diag(k), differences = 2))`.

- constraints:

  A matrix with one row per constraint and one column per evaluation
  point, whose row space the result is made empirically orthogonal to.
  Defaults to a constant and `x`, which is the separation of a linear
  from a nonlinear effect.

- scale:

  Whether to rescale so that \\\mathrm{tr}(Z^\top Z/n) = 1\\, which puts
  the bases of different terms on a common scale.

## Value

An object of class
[`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md).

## Details

The construction has three steps. The constraint matrix \\C =
(\mathbf{1}, x)^\top B\\ is formed and the basis restricted to its null
space \\V_0\\, which makes the remaining functions empirically
orthogonal to a constant and to \\x\\. The pencil \\(V_0^\top (B^\top
B/n) V_0,\\ V_0^\top P V_0)\\ is then diagonalised, and the transform is
\\T = V_0 A\\. The resulting design matrix \\Z = B T\\ has \\Z^\top Z\\
diagonal, has \\T^\top P T\\ equal to the identity, and satisfies
\\(\mathbf{1}, x)^\top Z = 0\\ exactly.

Three properties make this construction usable where others are not. It
factorises only a \\q \times K\\ and a \\(K-q) \times (K-q)\\ matrix,
never anything of the size of the sample. It tolerates a rank-deficient
\\B\\, which equally spaced knots produce whenever the data leave a knot
span empty, because the matrix inverted is the penalty and not the
design. And the transform is kept, so prediction at new points is the
parent's evaluation multiplied by it, like any other transformed basis.

The last property is what the separation of a linear from a nonlinear
effect needs: a reparametrisation that does not satisfy \\(\mathbf{1},
x)^\top Z = 0\\ estimates the sum of the two correctly and the split
between them with bias.

## References

Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
smoothing. *Numerische Mathematik* 24, 375-382.

## See also

[`constrain_basis`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
[`basis_gram`](https://statmodels7.github.io/basis7/reference/basis_gram.md)

## Examples

``` r
set.seed(1)
x <- sort(runif(200))
d <- dr_basis(bspline_basis(dimension = 12), x)

z <- basis_eval(d, x)
round(max(abs(crossprod(z)[upper.tri(crossprod(z))])), 10) # diagonal
#> [1] 0
round(max(abs(crossprod(cbind(1, x), z))), 10) # orthogonal to 1 and x
#> [1] 0
```
