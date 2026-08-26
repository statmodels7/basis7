# Demmler-Reinsch Basis

Returns the basis that simultaneously diagonalizes the empirical inner
product at the given points and a penalty, and is empirically orthogonal
to a constant and to `x`. Its columns are ordered from the smoothest to
the wiggliest, each carrying a known share of the empirical variance, so
a penalized fit against it is a shrinkage of independent coordinates and
the linear part of the effect is separated from the nonlinear part
exactly.

## Usage

``` r
dr_basis(basis, x, penalty = NULL, constraints = NULL, scale = TRUE)
```

## Arguments

- basis:

  The basis to transform, any object inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  The points the empirical inner product is taken at, normally the
  observed covariate. A numeric vector inside the basis interval, with
  no missing values. The construction depends on where they lie, so a
  basis built for one sample is not the basis for another.

- penalty:

  A square numeric penalty matrix with one row and column per basis
  function. `NULL`, the default, uses `basis_gram(basis, order = 2)`,
  the integrated squared second derivative. A discrete difference
  penalty is passed explicitly, for instance
  `crossprod(diff(diag(k), differences = 2))`.

- constraints:

  A matrix with one row per constraint and one column per element of
  `x`, whose row space the result is made empirically orthogonal to.
  `NULL`, the default, uses `t(cbind(1, x))`, which is the separation of
  a linear from a nonlinear effect.

- scale:

  Whether to rescale so that \\\mathrm{tr}(Z^\top Z/n) = 1\\, default
  `TRUE`. See above for what it does to \\T^\top P T\\.

## Value

An object of class
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
of dimension `basis@dimension - rank(C)`, with `basis_name`
`dr(<parent>)`, column names `dr1`, `dr2`, and so on ordered from
smoothest to wiggliest, and `basis_params` holding `empirical_variance`
and `constraint_rank`.

## The construction

Three steps. The constraint matrix \\C = (\mathbf{1}, x)^\top B\\ is
formed and the basis restricted to its null space \\V_0\\, which makes
the remaining functions empirically orthogonal to a constant and to
\\x\\. The pencil \\(V_0^\top (B^\top B/n) V_0,\\ V_0^\top P V_0)\\ is
then diagonalized, and the transform is \\T = V_0 A\\.

The resulting design matrix \\Z = B T\\ then has \\Z^\top Z\\ diagonal
and satisfies \\(\mathbf{1}, x)^\top Z = 0\\. Measured on a B-spline of
twelve functions at 200 random points: the largest off-diagonal entry of
\\Z^\top Z\\ is 8.3e-14 and the largest entry of \\(\mathbf{1}, x)^\top
Z\\ is 1.3e-14.

## What `scale` does to the penalty

\\T^\top P T\\ is **proportional to** the identity, and equal to it only
when `scale = FALSE`. With `scale = TRUE`, the default, \\T\\ is divided
by \\\sqrt{\sum_j \lambda_j}\\, so \\T^\top P T\\ is \\(\sum_j
\lambda_j)^{-1} I\\: on the example above, 481.6 times the identity,
with `tr(Z'Z/n)` exactly 1. At `scale = FALSE` the two swap, \\T^\top P
T\\ being the identity to 8.3e-14 and `tr(Z'Z/n)` being 0.00208.

The scaling is what puts the bases of different terms on a common
footing, so one smoothing parameter means the same thing across them. A
consumer that needs the penalty to be exactly \\I\\ passes
`scale = FALSE`.

`basis_params$empirical_variance` holds the eigenvalues, normalized to
sum to one when `scale = TRUE`, and they are exactly `diag(Z'Z/n)`: the
share of empirical variance each column carries, falling from 0.834 for
the smoothest to 4.4e-05 for the wiggliest on the example above.

## What this construction costs

It factorizes only a \\q \times K\\ matrix and a \\(K-q) \times (K-q)\\
one, never anything of the size of the sample. It tolerates a
rank-deficient \\B\\, which equally spaced knots produce whenever the
data leave a knot span empty, the matrix inverted being the penalty. And
the transform is kept, so prediction at new points is the parent's
evaluation multiplied by it, as for any other transformed basis.

That last property is what separating a linear from a nonlinear effect
needs: a reparametrization not satisfying \\(\mathbf{1}, x)^\top Z = 0\\
estimates the sum of the two correctly and the split between them with
bias.

## Errors

A missing value in `x`, a `penalty` that is not `K` by `K`, or a
`constraints` with the wrong number of columns each throw. A penalty
singular on the constrained space throws with the two remedies named:
there is then a direction neither penalized nor identified, and no
rotation can fix it.

## References

Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
smoothing. *Numerische Mathematik* **24**, 375-382.

## See also

[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
which does its first step alone;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
which supplies the default penalty;
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
for a rotation that diagonalizes one matrix.

## Examples

``` r
set.seed(1)
x <- sort(runif(200))
d <- dr_basis(bspline_basis(dimension = 12), x)

# Two columns went to the constraint, leaving ten.
c(d@dimension, d@basis_params$constraint_rank)
#> [1] 10  2

z <- basis_eval(d, x)

# Z'Z is diagonal, and Z is orthogonal to a constant and to x.
max(abs(crossprod(z)[upper.tri(crossprod(z))]))
#> [1] 9.114426e-14
max(abs(crossprod(cbind(1, x), z)))
#> [1] 1.953993e-14

# The columns run from smoothest to wiggliest, and the recorded shares of
# empirical variance are exactly the diagonal of Z'Z/n.
round(d@basis_params$empirical_variance, 6)
#>  [1] 0.833632 0.115301 0.030148 0.011133 0.004913 0.002493 0.001384 0.000875
#>  [9] 0.000077 0.000044
max(abs(d@basis_params$empirical_variance - diag(crossprod(z)) / length(x)))
#> [1] 6.661338e-16

# Scaled, the penalty is a multiple of the identity and the trace is one.
P <- basis_gram(bspline_basis(dimension = 12), order = 2)
round(unique(round(diag(crossprod(d@transform, P %*% d@transform)), 6)), 4)
#> [1] 481.5973
sum(diag(crossprod(z))) / length(x)
#> [1] 1

# Unscaled, the penalty is the identity instead.
du <- dr_basis(bspline_basis(dimension = 12), x, scale = FALSE)
max(abs(crossprod(du@transform, P %*% du@transform) - diag(du@dimension)))
#> [1] 2.975398e-14
```
