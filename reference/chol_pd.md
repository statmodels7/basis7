# Cholesky Factor with a Positive-Definiteness Test

Returns the upper triangular Cholesky factor of a symmetric matrix, or
`NULL` when the matrix is not positive definite to the given relative
tolerance. The test is made on the eigenvalues, so that its outcome is
the same on every platform.

## Usage

``` r
chol_pd(m, tol = 1e-12)
```

## Arguments

- m:

  A symmetric numeric matrix. Symmetry is assumed, never checked:
  `eigen(symmetric = TRUE)` reads the lower triangle.

- tol:

  The relative tolerance below which the smallest eigenvalue counts as
  zero, default `1e-12`. The test is `min(ev) > tol * max(ev)`, so it is
  scale-free.

## Value

The upper triangular Cholesky factor \\R\\ with \\m = R^\top R\\, or
`NULL` when `m` is empty, holds a value that is not finite, or is not
positive definite to `tol`.

## Details

On a matrix with an exactly zero eigenvalue the pivot that should be
zero comes out positive or negative according to rounding, so
[`base::chol()`](https://rdrr.io/r/base/chol.html) succeeds on some
platforms and fails on others. A construction that relied on
[`chol()`](https://rdrr.io/r/base/chol.html) to decide whether a Gram
matrix or a penalty is usable would therefore give different results on
different platforms.

Comparing the smallest eigenvalue with the largest depends only on the
matrix and gives the same result everywhere. The eigendecomposition
costs little beside the other work of the callers.

The [`chol()`](https://rdrr.io/r/base/chol.html) call is still wrapped,
so a matrix that passes the eigenvalue test and fails the factorization
returns `NULL` instead of signalling an error.

## See also

[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
and
[`smoother_reparam()`](https://statmodels7.github.io/basis7/reference/smoother_reparam.md),
its callers, which turn a `NULL` into an error stating what to change.
