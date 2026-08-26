# Cholesky Factorization, With the Rank Decided Before It

Returns the upper triangular Cholesky factor of a symmetric matrix, or
`NULL` when the matrix is not positive definite to the given relative
tolerance. The verdict comes from the eigenvalues and never from whether
[`base::chol()`](https://rdrr.io/r/base/chol.html) raised, so it is the
same on every platform.

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
`NULL` when `m` is empty, holds an `NA`, or is not positive definite to
`tol`.

## Details

On a matrix with an exactly zero eigenvalue the pivot that should be
zero comes out positive or negative according to rounding, so
[`chol()`](https://rdrr.io/r/base/chol.html) succeeds on some platforms
and fails on others. A construction that asks it whether a Gram matrix
or a penalty is usable therefore gets a different answer on different
machines:
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
once disagreed with themselves between this machine and the CI runners
for exactly that reason.

Comparing the smallest eigenvalue with the largest is a statement about
the matrix and gives the same answer everywhere. The eigendecomposition
costs little beside what both callers already do.

The [`chol()`](https://rdrr.io/r/base/chol.html) call is still wrapped,
so a matrix that passes the eigenvalue test and fails the factorization
anyway returns `NULL` instead of throwing.

## See also

[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
its two callers, which turn a `NULL` into an error naming what to do
about it.
