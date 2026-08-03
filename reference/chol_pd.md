# Cholesky Factorisation, With the Rank Decided Before It

The Cholesky factor of a symmetric matrix, or `NULL` when the matrix is
not positive definite to the given relative tolerance.

## Usage

``` r
chol_pd(m, tol = 1e-12)
```

## Arguments

- m:

  A symmetric numeric matrix.

- tol:

  The relative tolerance below which the smallest eigenvalue counts as
  zero.

## Value

The upper triangular Cholesky factor, or `NULL`.

## Details

The verdict comes from the eigenvalues rather than from whether
[`chol`](https://rdrr.io/r/base/chol.html) raises. On a matrix with an
exactly zero eigenvalue the pivot that should be zero comes out positive
or negative according to rounding, so
[`chol()`](https://rdrr.io/r/base/chol.html) succeeds on some platforms
and fails on others, and a construction that asks it whether a penalty
is usable gets a different answer on different machines. Comparing the
smallest eigenvalue with the largest is a statement about the matrix,
and gives the same answer everywhere.
