# Gauss-Legendre Nodes and Weights

The `n`-point Gauss-Legendre rule on \\\[-1, 1\]\\, which integrates
polynomials of degree up to \\2n - 1\\ exactly.

## Usage

``` r
gauss_legendre(n)
```

## Arguments

- n:

  The number of nodes, a positive integer.

## Value

A list with components `nodes` and `weights`.

## Details

The nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix
of the Legendre recurrence and the weights come from the first component
of each eigenvector, which is the Golub-Welsch construction. Computing
them rather than tabulating them keeps any node count available, which
the exact spline rules need: a rule per knot interval, sized from the
degree.
