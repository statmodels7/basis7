# Gauss-Legendre Nodes and Weights

Returns the `n`-point Gauss-Legendre rule on \\\[-1, 1\]\\: the nodes
and the weights that integrate every polynomial of degree up to \\2n -
1\\ exactly. At `n = 5` the rule reproduces \\\int t^9\\ as `0` and
first departs at \\t^{10}\\, returning 0.17889 against 0.18182.

## Usage

``` r
gauss_legendre(n)
```

## Arguments

- n:

  The number of nodes, a single positive whole number. `1` returns the
  midpoint rule directly. `n < 1` throws
  `'n' must be a positive integer.`

## Value

A list of two numeric vectors of length `n`: `nodes`, in increasing
order, and `weights`, positive and summing to 2.

## Details

The construction is Golub-Welsch: the nodes are the eigenvalues of the
symmetric tridiagonal Jacobi matrix of the Legendre recurrence, with
off-diagonal \\i/\sqrt{4i^2 - 1}\\, and the weights are twice the square
of the first component of each eigenvector. The weights sum to 2, the
length of the interval.

They are computed at call time, which keeps every node count available.
That is what the exact spline rules need: one rule per knot interval,
sized from the degree and the derivative order, where a table would
offer only the counts someone thought to tabulate.

## References

Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
rules. *Mathematics of Computation* **23**, 221-230.

## See also

[`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md),
which maps this rule onto a sequence of intervals.
