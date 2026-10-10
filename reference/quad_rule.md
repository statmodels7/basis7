# Map a Quadrature Rule onto Intervals

Places an `n`-point Gauss-Legendre rule on each interval between
consecutive breakpoints and returns the pooled nodes and weights, so
that `sum(weights * f(nodes))` is the integral over the whole span. The
composite rule every quadrature in the package is built from.

## Usage

``` r
quad_rule(breaks, n)
```

## Arguments

- breaks:

  A numeric vector of at least two increasing breakpoints. The first and
  last are the ends of the span.

- n:

  The number of nodes per interval, a single positive whole number.

## Value

A list of two numeric vectors of length `n * (length(breaks) - 1)`:
`nodes` and `weights`, in the order described above.

## Details

Each interval gets the same `n` nodes, affinely mapped from \\\[-1,
1\]\\, and its weights scaled by half its width. The result is exact for
any function that is a polynomial of degree at most \\2n - 1\\ **on each
interval separately**, which is why the callers choose their breaks with
care:
[`basis_gram.BsplineBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.BsplineBasis.md)
uses the knots, so no interval straddles the point where a spline's
derivative jumps.

The function does not validate its arguments. `breaks` must be
increasing and of length at least two. The nodes are ordered by their
position within the rule and then by interval: the first node of every
interval comes first, then the second node of every interval, and so on.

## See also

[`gauss_legendre()`](https://statmodels7.github.io/basis7/reference/gauss_legendre.md),
which supplies the rule on \\\[-1, 1\]\\;
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
and
[`basis_int.basis()`](https://statmodels7.github.io/basis7/reference/basis_int.basis.md),
which consume it.
