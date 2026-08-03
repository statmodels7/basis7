# Finite-Difference Weights for an Arbitrary Stencil

The weights that turn function values at `x + offsets * h` into the
`order`-th derivative at `x`, divided by `h^order`.

## Usage

``` r
fd_weights(offsets, order)
```

## Arguments

- offsets:

  A numeric vector of stencil offsets, in units of the step.

- order:

  The derivative order.

## Value

A numeric vector of weights, the same length as `offsets`.

## Details

The weights solve the linear system that makes the combination exact on
polynomials up to the degree the stencil can carry. Writing \\s_j\\ for
the offsets, the Taylor expansion of \\\sum_j w_j f(x + s_j h)\\ has
\\f^{(i)}(x)\\ multiplied by \\h^i/i! \sum_j w_j s_j^i\\, so requiring
\\\sum_j w_j s_j^i = 0\\ for \\i \neq d\\ and \\= d!\\ for \\i = d\\
leaves exactly \\h^d f^{(d)}(x)\\. That is a Vandermonde system in the
offsets, solved once per stencil shape.

Building the weights this way, rather than composing lower-order
differences, is what keeps a high order usable: each numerical
differentiation multiplies the error of the one before it, so a fourth
derivative reached by four nested first differences is noise. One
stencil, never nested.
