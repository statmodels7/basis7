# Numerically Differentiate a Matrix-Valued Function

The `order`-th derivative of `f` at each point of `x`, by a single
finite-difference stencil, symmetric where the interval leaves room for
it and one-sided at the endpoints.

## Usage

``` r
numerical_deriv_matrix(f, x, order, lower, upper, step_scale = 1)
```

## Arguments

- f:

  A function of a numeric vector returning a matrix with one row per
  element.

- x:

  A numeric vector of evaluation points.

- order:

  The derivative order.

- lower, upper:

  The endpoints of the interval `f` is defined on.

- step_scale:

  A factor applied to the step. Halving it is how
  [`fd_reference`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  measures its own uncertainty.

## Value

A numeric matrix with `length(x)` rows.

## Details

A basis is evaluated at its endpoints as readily as anywhere else, and a
symmetric stencil centred on an endpoint would ask for points outside
the interval, where the basis is not defined. Those points therefore get
a one-sided stencil of the same order and the same number of nodes,
built by
[`fd_weights`](https://statmodels7.github.io/basis7/reference/fd_weights.md)
from shifted offsets.

The step is \\\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)\\, which
balances truncation against rounding for order \\d\\, capped so that the
whole stencil fits inside the interval.
