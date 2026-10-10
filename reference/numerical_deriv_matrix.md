# Numerically Differentiate a Matrix-Valued Function

Estimates the `order`-th derivative of `f` at each point of `x` by a
single finite-difference stencil, symmetric where the interval leaves
room for it and one-sided where it does not. One stencil of the order
wanted, never a composition of lower-order differences, so the error is
the truncation of that one formula. It computes the numerical
derivatives of every basis; the numerical integral and Gram matrix use
Gauss-Legendre quadrature through
[`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md)
instead.

## Usage

``` r
numerical_deriv_matrix(f, x, order, lower, upper, step_scale = 1)
```

## Arguments

- f:

  A function of one numeric vector returning a numeric matrix with one
  row per element. Called `max(2 * reach + 1, order + 2)` times, from
  three times at order 1 to six at order 4.

- x:

  A numeric vector of evaluation points. `NA` entries are given the
  central stencil and propagate to an `NA` row.

- order:

  The derivative order, a single positive whole number.

- lower, upper:

  The endpoints of the interval on which `f` is defined, used to choose
  each point's stencil and to cap the step.

- step_scale:

  A multiplier on the step, default `1`.
  [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  passes `0.5` to measure the uncertainty of the reference by the change
  in the result.

## Value

A numeric matrix with `length(x)` rows and as many columns as `f`
returns, with no dimnames; callers add them through
[`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md).

## Choice of the stencil

A basis is evaluated at its endpoints as readily as anywhere else, and a
symmetric stencil centered on an endpoint would need points outside the
interval, where the basis signals an error. Each point therefore gets
the central stencil when both sides have room, and otherwise the
one-sided stencil that points inward. The one-sided stencils have
`order + 2` nodes, so that they keep the second-order accuracy of the
central one, which at an even order has one node fewer; the central
stencil is padded with zero weights to the same length. The weights come
from
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html),
and the central offsets from
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html).

The one-sided stencils carry a larger error constant. Measured against
exact Legendre derivatives of dimension 5 on \\\[0, 1\]\\, relative to
the scale of the result, the first derivative agrees to 5.1e-10 at the
two endpoints and to 3.9e-10 in the interior, the second to 1.3e-07 and
1.5e-08.

## The step

The step is that of
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
\\\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)\\, which balances
truncation against rounding for order \\d\\. The rule is not repeated in
this package, so that it has a single definition in numericals7.

It is then capped at `0.4 * (upper - lower) / (2 * reach)`, `reach`
being the half-width of the central stencil, so that every stencil fits
inside the interval: `0.2` of the width at orders 1 and 2, and `0.1` at
orders 3 and 4.

## What it costs in accuracy

Measured against exact Legendre derivatives of dimension 5 at 21
interior points of \\\[0, 1\]\\, relative to the scale of the result:
3.9e-10 at order 1, 1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08
at order 4, and about ten times more at the endpoints from order 2 on. A
closed form, where one exists, is exact.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html),
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
which supply the weights, the central offsets and the step;
[`basis_deriv.basis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.basis.md),
its main caller.
