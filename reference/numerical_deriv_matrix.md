# Numerically Differentiate a Matrix-Valued Function

Estimates the `order`-th derivative of `f` at each point of `x` by a
single finite-difference stencil, symmetric where the interval leaves
room for it and one-sided where it does not. One stencil of the order
wanted, never a composition of lower-order differences, so the error is
the truncation of that one formula. The engine behind every numerical
fallback in the package.

## Usage

``` r
numerical_deriv_matrix(f, x, order, lower, upper, step_scale = 1)
```

## Arguments

- f:

  A function of one numeric vector returning a numeric matrix with one
  row per element. Called `2 * reach + 1` times, so a costly `f` is
  evaluated three or five times over.

- x:

  A numeric vector of evaluation points. `NA` entries are given the
  central stencil and propagate to an `NA` row.

- order:

  The derivative order, a single positive whole number.

- lower, upper:

  The endpoints of the interval `f` is defined on, used to choose each
  point's stencil and to cap the step.

- step_scale:

  A multiplier on the step, default `1`.
  [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  passes `0.5` to measure the reference's own uncertainty by how much
  the answer moves.

## Value

A numeric matrix with `length(x)` rows and as many columns as `f`
returns, with no dimnames; callers add them through
[`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md).

## One stencil per point, chosen by the room available

A basis is evaluated at its endpoints as readily as anywhere else, and a
symmetric stencil centered on an endpoint would ask for points outside
the interval, where the basis throws. Each point therefore gets the
central stencil when both sides have room, and otherwise the one-sided
stencil that points inward, with the same order and the same number of
nodes. All three weight vectors come from
[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html)
on the offsets
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
supplies.

Accuracy is not lost at the ends. Measured against exact Legendre
derivatives on \\\[0, 1\]\\, the first derivative agrees to 5.1e-10
relative at the two endpoints and to 3.9e-10 in the interior.

## The step

The step is
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)'s,
\\\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)\\, which balances
truncation against rounding for order \\d\\. It is written nowhere in
this package, a second copy of a rule with one home being how two
packages come to disagree.

It is then capped at `0.4 * (upper - lower) / (2 * reach)` so that the
whole stencil fits inside the interval: `0.2` of the width at orders 1
and 2, where the stencil has three nodes, and `0.1` at orders 3 and 4,
where it has five.

## What it costs in accuracy

Measured against exact Legendre derivatives at 21 interior points of
\\\[0, 1\]\\, relative to the scale of the answer: 3.9e-10 at order 1,
1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08 at order 4. That is
what a family which registers no derivative method gets, and the reason
to write a closed form where one exists.

## See also

[`numericals7::fd_weights()`](https://statmodels7.github.io/numericals7/reference/fd_weights.html),
[`numericals7::fd_offsets()`](https://statmodels7.github.io/numericals7/reference/fd_offsets.html)
and
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html),
which supply the three pieces;
[`basis_deriv.basis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.basis.md),
its main caller.
