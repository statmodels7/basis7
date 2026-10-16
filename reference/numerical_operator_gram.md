# The Roughness Matrix of an Operator by Quadrature

Integrates \\(Lb)(Lb)^\top\\ numerically: over Gauss-Legendre panels for
the length measure, over the covariate values for the empirical measure,
and against a density where one is given. It is the base method of
[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
and the route that a family with no closed form takes.

## Usage

``` r
numerical_operator_gram(
  basis,
  op,
  at = NULL,
  weight = NULL,
  panels = 50L,
  nodes = 12L,
  ...
)
```

## Arguments

- basis:

  A [basis](https://statmodels7.github.io/basis7/reference/basis.md) of
  one variable.

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md)
  whose period is resolved.

- at:

  The covariate values for the empirical measure, or `NULL`.

- weight:

  A density with respect to which the integral is taken, or `NULL`.

- panels, nodes:

  The quadrature rule: the number of equal panels into which the
  interval is split and the number of Gauss-Legendre nodes that each
  panel carries.

- ...:

  Ignored.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns.

## Details

The nodes are those of the rule used by
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md),
and the two functions differ only in what they evaluate at the nodes:
this function evaluates \\Lb\\ through
[`operator_eval()`](https://statmodels7.github.io/basis7/reference/operator_eval.md)
instead of one derivative. For a piecewise polynomial basis this
function therefore returns an approximation, with the accuracy of the
quadrature, where the derivative Gram matrix of the basis itself is
exact.

## See also

[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md),
the generic.
