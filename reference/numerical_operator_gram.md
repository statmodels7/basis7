# The Roughness Matrix of an Operator by Quadrature

Integrates \\(Lb)(Lb)^\top\\ numerically: over Gauss-Legendre panels for
the length measure, over the covariate values for the empirical measure,
and against a density where one is given. It is the base method of
[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
and the route a family with no closed form takes.

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
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- at:

  The covariate values for the empirical measure, or `NULL`.

- weight:

  A density to integrate against, or `NULL`.

- panels, nodes:

  The quadrature: how many equal panels the interval is split into and
  how many Gauss-Legendre nodes each carries.

- ...:

  Ignored.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns.

## Details

The nodes are the same rule
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
uses, and the only difference is what is evaluated at them: \\Lb\\ from
[`operator_eval()`](https://statmodels7.github.io/basis7/reference/operator_eval.md)
rather than one derivative. A piecewise polynomial basis therefore gets
an approximation here where its own derivative Gram matrix is exact, and
the accuracy is the quadrature's.

## See also

[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md),
the generic.
