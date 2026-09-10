# The Roughness Matrix of a Differential Operator

The Gram matrix of \\Lb\\, that is \$\$R = \int_a^b (Lb)(Lb)^\top \\
\mathrm{d}\mu,\$\$ the matrix for which \\\lVert Lx \rVert^2 = c^\top R
c\\ when \\x = b^\top c\\. It is what
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
returns when its `order` is a
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
and the generic exists so that a family with a closed form can declare
one.

## Usage

``` r
basis_operator_gram(basis, op, at = NULL, weight = NULL, ...)
```

## Arguments

- basis:

  A [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- at:

  The covariate values, for the empirical measure, or `NULL`.

- weight:

  A density to integrate against, or `NULL`.

- ...:

  Passed to methods, and on to the quadrature (`panels`, `nodes`).

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns.

## Details

The base method integrates numerically, evaluating \\Lb\\ through
[`operator_eval()`](https://statmodels7.github.io/basis7/reference/operator_eval.md)
at Gauss-Legendre nodes, at the covariate values for the empirical
measure, or against a weight function. A
[FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
over a full period overrides it with an exact diagonal form; see
[`basis_operator_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.FourierBasis.md).

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
the generic that routes to it;
[`operator_eval()`](https://statmodels7.github.io/basis7/reference/operator_eval.md)
for \\Lb\\.

## Examples

``` r
b <- fourier_basis(lower = 0, upper = 365, dimension = 7)
round(diag(basis_gram(b, order = harmonic_operator(365))), 6)
#> const  sin1  cos1  sin2  cos2  sin3  cos3 
#> 0e+00 0e+00 0e+00 0e+00 0e+00 3e-06 3e-06 
```
