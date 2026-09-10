# Apply an Operator to a Basis

Evaluates \\Lb(x)\\, the operator applied to every function of a basis:
the matrix whose Gram matrix is the roughness the penalty measures.

## Usage

``` r
operator_eval(basis, x, op)
```

## Arguments

- basis:

  A [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- x:

  The points to evaluate at.

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

## Value

A numeric matrix of `length(x)` rows and `basis@dimension` columns.

## Details

It is one weighted sum of derivative evaluations, \\\sum_j w_j D^j b +
D^m b\\, and every term of it comes from
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
so a basis that answers its derivatives answers this. Terms whose weight
is exactly zero are skipped, which is why the pure derivative operator
costs one call and not \\m + 1\\.
