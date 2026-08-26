# Gram Matrix Against the Empirical Measure

Computes \\B^{(d)\top} B^{(d)} / n\\ at the given points: the inner
products a design matrix produces, in place of those of the functions on
their interval. Called from the body of
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
when `at` is supplied, so no method ever sees this case.

## Usage

``` r
empirical_gram(basis, order, at)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- order:

  The derivative order, already validated by
  [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md).

- at:

  A numeric vector of points, or a matrix of
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  columns. Not range-checked here;
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  does that and throws for a point outside the interval.

## Value

A symmetric numeric matrix of `basis@dimension` rows and columns, with
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## Details

Missing points are dropped **before** the basis is evaluated. A basis is
entitled to refuse a vector that is entirely missing, and its refusal
would name the wrong thing here, so `at` with no usable point throws
`'at' has no usable points.` instead. For a basis of several variables
`at` is coerced to a matrix and rows with any missing entry are dropped
whole.

The result is symmetrized as `(G + t(G))/2` before it is returned, the
two triangles of a crossproduct differing in their last bits, and given
[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
on both margins.

## See also

[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
its only caller, and
[`weighted_gram()`](https://statmodels7.github.io/basis7/reference/weighted_gram.md)
for the other alternative measure.
