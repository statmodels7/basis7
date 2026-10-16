# Is This a Differential Operator?

Reports whether an object is a
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).
The `order` argument of a smoother and of
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
uses this test to distinguish its two spellings, a whole number and an
operator.

## Usage

``` r
is_operator(x)
```

## Arguments

- x:

  Any object.

## Value

A single logical.
