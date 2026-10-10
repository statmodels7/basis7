# Read an Order Argument as an Operator

Normalizes the `order` argument every smoother family takes: a whole
number `m` becomes `deriv_operator(m)` and an operator is returned
unchanged, so that `order = 2` is a shorthand for `deriv_operator(2)`.

## Usage

``` r
as_operator(order, nm = "order")
```

## Arguments

- order:

  A whole number or a
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

- nm:

  The argument's name, for the error message.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).
