# Read an Order Argument as an Operator

Normalizes the `order` argument every smoother family takes: a whole
number `m` becomes `deriv_operator(m)` and an operator is returned
unchanged. It is what makes `order = 2` the shorthand rather than a
second way of saying the same thing.

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
