# Check That a Periodic Basis Is Wider Than What the Operator Removes

Signals an error at construction when the constant and the functions of
the operator's null space take all `k` functions of the basis, leaving
none to penalize. The count is the order of the operator, plus one where
the constant is not in its null space, as for
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md).

## Usage

``` r
check_periodic_size(op, k)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  resolved or not.

- k:

  The number of basis functions.

## Value

`NULL`, invisibly; called for the error.
