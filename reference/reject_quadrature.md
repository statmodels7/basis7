# Reject Quadrature Settings on an Exact Route

Signals an error when `panels` or `nodes` reaches a method that computes
its Gram matrix exactly, where the two would otherwise be accepted and
have no effect.

## Usage

``` r
reject_quadrature(...)
```

## Arguments

- ...:

  The arguments a method received through its own `...`.

## Value

`NULL`, invisibly; called for the error.
