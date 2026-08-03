# Gram Matrix by Composite Quadrature

The inner products of the order-`order` derivatives, by Gauss-Legendre
on equal subintervals. Shared by the default method and by any basis
whose closed form does not apply to the arguments it was given.

## Usage

``` r
numerical_gram(basis, order = 0L, panels = 50L, nodes = 12L)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order.

- panels:

  The number of subintervals.

- nodes:

  The number of quadrature nodes per subinterval.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.
