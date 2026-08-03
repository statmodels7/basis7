# Gram Matrix Against a Weighted Lebesgue Measure

\\\int B^{(d)} B^{(d)\top} w(t)\\ \mathrm{d}t\\, by composite
Gauss-Legendre. A weight is an arbitrary function, so no family has a
closed form for it and the quadrature is always used.

## Usage

``` r
weighted_gram(basis, order, weight, panels = 50L, nodes = 12L, ...)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order.

- weight:

  A function of one numeric vector.

- panels:

  The number of subintervals.

- nodes:

  The number of quadrature nodes per subinterval.

- ...:

  Unused.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.
