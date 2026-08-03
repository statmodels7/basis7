# Numerical Gram Matrix of a Basis

The default inner-product method: composite Gauss-Legendre over the
basis interval, applied to the outer product of the requested
derivatives.

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order.

- panels:

  The number of subintervals.

- nodes:

  The number of quadrature nodes per subinterval.

- ...:

  Unused.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

Numerical Gram Matrix of a Basis

The rule is placed on `panels` equal subintervals. A basis with
closed-form inner products, or one that is piecewise polynomial and so
can be integrated exactly by a rule sized from its degree, registers its
own method and this one is not used.
