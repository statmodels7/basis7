# Gram Matrix of a B-Spline Basis

Exact inner products, integrated knot interval by knot interval.

## Arguments

- basis:

  A
  [`BsplineBasis`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  object.

- order:

  The derivative order.

- ...:

  Unused.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.

## Details

On each knot interval the order-\\d\\ derivative of a spline of degree
\\m\\ is a polynomial of degree \\m - d\\, so the integrand of the Gram
matrix has degree \\2(m - d)\\. A Gauss-Legendre rule with \\m - d + 1\\
nodes integrates degree \\2(m - d) + 1\\ exactly, so the result carries
no quadrature error, only floating-point error.
