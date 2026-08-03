# The Trigonometric Columns of a Fourier Basis

The sine and cosine columns at derivative order `d`, from the shift
identity. Order \\-1\\ gives an antiderivative.

## Usage

``` r
fourier_trig(basis, x, d)
```

## Arguments

- basis:

  A
  [`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- d:

  The order, which may be negative.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension - 1`
columns.
