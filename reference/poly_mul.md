# Multiply Two Polynomials Given by Their Coefficients

Convolves two coefficient vectors written in increasing degree. The
convolution multiplies the two polynomials, as composing two
differential operators multiplies their characteristic polynomials.

## Usage

``` r
poly_mul(a, b)
```

## Arguments

- a, b:

  Numeric vectors of coefficients in increasing degree.

## Value

A numeric vector of `length(a) + length(b) - 1` coefficients, in
increasing degree.
