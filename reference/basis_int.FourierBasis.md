# Integral of a Fourier Basis

The exact integral from the lower endpoint.

## Arguments

- basis:

  A
  [`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points.

- ...:

  Unused.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns.

## Details

The shift identity at order \\-1\\ gives an antiderivative, not the one
this package's convention asks for, and the two differ by a constant
that is not the same for every column: at the lower endpoint the sine
columns of the raw antiderivative are \\-\omega/(2\pi j)\\ while the
cosine columns are already zero. Subtracting the value at the lower
endpoint fixes every column at once and makes the convention hold
exactly.
