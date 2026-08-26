# The Trigonometric Columns of a Fourier Basis

Builds the sine and cosine columns at derivative order `d` from the
phase-shift identity \$\$\frac{\mathrm{d}^{d}}{\mathrm{d}x^{d}} \sin(jz)
= \left(\frac{2\pi j}{\omega}\right)^{d} \sin\\\left(jz +
\frac{d\pi}{2}\right),\$\$ the one routine behind
[`basis_eval.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.FourierBasis.md)
at `d = 0`,
[`basis_deriv.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.FourierBasis.md)
above it and
[`basis_int.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.FourierBasis.md)
at `d = -1`. The constant column is not included; each caller prepends
its own.

## Usage

``` r
fourier_trig(basis, x, d)
```

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points. Not range-checked here.

- d:

  The order, a single whole number that may be negative. `-1` gives an
  antiderivative, `0` the functions themselves.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension - 1`
columns, no dimnames, sine and cosine interleaved by frequency.

## Details

At `d = -1` the identity gives an antiderivative whose constant is not
the one
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
promises; the caller subtracts the row at the lower endpoint. Columns
are interleaved sine-then-cosine per frequency, matching
[`basis_colnames.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.FourierBasis.md).

At `n_pairs == 0` the result is a `length(x)` by 0 matrix, and the loop
is skipped.

## See also

[`basis_eval.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.FourierBasis.md),
[`basis_deriv.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.FourierBasis.md)
and
[`basis_int.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.FourierBasis.md),
its three callers.
