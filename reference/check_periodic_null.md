# Refuse an Operator Whose Null Space a Periodic Basis Cannot Carry

Every function of the operator's null space other than the constant is
restored as a free column of a periodic block, so each must itself be
periodic on the basis's period: a pure sine or cosine at a whole
multiple of the fundamental frequency, with no power of \\t\\ and no
exponential in front of it.

## Usage

``` r
check_periodic_null(sm, op, x)
```

## Arguments

- sm:

  A
  [FourierSmoother](https://statmodels7.github.io/basis7/reference/FourierSmoother.md).

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- x:

  The covariate, from which the interval is taken where the smoother
  does not fix it.

## Value

`NULL`, invisibly; called for the error.

## Details

Restoring anything else gives a block whose columns are not all
periodic, and a fit on it no longer takes the same value at the two ends
of the interval, which is the one property a Fourier basis is chosen
for. The check is on the operator's own null space, read analytically,
rather than on the rank of the assembled penalty: the two are different
questions and the second gives an answer that moves with the tolerance.
