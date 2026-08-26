# Validate a Derivative Order

Checks that `order` is a non-negative whole number, or a vector of them
of length `nvar`, and returns it as an integer vector of length `nvar`.
Called from the bodies of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
so both report the same errors in the same words.

## Usage

``` r
check_order(order, nvar = 1L)
```

## Arguments

- order:

  The value supplied by the caller: a single non-negative whole number,
  or a vector of `nvar` of them.

- nvar:

  The number of variables the basis takes, from
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md).
  Default `1`.

## Value

`order` as an integer vector of length `nvar`.

## Details

Anything failing the first test throws
`'order' must be a non-negative integer.`; that covers a negative value,
a fraction, an infinity, an `NA` and a non-numeric.

A vector of length `nvar` is returned as it stands, and a single `0` is
repeated to that length. For a basis of several variables a single
**non-zero** order throws with a longer message, because a scalar has
two readings there: that order in every coordinate, or that total order.
Choosing one silently would fit a different model from the one the
caller wrote. Zero is exempt, meaning no derivative under either
reading.

## See also

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
its two callers.
