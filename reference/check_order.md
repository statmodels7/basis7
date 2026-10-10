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

Anything failing the first test signals the error
`'order' must be a non-negative integer.`; that covers a negative value,
a fraction, an infinity, an `NA` and a non-numeric.

A vector of length `nvar` is returned as it stands, and a single `0` is
repeated to that length. For a basis of several variables a single
**non-zero** order signals an error, because a scalar has two readings
there: that order in every coordinate, or that total order. Zero is
exempt, meaning no derivative under either reading. A vector of any
other length signals an error naming the length expected.

## See also

[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
and
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
its callers, together with
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
and the Gram method of a tensor basis.
