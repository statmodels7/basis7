# Fourier Basis

The S7 class of Fourier bases. Constructed by
[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md).

## Usage

``` r
FourierBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list()
)
```

## Arguments

- basis_name:

  A short name for the family, used when printing.

- dimension:

  The number of functions in the basis, a positive integer.

- lower, upper:

  The endpoints of the interval the basis lives on.

- basis_params:

  A named list of whatever else the subclass needs.

## Value

An object of class `FourierBasis`. Use
[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
rather than calling the class directly, so that the dimension is checked
and the period recorded.

## Details

The basis holds a constant function and pairs of sines and cosines of
increasing frequency. Every derivative and the antiderivative come from
one identity, so no order is special:
\$\$\frac{\mathrm{d}^{k}}{\mathrm{d}x^{k}} \sin(j z) = \left(\frac{2\pi
j}{\omega}\right)^{k} \sin\\\left(j z + \frac{k\pi}{2}\right), \qquad z
= \frac{2\pi (x - \ell)}{\omega},\$\$ and likewise for the cosine. The
identity holds for negative \\k\\ as well, which is where the
antiderivative comes from.

## See also

[`fourier_basis`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)

## Examples

``` r
f <- fourier_basis(dimension = 5)
S7::S7_inherits(f, FourierBasis)
#> [1] TRUE
```
