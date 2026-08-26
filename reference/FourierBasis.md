# Fourier Basis

The S7 class of Fourier bases, the objects
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
returns. It adds no property to
[basis](https://statmodels7.github.io/basis7/reference/basis.md) and
exists as the class the trigonometric methods dispatch on. A Fourier
basis holds a constant function and pairs of sines and cosines of
increasing frequency, and answers every derivative and its integral from
one identity.

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

  A single string naming the family, printed by
  [`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  and used by wrappers to build their own name. Not read by any
  computation.

- dimension:

  The number of functions in the basis, a single integer of at least 1.
  Must be of storage mode integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several. Both finite, of the same length, and
  `lower[j] < upper[j]` for every `j`. Their length is what
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  reports.

- basis_params:

  A named list of whatever else the subclass needs: the knots and degree
  of a B-spline, the frequency of a Fourier basis, the marginal
  dimensions of a product.
  [`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  shows it, abbreviating any numeric entry of more than four values.
  Defaults to an empty list.

## Value

An object of class `FourierBasis`, inheriting from
[basis](https://statmodels7.github.io/basis7/reference/basis.md), with
the same five properties and `basis_params` holding `omega`, `n_pairs`
and `full_period`. Call
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
instead of this class directly; it rejects an even dimension, defaults
the period and records all three.

## One identity for every order

Writing \\z = 2\pi (x - \ell)/\omega\\ for the phase and \\\omega\\ for
the period, \$\$\frac{\mathrm{d}^{k}}{\mathrm{d}x^{k}} \sin(j z) =
\left(\frac{2\pi j}{\omega}\right)^{k} \sin\\\left(j z +
\frac{k\pi}{2}\right),\$\$ and the same for the cosine. Differentiating
a sinusoid shifts its phase by a quarter turn and multiplies it by its
frequency, so no order is a special case and the fourth derivative costs
exactly what the first does.

The identity holds for negative \\k\\, which is where the antiderivative
comes from:
[`basis_int.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.FourierBasis.md)
evaluates it at \\k = -1\\ and subtracts the value at the lower
endpoint.

## Its `basis_params`

Three entries. `omega` is the period, `n_pairs` is
`(dimension - 1) %/% 2`, and `full_period` records whether `omega`
equals the width of the interval. The last decides which route
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
takes: a closed diagonal matrix when it is `TRUE`, a quadrature when it
is not.

## See also

[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
the constructor;
[`basis_eval.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_eval.FourierBasis.md),
[`basis_deriv.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_deriv.FourierBasis.md),
[`basis_int.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_int.FourierBasis.md)
and
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
for the methods registered on it.

## Examples

``` r
f <- fourier_basis(dimension = 5)
S7::S7_inherits(f, FourierBasis)
#> [1] TRUE
f@basis_params
#> $omega
#> [1] 1
#> 
#> $n_pairs
#> [1] 2
#> 
#> $full_period
#> [1] TRUE
#> 

# Differentiating shifts the phase by a quarter turn: the first derivative
# of sin1 at the start of the period is its own frequency, 2*pi.
basis_deriv(f, 0, order = 1)
#>      const     sin1         cos1     sin2         cos2
#> [1,]     0 6.283185 3.847341e-16 12.56637 7.694683e-16
2 * pi
#> [1] 6.283185
```
