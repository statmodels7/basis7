# Construct a Fourier Basis

A basis of a constant function and \\(K-1)/2\\ sine-cosine pairs of
increasing frequency on \\\[\ell, u\]\\.

## Usage

``` r
fourier_basis(lower = 0, upper = 1, dimension = 5, omega = NULL)
```

## Arguments

- lower, upper:

  The endpoints of the interval.

- dimension:

  The number of basis functions, an odd positive integer.

- omega:

  The period. Defaults to `upper - lower`.

## Value

An object of class
[`FourierBasis`](https://statmodels7.github.io/basis7/reference/FourierBasis.md).

## Details

The dimension must be odd. A sine without its cosine is a basis that can
represent a wave at one phase and not at another, so a dimension that
would leave a half pair is refused rather than adjusted: growing it
silently would return a basis of a size the caller did not ask for, and
the constructor is the only place the inconsistency can be caught.

The period `omega` defaults to the width of the interval, which is what
makes the basis functions orthogonal on it. A different period is
accepted and everything still works, but the interval is then no longer
a whole number of periods, so the Gram matrix stops being diagonal and
is computed by quadrature instead of in closed form.

## See also

[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
[`check_basis`](https://statmodels7.github.io/basis7/reference/check_basis.md)

## Examples

``` r
b <- fourier_basis(dimension = 5)
b
#> Basis: fourier
#> Functions: 5   Interval: [0, 1]
#> Parameters:
#>   omega        1
#>   n_pairs      2
#>   full_period  TRUE
#> Numerical: none
basis_colnames(b)
#> [1] "const" "sin1"  "cos1"  "sin2"  "cos2" 

# orthogonal on a whole period, so the Gram matrix is diagonal
round(basis_gram(b), 10)
#>       const sin1 cos1 sin2 cos2
#> const     1  0.0  0.0  0.0  0.0
#> sin1      0  0.5  0.0  0.0  0.0
#> cos1      0  0.0  0.5  0.0  0.0
#> sin2      0  0.0  0.0  0.5  0.0
#> cos2      0  0.0  0.0  0.0  0.5

# an even dimension would leave half a pair, and is refused
try(fourier_basis(dimension = 4))
#> Error : 'dimension' must be odd: a Fourier basis holds a constant plus complete sine-cosine pairs, so 4 would leave half a pair. Use 3 or 5.
```
