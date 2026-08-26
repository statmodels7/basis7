# Construct a Fourier Basis

Returns a basis of a constant function and \\(K-1)/2\\ sine-cosine pairs
of increasing frequency on \\\[\ell, u\]\\. It is the basis to reach for
when the function being modeled is periodic or nearly so: over a whole
period the functions are mutually orthogonal, and stay orthogonal after
differentiation, so the Gram matrix is diagonal at every order and a
roughness penalty is one vector of numbers.

## Usage

``` r
fourier_basis(lower = 0, upper = 1, dimension = 5, omega = NULL)
```

## Arguments

- lower, upper:

  The endpoints of the interval, each a single finite number with
  `lower < upper`. Default \\\[0, 1\]\\. Evaluating outside throws.

- dimension:

  The number of basis functions, a single **odd** whole number of at
  least 1, default `5`. `1` is the constant alone, `5` is the constant
  and two pairs. An even value throws.

- omega:

  The period, a single positive finite number. `NULL`, the default, uses
  `upper - lower`, the value that makes the basis orthogonal on its
  interval. Anything not a single positive number throws.

## Value

An object of class
[FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md),
with `basis_name` `"fourier"`, `basis_params` holding `omega`, `n_pairs`
and `full_period`, and column names `const`, `sin1`, `cos1`, and so on.

## What the columns are

With \\z = 2\pi (x - \ell)/\omega\\, the columns are \\1\\, \\\sin z\\,
\\\cos z\\, \\\sin 2z\\, \\\cos 2z\\, and so on up to `n_pairs`
frequencies, named `const`, `sin1`, `cos1`, `sin2`, `cos2`.

## Why the dimension must be odd

A sine without its cosine represents a wave at one phase but not at the
next, so the fitted function would depend on where the interval was cut.
An even dimension throws, naming the two odd numbers on either side. It
is not adjusted silently: growing it would return a basis of a size the
caller did not ask for, and the constructor is the only place the
mismatch can be caught.

## The period, and what changes when it is not the interval

`omega` defaults to `upper - lower`, so the interval is exactly one
period and the functions are orthogonal on it. The order-\\d\\ Gram
matrix is then diagonal with entries \\\omega\\ for the constant at
order 0, zero for it above, and \\(\omega/2)(2\pi j/\omega)^{2d}\\ for
both members of pair \\j\\, written in closed form.

Any other positive period is accepted and every generic still answers,
but the interval is no longer a whole number of periods, the
orthogonality fails, and the Gram matrix is computed by composite
Gauss-Legendre instead. `basis_params$full_period` records which case
the object is in, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
reports `basis_gram` as `TRUE` there: the family says so itself through
[`basis_numerical_route.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.FourierBasis.md),
where reading which class the method is registered on would answer
`FALSE`, the owner being `FourierBasis` either way.

## Periodic by construction

Every column except the constant takes the same value at both ends of a
full period, so a fitted curve joins up. That is the property to want
here, and the reason not to reach for
[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
whose ends are free.

## See also

[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
for a local basis with free ends and
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
for a global polynomial one;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the roughness penalty this basis diagonalizes.

## Examples

``` r
b <- fourier_basis(dimension = 5)
b
#> Basis: fourier
#> Functions: 5   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   omega        1
#>   n_pairs      2
#>   full_period  TRUE
#> Numerical: none
basis_colnames(b)
#> [1] "const" "sin1"  "cos1"  "sin2"  "cos2" 

# Orthogonal on a whole period, so the Gram matrix is diagonal: omega for
# the constant, omega/2 for every sine and cosine.
round(basis_gram(b), 10)
#>       const sin1 cos1 sin2 cos2
#> const     1  0.0  0.0  0.0  0.0
#> sin1      0  0.5  0.0  0.0  0.0
#> cos1      0  0.0  0.5  0.0  0.0
#> sin2      0  0.0  0.0  0.5  0.0
#> cos2      0  0.0  0.0  0.0  0.5

# It stays diagonal after differentiating, with entries (2 pi j)^(2d) / 2.
round(diag(basis_gram(b, order = 2)), 4)
#>      const       sin1       cos1       sin2       cos2 
#>     0.0000   779.2727   779.2727 12468.3637 12468.3637 
c(0, rep((2 * pi * (1:2))^4 / 2, each = 2))
#> [1]     0.0000   779.2727   779.2727 12468.3637 12468.3637

# Periodic: the ends of a full period agree.
max(abs(basis_eval(b, 0) - basis_eval(b, 1)))
#> [1] 4.898587e-16

# An even dimension would leave half a pair, and is rejected.
try(fourier_basis(dimension = 4))
#> Error : 'dimension' must be odd: a Fourier basis holds a constant plus complete sine-cosine pairs, so 4 would leave half a pair. Use 3 or 5.

# A period that is not the interval width gives up the orthogonality, and
# the Gram matrix is then a quadrature and no longer diagonal.
f <- fourier_basis(dimension = 5, omega = 0.7)
f@basis_params$full_period
#> [1] FALSE
round(basis_gram(f), 4)
#>         const    sin1   cos1   sin2    cos2
#> const  1.0000  0.2118 0.0483 0.0210 -0.0436
#> sin1   0.2118  0.5218 0.0105 0.0061 -0.0832
#> cos1   0.0483  0.0105 0.4782 0.1286  0.0423
#> sin2   0.0210  0.0061 0.1286 0.5136  0.0170
#> cos2  -0.0436 -0.0832 0.0423 0.0170  0.4864
```
