# Construct a Legendre Polynomial Basis

Returns a basis of the polynomials of degree below `dimension` on
\\\[\ell, u\]\\, written in the Legendre polynomials and not in raw
powers. The two span the same space, so a fit through either gives the
same curve; the Legendre form is the one that survives being fitted, its
Gram matrix being diagonal where the raw powers' is a Hilbert matrix.

## Usage

``` r
poly_basis(lower = 0, upper = 1, dimension = 4)
```

## Arguments

- lower, upper:

  The endpoints of the interval, each a single finite number with
  `lower < upper`. Default \\\[0, 1\]\\. Evaluating outside throws.

- dimension:

  The number of polynomials, a single whole number of at least 1,
  default `4`. The highest degree is `dimension - 1`, so `1` gives the
  constant alone.

## Value

An object of class
[PolyBasis](https://statmodels7.github.io/basis7/reference/PolyBasis.md),
with `basis_name` `"legendre"`, `basis_params` holding
`degree = dimension - 1`, and column names `P0`, `P1`, and so on.

## Why not raw powers

Fitting \\1, x, x^2, \ldots\\ on \\\[0, 1\]\\ gives the Gram matrix
\\H\_{mn} = 1/(m + n + 1)\\, the Hilbert matrix, whose condition number
grows geometrically. The condition numbers of the two order-0 Gram
matrices, measured on \\\[0, 1\]\\:

|             |          |            |
|-------------|----------|------------|
| `dimension` | Legendre | raw powers |
| 5           | 9        | 4.8e+05    |
| 10          | 19       | 1.6e+13    |
| 15          | 29       | 2.5e+17    |

At ten raw powers a least-squares solve has already lost most of its
digits; at fifteen the matrix is numerically singular. The Legendre
condition number is \\2K - 1\\, growing linearly.

## The Gram matrix

The polynomials are orthogonal on their own interval, \\\int\_{-1}^{1}
P_m P_n \\\mathrm{d}t = 2\delta\_{mn}/(2n+1)\\, so `basis_gram(p)` is
diagonal with entries \\(u - \ell)/(2n+1)\\ for \\n = 0, \ldots, K-1\\.
They are **not** orthonormal: pass the basis through
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
for that.

## Degree and dimension

`dimension` counts functions and the degrees start at zero, so
`poly_basis(dimension = 4)` holds \\P_0, P_1, P_2, P_3\\, a cubic. The
basis contains the constant, so it is collinear with an intercept in the
same design; drop one or the other, or use
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md).

## References

Abramowitz, M. and Stegun, I. A. (1964). *Handbook of Mathematical
Functions*, chapter 22. National Bureau of Standards.

## See also

[`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
for a local basis and
[`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
for a periodic one;
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
to make the Gram matrix the identity;
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
for the roughness penalty of this basis.

## Examples

``` r
p <- poly_basis(dimension = 5)
p
#> Basis: legendre
#> Functions: 5   Variables: 1
#> Domain: [0, 1]
#> Parameters:
#>   degree  4
#> Numerical: none

# Orthogonal, so the Gram matrix is diagonal, with entries (u - l)/(2n + 1).
round(diag(basis_gram(p)), 8)
#>        P0        P1        P2        P3        P4 
#> 1.0000000 0.3333333 0.2000000 0.1428571 0.1111111 
1 / (2 * (0:4) + 1)
#> [1] 1.0000000 0.3333333 0.2000000 0.1428571 0.1111111

# And well conditioned where the raw powers are not.
H <- outer(0:9, 0:9, function(a, b) 1 / (a + b + 1))
c(legendre = kappa(basis_gram(poly_basis(dimension = 10)), exact = TRUE),
  raw_powers = kappa(H, exact = TRUE))
#>     legendre   raw_powers 
#> 1.900000e+01 1.602442e+13 

# It spans the same space as the raw powers, so a cubic is fitted exactly.
x <- seq(0, 1, length.out = 40)
max(abs(lm.fit(basis_eval(p, x), x^3)$fitted.values - x^3))
#> [1] 1.193503e-16

# P_n(1) = 1 and P_n(-1) = (-1)^n, at the two ends of the interval.
basis_eval(p, c(0, 1))
#>      P0 P1 P2 P3 P4
#> [1,]  1 -1  1 -1  1
#> [2,]  1  1  1  1  1

# Derivatives above the highest degree are exactly zero.
range(basis_deriv(p, c(0.2, 0.8), order = 5))
#> [1] 0 0
```
