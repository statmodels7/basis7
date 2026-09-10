# The Harmonic Acceleration Operator

The operator for periodic data, \$\$L x = D \prod\_{i=1}^{h} \left(D^2 +
(i\nu)^2\right) x, \qquad \nu = \frac{2\pi}{T},\$\$ whose null space is
the constant together with the first \\h\\ harmonics of the period
\\T\\. At `harmonics = 1` it is the harmonic acceleration operator \\L x
= \nu^2 Dx + D^3 x\\ of Ramsay and Silverman, which is what
`fda::vec2Lfd(c(0, (2*pi/T)^2, 0))` builds.

## Usage

``` r
harmonic_operator(period = NULL, harmonics = 1)
```

## Arguments

- period:

  The length of one cycle, a positive number, or `NULL` to take the
  width of the smoother's interval.

- harmonics:

  How many harmonics the penalty leaves unpenalized, a whole number of
  at least 1. The order of the operator is `2 * harmonics + 1` and the
  null space has that dimension.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Why a periodic basis wants it

The derivative operator asks a fit to contract toward a straight line,
and a straight line is not what a cyclic phenomenon simplifies to: it is
not even periodic. What a seasonal series contracts to is its
fundamental harmonic, a constant level plus one sine and one cosine of
the period, and this is the operator that leaves exactly that alone.

The difference is visible in the fit rather than only in the algebra. On
300 observations of a truth \\2 + 1.5\sin\nu t + \cos\nu t + 0.35\sin
3\nu t\\, at a smoothing parameter large enough to flatten the higher
harmonics, the fit under \\D^2\\ has fallen to a standard deviation of
0.745 and a fundamental amplitude of 1.053 against a true 1.803, while
the fit under this operator keeps 1.271 and 1.795 and is a pure sinusoid
to 3.9e-08.

## The period

`period` is \\T\\, the length of one cycle, in the units of the
covariate: 365 for a day of the year, 12 for a month, 24 for an hour. It
is **not** \\\nu = 2\pi/T\\, which some of the literature also calls
omega. `NULL` leaves it to be resolved from the interval of whichever
smoother the operator is given to, which is what
[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
does by default;
[`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md)
is what fills it in, and the accessors report an unresolved operator as
such rather than guessing.

## References

Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*,
second edition. Springer, chapter 5.

## See also

[`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md),
whose default it is;
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md)
for the non-periodic one;
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for the space it leaves.

## Examples

``` r
harmonic_operator(365)
#> <linear differential operator>  order 3
#>   L x = 0.0002963 Dx + D^3 x
#>   null space: 1, sin(0.0172142 t), cos(0.0172142 t)

# the null space is the constant and the fundamental
operator_null(harmonic_operator(365))
#>              label         rate       freq degree part
#> 1                1 0.000000e+00 0.00000000      0     
#> 2 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 3 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# two harmonics left free instead of one
harmonic_operator(365, harmonics = 2)
#> <linear differential operator>  order 5
#>   L x = 3.512e-07 Dx + 0.001482 D^3 x + D^5 x
#>   null space: 1, sin(0.0172142 t), cos(0.0172142 t), sin(0.0344284 t), cos(0.0344284 t)

# the weights are fda's: c(w0, w1, w2) = c(0, (2*pi/365)^2, 0)
round(operator_weights(harmonic_operator(365)), 8)
#> [1] 0.00000000 0.00029633 0.00000000

# the period may be left to the smoother's interval
operator_weights(harmonic_operator())
#> [1] NA NA NA
```
