# The Oscillator Operator

The periodic factor on its own, \$\$L x = \prod\_{i=1}^{h}\left(D^2 +
(i\nu)^2\right) x, \qquad \nu = \frac{2\pi}{T},\$\$ whose null space is
the first \\h\\ harmonics of the period \\T\\ **without** the constant.
At `harmonics = 1` it is the equation of simple harmonic motion, \\L x =
\nu^2 x + D^2 x\\.

## Usage

``` r
oscillator_operator(period = NULL, harmonics = 1)
```

## Arguments

- period:

  The length of one cycle, a positive number, or `NULL` to take the
  width of the smoother's interval.

- harmonics:

  How many harmonics the penalty leaves unpenalized, a whole number of
  at least 1. The order of the operator is `2 * harmonics`.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Details

It exists to be composed.
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
is this operator with a leading \\D\\, which is what puts the constant
into the null space, so the two are related by `harmonic_operator(T, h)`
being `deriv_operator(1) * oscillator_operator(T, h)` exactly, and a
penalty that should leave a **linear trend** and a cycle alone rather
than a level and a cycle is
`deriv_operator(2) * oscillator_operator(T)`.

Composing with
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
instead would raise the order by one and put an extra power of \\t\\ in
the null space, which is why the factor is offered separately rather
than left to be written out by hand.

## See also

[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
which is this with a leading derivative;
[operator-compose](https://statmodels7.github.io/basis7/reference/operator-compose.md)
for what composing does to the null space.

## Examples

``` r
oscillator_operator(365)
#> <linear differential operator>  order 2
#>   L x = 0.0002963 x + D^2 x
#>   null space: sin(0.0172142 t), cos(0.0172142 t)

# the null space is the fundamental alone: no constant
operator_null(oscillator_operator(365))
#>              label         rate       freq degree part
#> 1 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 2 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# a linear trend and a yearly cycle, and nothing else
op <- deriv_operator(2) * oscillator_operator(365)
operator_null(op)
#>              label         rate       freq degree part
#> 1                1 0.000000e+00 0.00000000      0     
#> 2                t 0.000000e+00 0.00000000      1     
#> 3 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 4 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# the harmonic operator is this one with a leading derivative
a <- operator_weights(deriv_operator(1) * oscillator_operator(365))
b <- operator_weights(harmonic_operator(365))
all.equal(a, b)
#> [1] TRUE
```
