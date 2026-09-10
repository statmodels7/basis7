# Compose Two Differential Operators

`L1 * L2` is the operator that applies one after the other. Composition
multiplies the characteristic polynomials, so the order of the product
is the sum of the orders and its null space is the union of the two null
spaces.

## Arguments

- e1, e2:

  Two
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md)
  objects.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Details

The union is what makes composition worth having: an operator whose
penalty leaves both a linear trend and a seasonal cycle alone is
`deriv_operator(2) * oscillator_operator(365)`, and there is no reason
to write its four weights by hand. Composition commutes, these operators
having constant coefficients.

⚠️ Compose with
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
and not with
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
unless the extra power of \\t\\ is wanted. The harmonic operator already
carries a leading \\D\\, so `deriv_operator(2) * harmonic_operator(365)`
is \\D^3(D^2 + \nu^2)\\, of order 5, whose null space is \\1, t, t^2,
\sin\nu t, \cos\nu t\\. The arithmetic is right either way and the
dimension of the null space is always the order; which of the two is
meant is the thing to read off
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
before fitting.

An operator whose period has not been resolved cannot be composed: the
weights of the product depend on the period, and resolving afterwards
would have to know which factor it came from. Give the period, or
compose after
[`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md).

## See also

[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md),
which reports the union.

## Examples

``` r
# a linear trend and a yearly cycle, both left unpenalized
op <- deriv_operator(2) * oscillator_operator(365)
op
#> <linear differential operator>  order 4
#>   L x = 0.0002963 D^2 x + D^4 x
#>   null space: 1, t, sin(0.0172142 t), cos(0.0172142 t)
operator_null(op)
#>              label         rate       freq degree part
#> 1                1 0.000000e+00 0.00000000      0     
#> 2                t 0.000000e+00 0.00000000      1     
#> 3 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 4 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# the order is the sum, the null space the union
c(operator_order(deriv_operator(2)), operator_order(oscillator_operator(365)),
  operator_order(op))
#> [1] 2 2 4

# composing with the harmonic operator instead raises the order by one
# and adds a power of t: read the null space, not the name
operator_null(deriv_operator(2) * harmonic_operator(365))$label
#> [1] "1"                "t"                "t^2"              "sin(0.0172142 t)"
#> [5] "cos(0.0172142 t)"
```
