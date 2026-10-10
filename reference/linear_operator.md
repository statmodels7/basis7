# A Linear Differential Operator From Its Weights

The general form: `linear_operator(w)` is \\L x = w_0 x + w_1 Dx +
\cdots + w\_{m-1} D^{m-1} x + D^m x\\, with the leading coefficient 1
and `w` given in increasing order of differentiation. It is the vector
`fda::vec2Lfd()` takes.

## Usage

``` r
linear_operator(w)
```

## Arguments

- w:

  The weights \\w_0, \ldots, w\_{m-1}\\, a numeric vector of at least
  one finite entry. Its length is the order of the operator.

## Value

An S7 object of class
[LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

## Details

The null space is read from the roots of the characteristic polynomial
\\r^m + w\_{m-1} r^{m-1} + \cdots + w_0\\: a real root \\a\\ of
multiplicity \\\mu\\ contributes \\t^i e^{at}\\ for \\i \< \mu\\, and a
complex pair \\a \pm bi\\ contributes \\t^i e^{at}\cos(bt)\\ and \\t^i
e^{at}\sin(bt)\\.
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
reports it.

Where they apply,
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
and
[`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
are preferable: they record what the operator is for, and the last two
carry a period that can be resolved from the interval of a smoother,
which a bare weight vector cannot.

## See also

[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md)
and
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
for the named instances,
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for the null space.

## Examples

``` r
# the harmonic acceleration operator written out, as in fda
nu2 <- (2 * pi / 365)^2
linear_operator(c(0, nu2, 0))
#> <linear differential operator>  order 3
#>   L x = 0.0002963 Dx + D^3 x
#>   null space: 1, sin(0.0172142 t), cos(0.0172142 t)

# and it is the same operator harmonic_operator() builds
all.equal(operator_weights(linear_operator(c(0, nu2, 0))),
          operator_weights(harmonic_operator(365)))
#> [1] TRUE

# exponential growth left unpenalized: L x = Dx - r x
operator_null(linear_operator(-0.5))
#>        label rate freq degree part
#> 1 exp(0.5 t)  0.5    0      0     
```
