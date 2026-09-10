# The Null Space of a Differential Operator

The functions a maximal penalty leaves untouched, that is the solutions
of \\L x = 0\\. They are what a strongly penalized fit contracts to, so
this is the one thing to look at before choosing an operator.

## Usage

``` r
operator_null(op, tol = 1e-06)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md).

- tol:

  The relative tolerance at which two roots count as one, a positive
  number.

## Value

A data frame of one row per function in the null space, with columns
`label` (how the function reads), `rate` (the real part of the root),
`freq` (the imaginary part) and `degree` (the power of \\t\\ multiplying
it).

## Details

The null space is read from the roots of the characteristic polynomial
\\r^m + w\_{m-1}r^{m-1} + \cdots + w_0\\, computed by
[`base::polyroot()`](https://rdrr.io/r/base/polyroot.html) and grouped
by multiplicity. A real root \\a\\ of multiplicity \\\mu\\ contributes
\\t^i e^{at}\\ for \\i \< \mu\\, and a complex pair \\a \pm bi\\ of
multiplicity \\\mu\\ contributes both \\t^i e^{at}\cos(bt)\\ and \\t^i
e^{at}\sin(bt)\\.

The null space is a property of the **operator** and is computed from
it, never from the rank of an assembled penalty matrix. The two are
different questions and they give different answers: the same
two-harmonic operator has a penalty of null dimension 5 on a Fourier
basis and 3 on a cubic B-spline at a relative tolerance of 1e-10,
because a spline represents a sine only approximately, while the
operator's own null space is five-dimensional in both cases. Which of
those functions a **basis** can carry is a separate question, and
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
is where it is asked.

## See also

[`operator_null_design()`](https://statmodels7.github.io/basis7/reference/operator_null_design.md)
for those functions evaluated,
[`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md)
and
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
for the two named instances.

## Examples

``` r
# a straight line
operator_null(deriv_operator(2))
#>   label rate freq degree part
#> 1     1    0    0      0     
#> 2     t    0    0      1     

# a constant plus the fundamental of a yearly cycle
operator_null(harmonic_operator(365))
#>              label         rate       freq degree part
#> 1                1 0.000000e+00 0.00000000      0     
#> 2 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 3 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# both at once
operator_null(deriv_operator(2) * harmonic_operator(365))
#>              label         rate       freq degree part
#> 1                1 0.000000e+00 0.00000000      0     
#> 2                t 0.000000e+00 0.00000000      1     
#> 3              t^2 0.000000e+00 0.00000000      2     
#> 4 sin(0.0172142 t) 1.577722e-30 0.01721421      0  sin
#> 5 cos(0.0172142 t) 1.577722e-30 0.01721421      0  cos

# the dimension is always the order
nrow(operator_null(harmonic_operator(365, harmonics = 3)))
#> [1] 7
```
