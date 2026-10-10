# The Roughness Matrix of an Operator on a Fourier Basis

Exact and diagonal, at any operator, when the interval is a full period
and neither `at` nor `weight` is given; otherwise the matrix is computed
by
[`numerical_operator_gram()`](https://statmodels7.github.io/basis7/reference/numerical_operator_gram.md)
and is in general full. The diagonal entry of the pair at frequency
\\j\\ is \$\$\lvert P(i\nu_j)\rvert^2 \\ \frac{T}{2}, \qquad \nu_j =
\frac{2\pi j}{T},\$\$ with \\P(r) = r^m + \sum\_{k\<m} w_k r^k\\ the
operator's characteristic polynomial, and the entry of the constant
column is \\w_0^2 T\\.

## Arguments

- basis:

  A [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- at:

  The covariate values, for the empirical measure, or `NULL`.

- weight:

  A density to integrate against, or `NULL`.

- ...:

  Passed on to the quadrature (`panels`, `nodes`) on the numerical
  route. On the exact route `panels` or `nodes` signals an error.

## Value

A numeric matrix of `basis@dimension` rows and columns, diagonal on the
exact route.

## Details

The identity behind it is the phase shift: \\D^k \sin(\nu_j t) = \nu_j^k
\sin(\nu_j t + k\pi/2)\\, so the even derivatives return a sine and the
odd ones a cosine, and collecting them gives \$\$L \sin(\nu_j t) =
\mathrm{Re}\\P(i\nu_j)\\\sin(\nu_j t) +
\mathrm{Im}\\P(i\nu_j)\\\cos(\nu_j t),\$\$ with \\L\cos\\ the same
rotation the other way. The two have equal squared norms and are
orthogonal to each other, and over a full period everything at one
frequency is orthogonal to everything at another, so the matrix is
diagonal.

The matrix is diagonal for a harmonic operator, as it is for a
derivative penalty, at any number of harmonics and for a composed
operator. The two differ in the null space, which is the constant alone
for a derivative penalty and the constant together with the retained
harmonics for the harmonic operator.

The period of the operator need not be that of the basis. Only the
period of the basis enters the orthogonality, and the period of the
operator enters through its weights. An operator tuned to another cycle
gives a diagonal with no zero at the sines and cosines; the entry of the
constant is zero whenever \\w_0 = 0\\, as for
[`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md).

## See also

[`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
for the generic and the numerical route,
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
for the derivative case.

## Examples

``` r
b <- fourier_basis(lower = 0, upper = 365, dimension = 9)
g <- basis_gram(b, order = harmonic_operator(365))

# diagonal, and zero on the constant and the fundamental
max(abs(g - diag(diag(g)))) / max(abs(g))
#> [1] 0
round(diag(g)[1:3], 10)
#> const  sin1  cos1 
#>     0     0     0 

# and it agrees with integrating (Lb)(Lb)' numerically
gn <- basis7:::numerical_operator_gram(b, harmonic_operator(365),
                                       panels = 200L)
max(abs(g - gn)) / max(abs(g))
#> [1] 1.98186e-16
```
