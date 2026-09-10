# The Roughness Matrix of an Operator on a Fourier Basis

Exact and diagonal, at any operator, whenever the interval is a full
period. The diagonal entry of the pair at frequency \\j\\ is \$\$\lvert
P(i\nu_j)\rvert^2 \\ \frac{T}{2}, \qquad \nu_j = \frac{2\pi j}{T},\$\$
with \\P(r) = r^m + \sum\_{k\<m} w_k r^k\\ the operator's characteristic
polynomial, and the entry of the constant column is \\w_0^2 T\\.

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

  Passed on to the quadrature (`panels`, `nodes`).

## Value

A diagonal numeric matrix of `basis@dimension` rows and columns.

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

⚠️ The book this operator comes from states that the harmonic penalty
makes \\R\\ structurally different and **more complex** than the
diagonal matrix the derivative penalty gives. The first half is right
and the second is not: measured on a nine-function basis, the largest
off-diagonal entry is 1.7e-16 of the largest entry, and it stays that
way for a two-harmonic operator and for a composed one. What changes is
the null space, which goes from the constant alone to the constant and
the harmonics the operator keeps.

The operator's own period does not have to be the basis's. Only the
basis's enters the orthogonality; the operator's enters through its
weights, and an operator tuned to another cycle simply gives a diagonal
with no exact zeros.

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
