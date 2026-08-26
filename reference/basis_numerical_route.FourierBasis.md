# The Gram Route of a Fourier Basis

Reports `basis_gram` as `TRUE` when `basis_params$full_period` is
`FALSE`, where the owner test would read `FourierBasis` and answer
`FALSE`.
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
delegates to
[`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
in that case, so the matrix is a composite Gauss-Legendre quadrature and
carries its error. The evaluation, the derivatives and the anchored
integral are closed form at any period and are left as the owner test
finds them.

What this buys is that
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
holds the Gram matrix to the tolerance a quadrature deserves rather than
the one meant for a closed form, and that
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
names the route in use.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

The named logical vector
[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md)
describes, with `basis_gram` `TRUE` for a basis whose period is not the
interval width.

## See also

[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
for the two branches, and
[`basis_numerical_route.basis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.basis.md)
for the default this starts from.
