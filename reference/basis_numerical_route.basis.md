# The Owner Test, Which Is the Default Route

Answers with
[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md):
which class each of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
is registered on. That is right for every family whose methods take one
route, which is all of them but
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md).

## Arguments

- basis:

  A basis object.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

The named logical vector
[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md)
describes.

## See also

[`basis_numerical_route.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.FourierBasis.md),
the one family that needs more than this.
