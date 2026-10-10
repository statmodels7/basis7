# Default Numerical Route by Owner Test

Returns
[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md):
the class on which each of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
is registered. That is correct for every family whose methods take one
route, which excludes the Fourier family, and for a class that is not a
wrapper;
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
handles the wrappers itself.

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
