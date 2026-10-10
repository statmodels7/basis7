# Which Route a Basis's Methods Take

Reports, for each of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
whether this basis computes the quantity numerically.
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
calls it once it has handled the two wrapper classes, and a family
overrides it when its route depends on its own parameters instead of the
class on which its method is registered.

## Usage

``` r
basis_numerical_route(basis, ...)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- ...:

  Unused, and accepted so a method's signature can match.

## Value

A named logical vector of length 3, with elements `basis_deriv`,
`basis_int` and `basis_gram`, `TRUE` where the quantity is computed
numerically.

## Purpose of the generic

The default method reads the class on which each method is registered,
which records where a method came from and not what it does. A method
registered on a concrete class may still call the fallback:
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
does, whenever the period is not the interval width, and the owner is
`FourierBasis` in both branches. Such a family reports its own route by
registering a method here. The generic is exported so that a basis
written outside the package can do the same.

The two wrapper classes,
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
and
[TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md),
are handled by
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
before this generic is reached; called directly on a wrapper, the
default method reports the wrapper's own methods.

## Writing one

A method takes the default through
[`S7::super()`](https://rconsortium.github.io/S7/reference/super.html)
and sets the entries decided by its own branching:

    S7::method(basis_numerical_route, MyBasis) <- function(basis, ...) {
      out <- basis_numerical_route(S7::super(basis, basis7::basis))
      out[["basis_gram"]] <- !my_closed_form_applies(basis)
      out
    }

The package's own override calls
[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md)
instead, which is that default under a name: inside a method the formal
`basis` shadows the class of the same name, so reaching the class
through
[`S7::super()`](https://rconsortium.github.io/S7/reference/super.html)
would require naming the package inside itself.

## See also

[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md),
the predicate consumers call;
[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md)
for the default computation; and
[`basis_numerical_route.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.FourierBasis.md)
for the one family that overrides.

## Examples

``` r
# A B-spline basis uses the default method.
basis_numerical_route(bspline_basis(dimension = 5))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 

# A Fourier basis whose period is not the interval width computes its Gram
# matrix by quadrature, which its own method reports and the owner test
# would not.
basis_numerical_route(fourier_basis(dimension = 5, omega = 0.7))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE        TRUE 
```
