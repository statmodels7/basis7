# Which Route a Basis's Methods Take

Reports, for each of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
whether this basis computes the quantity numerically. It is what
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
asks once it has dealt with the two wrapper classes, and it is the
generic a family overrides when its route depends on its own parameters
rather than on which class its method is registered on.

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

## Why it is a generic

The default answer reads which class each method is registered on, and
that says where a method came from rather than what it does. A method
registered on a concrete class may still call the fallback:
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
does, whenever the period is not the interval width, and the owner is
`FourierBasis` in both branches. Such a family answers for itself by
registering a method here.

distributions7 met the same defect in its expected information and
resolved it the same way. What differs is that this generic is exported,
because a basis written outside the package is an ordinary thing to
write and has the same need.

## Writing one

Take the default through
[`S7::super()`](https://rconsortium.github.io/S7/reference/super.html)
and set the entries your own branching decides:

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
would mean naming its own package.

## See also

[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md),
the predicate consumers call;
[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md)
for the default computation; and
[`basis_numerical_route.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.FourierBasis.md)
for the one family that overrides.

## Examples

``` r
# Every shipped family answers through the default method.
basis_numerical_route(bspline_basis(dimension = 5))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 

# A Fourier basis whose period is not the interval width computes its Gram
# matrix by quadrature, and says so, where reading the owner would not.
basis_numerical_route(fourier_basis(dimension = 5, omega = 0.7))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE        TRUE 
```
