# The Owner Test Behind the Default Route

Reports which of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
this basis takes from the numerical fallback, by asking which class each
method is registered on through `attr(m, "signature")[[1]]` and testing
it with
[`is_base_basis_class()`](https://statmodels7.github.io/basis7/reference/is_base_basis_class.md).
A generic with no method at all counts as numerical, the base class
being where the fallback lives.

## Usage

``` r
route_by_owner(basis)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

## Value

A named logical vector of length 3, elements `basis_deriv`, `basis_int`
and `basis_gram`.

## Details

It is the body of
[`basis_numerical_route.basis()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.basis.md),
kept under a name so that the package's own override can start from it.
Inside a method the formal `basis` shadows the class of the same name,
so reaching that default through
[`S7::super()`](https://rconsortium.github.io/S7/reference/super.html)
would mean naming the package inside itself; a basis written elsewhere
has no such difficulty and does use
[`S7::super()`](https://rconsortium.github.io/S7/reference/super.html).

What it cannot see is which branch a method takes once it has been
reached: that is the question
[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md)
exists to let a family answer.

## See also

[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md),
the generic it is the default of, and
[`is_base_basis_class()`](https://statmodels7.github.io/basis7/reference/is_base_basis_class.md),
which it tests each owner with.
