# Is This the Package's Own Base Class?

Reports whether an S7 class is the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class.
That is how
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
tells a method registered on the base class, which is a numerical
fallback, from one a subclass supplied.

## Usage

``` r
is_base_basis_class(cls)
```

## Arguments

- cls:

  An S7 class object, as returned by
  [`S7::S7_class()`](https://rconsortium.github.io/S7/reference/S7_class.html).

## Value

A single `TRUE` or `FALSE`.

## Details

Identity is tried first, being the usual case and costing nothing, and
the name and package are compared when it fails. The second test is
needed because [`identical()`](https://rdrr.io/r/base/identical.html)
can return `FALSE` for a class re-created from the same definition,
which happens when the code of a package is evaluated instead of loaded,
as coverage tools do.

## See also

[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md),
its only caller, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md),
which relies on it.
