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
necessary: [`identical()`](https://rdrr.io/r/base/identical.html) on an
S7 class is object identity, so it is `FALSE` for a class re-created
from the same definition, and that is what happens whenever a package's
code is evaluated instead of loaded. Coverage tools do exactly that, so
an identity-only test passes every ordinary check and fails under `covr`
alone.

The same defect in `linkfunctions7` made every fallback differentiate
the order below it, and the log link's fourth derivative came back wrong
by a factor of 900 while the whole five-platform check matrix stayed
green.

## See also

[`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md),
its only caller, and
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md),
which that answers for.
