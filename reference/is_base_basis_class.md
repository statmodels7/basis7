# Is This the Package's Own Base Class?

Asks whether an S7 class is the abstract
[`basis`](https://statmodels7.github.io/basis7/reference/basis.md)
class, which is how a method registered on it is told apart from one a
subclass supplied.

## Usage

``` r
is_base_basis_class(cls)
```

## Arguments

- cls:

  An S7 class.

## Value

`TRUE` or `FALSE`.

## Details

Identity is tried first because it is the usual case and costs nothing,
then the name and package, because identity is not preserved when a
package's code is re-evaluated rather than loaded. Coverage tools do
exactly that, so an identity-only test passes every ordinary check and
fails there.
