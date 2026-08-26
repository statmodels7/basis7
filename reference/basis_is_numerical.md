# Which of a Basis's Methods Are Numerical

Reports, for each of
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
and
[`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md),
whether the basis supplies its own method or falls back to the numerical
one on the abstract
[basis](https://statmodels7.github.io/basis7/reference/basis.md) class.
`TRUE` means the values come from finite differences or quadrature, so
they carry that method's error and cannot be verified against a
numerical reference.

## Usage

``` r
basis_is_numerical(basis)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

## Value

A named logical vector of length 3, with elements `basis_deriv`,
`basis_int` and `basis_gram`, `TRUE` where the numerical fallback is in
force.

## What it decides

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
reads it to decide which of its checks can mean anything. Comparing a
finite difference against a finite difference is the same arithmetic
twice, agreeing however wrong the basis is, so the `deriv` and
`integral` checks are not run at all where this reports `TRUE`.
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
shows the same information on its `Numerical:` line.

## How it is answered

For an ordinary class, by asking the family through
[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md),
whose default method is the owner test: which class each method is
registered on, through `attr(m, "signature")[[1]]`, tested with
[`is_base_basis_class()`](https://statmodels7.github.io/basis7/reference/is_base_basis_class.md).
A generic with no method at all counts as numerical.

Two classes are answered by delegation instead. A
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
reports its parent's: all three of its methods are registered, but each
one calls the parent and multiplies, so what is exact about it is
whatever was exact about the parent. A
[TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
reports `TRUE` for a quantity that is numerical in **any** margin, every
one of its methods being a product of its margins'.

## A method that branches

The owner test says where a method came from and not what it does, so a
method registered on a concrete class that itself calls the fallback
would be reported as exact.
[`basis_gram.FourierBasis()`](https://statmodels7.github.io/basis7/reference/basis_gram.FourierBasis.md)
does that when the period is not the interval width, and it is why
[`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md)
is a generic: a family whose route depends on its own parameters
registers a method and is believed over the owner test. A basis written
outside the package does the same.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
which reads it to decide what to test;
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md),
which prints it;
[`vignette("defining-a-basis")`](https://statmodels7.github.io/basis7/articles/defining-a-basis.md),
where a basis is taken from all three `TRUE` to all three `FALSE`.

## Examples

``` r
# Every shipped family, and every wrapper over one, is exact throughout.
basis_is_numerical(bspline_basis(dimension = 5))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 
basis_is_numerical(orthonorm_basis(fourier_basis(dimension = 5)))
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 

# A basis defined from its evaluation alone answers TRUE to all three.
Bumps <- S7::new_class("Bumps", parent = basis)
S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
  out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
                          "-")^2 / 0.12^2)
  colnames(out) <- basis_colnames(basis)
  out
}
#> Overwriting method basis_eval(<Bumps>)
bump <- Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)
basis_is_numerical(bump)
#> basis_deriv   basis_int  basis_gram 
#>        TRUE        TRUE        TRUE 

# A product is exact only where every margin is.
basis_is_numerical(tensor_basis(bump, poly_basis(dimension = 3)))
#> basis_deriv   basis_int  basis_gram 
#>        TRUE        TRUE        TRUE 
```
