# Build a Transformed Basis

Wraps a basis in a
[`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md),
collapsing the transform into the parent's when the parent is already
one.

## Usage

``` r
new_transformed(basis, transform, name, prefix, params = list())
```

## Arguments

- basis:

  The basis to transform.

- transform:

  The matrix \\T\\.

- name:

  The name of the resulting basis.

- prefix:

  The prefix for its column names.

- params:

  Extra entries for `basis_params`.

## Value

An object of class
[`TransformedBasis`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md).
