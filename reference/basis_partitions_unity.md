# Does This Basis Sum to One?

Reports whether the family is a partition of unity, so that
[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
tests the row sums only where the property is claimed. A family that is
not one would fail a check it never promised to pass.

## Usage

``` r
basis_partitions_unity(basis)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

## Value

A single `TRUE` or `FALSE`.

## Details

`TRUE` for a
[BsplineBasis](https://statmodels7.github.io/basis7/reference/BsplineBasis.md),
which sums to one by construction, and for a
[TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
all of whose margins do: the row sums of a Kronecker product are the
products of the row sums.

`FALSE` for everything else, including a
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
over a B-spline. That is deliberate: an orthonormalization or a
constraint takes linear combinations of the columns, and the sum of the
new columns is generally not one. It is also conservative, so a
transformation that happens to preserve the property is untested rather
than wrongly failed.

## See also

[`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
its only caller.
