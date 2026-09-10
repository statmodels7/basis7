# The Null Basis of a Constraint

An orthonormal basis of the null space of `cm`, the columns a
constrained basis is built from. It is the transform **relative to the
basis being constrained**, which is not always the one the resulting
object stores.

## Usage

``` r
constraint_null(cm, dimension, tol = 1e-10)
```

## Arguments

- cm:

  The constraint, one row per direction removed and one column per basis
  function.

- dimension:

  The number of basis functions.

- tol:

  The relative tolerance at which a singular value counts as zero.

## Value

A numeric matrix of `dimension` rows and one column per direction that
survives the constraint.

## Details

[`new_transformed()`](https://statmodels7.github.io/basis7/reference/new_transformed.md)
flattens a nested transform, so a basis that is itself transformed comes
back carrying the product against its own parent: constraining a cyclic
smoother's twelve periodic functions, which are a transform of fifteen
B-splines, gives an object whose transform is 15 by 11 and not 12 by 11.
That is right for evaluating it and wrong for carrying a matrix defined
on the twelve, which is why the local transform is available here rather
than read off the object.

## See also

[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
and
[`smoother_reparam()`](https://statmodels7.github.io/basis7/reference/smoother_reparam.md),
the two callers.
