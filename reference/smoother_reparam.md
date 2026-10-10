# The Coordinates of the Coefficients of a Smoother

Applies the constraint and the reparametrization, returning the penalty
matrix on the block and the basis object that
[`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
reapplies.

## Usage

``` r
smoother_reparam(sm, b, x, g, cons)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- b:

  The basis returned by
  [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md).

- x:

  The covariate.

- g:

  The roughness matrix returned by
  [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md).

- cons:

  The constraint, one column per direction removed, or `NULL`.

## Value

A list of `S`, the penalty on the reparametrized block, and `basis`, the
object whose evaluation **is** that block and which
[`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
evaluates at new values.

## Details

The three routes differ in what the coefficients mean, and the fit they
express is the same span in each.

- `"dr"`:

  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  diagonalizes the pencil of the empirical Gram matrix against the
  roughness matrix, so the columns are orthogonal over the data, ordered
  from the smoothest to the most oscillatory, and the penalty is the
  identity. A separable penalty is available under a diagonal map and
  not under a general one, so the rotation makes a sparse or
  heavy-tailed prior on a smooth computationally feasible.

- `"none"`:

  The constrained basis as it stands, with the penalty the congruence of
  the roughness matrix by the transform of the constraint. The
  coefficients are those of the constrained basis, one fewer than the
  basis functions for each direction removed.

- `"orthonorm"`:

  The constrained basis rotated so that it satisfies \\X'X = I\\ over
  the observed covariate. The orthonormality is against the
  **empirical** measure, which makes the design orthonormal;
  [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  orthonormalizes against the \\L^2\\ inner product instead and remains
  available as an operation on a basis. It is the **reparametrized
  part** that is orthonormal: with `null_space = "keep"` a free column
  is prepended afterwards and the whole block is then not orthonormal,
  while `null_space = "drop"` gives \\X'X = I\\ for the block itself, up
  to rounding.

The three describe the same space, so the fitted values of an
unpenalized fit agree to rounding. They differ in what a coefficient
means, and therefore in the penalty.
