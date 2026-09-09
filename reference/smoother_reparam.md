# The Coordinates a Smoother's Coefficients Live In

Applies the constraint and the reparametrization, returning the block,
the penalty matrix on it, and the basis object
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

  The basis
  [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md)
  returned.

- x:

  The covariate.

- g:

  The roughness matrix
  [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
  returned.

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
  not under a general one, so a sparse or heavy-tailed prior on a smooth
  is computationally reachable precisely because the basis is rotated
  this way.

- `"none"`:

  The constrained basis as it stands, with the penalty the congruence of
  the roughness matrix. The coefficients are the basis coefficients,
  which is what a difference penalty is written on.

- `"orthonorm"`:

  The constrained basis rotated so that it satisfies \\X'X = I\\ over
  the observed covariate. The orthonormality is against the
  **empirical** measure, which is what makes the design orthonormal;
  [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  orthonormalizes against the \\L^2\\ inner product instead and remains
  available as an operation on a basis. It is the **reparametrized
  part** that is orthonormal: with `null_space = "keep"` a free column
  is prepended afterwards and the whole block is then not orthonormal,
  while `null_space = "drop"` gives \\X'X = I\\ for the block itself,
  measured at 3.1e-15.

The three describe the same space, so an unpenalized fit cannot tell
them apart: measured at `k = 12` over 300 observations, the fitted
values of the three agree to 1.8e-15. What differs is what a coefficient
means and therefore what the penalty is.
