# Integral of a Tensor Product Basis

Returns the integral over the box from the lower corner to each point,
one iterated integral per variable, in closed form wherever the margins
have one. It is the only integral over a box in the package, the
numerical fallback taking one variable alone.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable, the upper corners of
  the boxes.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `nrow(x)` rows and `basis@dimension` columns,
exactly zero in the row at the lower corner.

## Details

The integrand separates, so the multiple integral is the product of the
marginal integrals. The anchoring convention of
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
survives without a correction: a product in which every factor is zero
at the corner is zero at the corner, so the row at `basis@lower` is
exactly zero.

## See also

[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md),
which does the work;
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the generic and the anchoring convention.
