# Contract an Ordinary Basis Against Coefficients

Returns `basis_eval(basis, x) %*% coef`. For a basis of one variable the
design matrix has `basis@dimension` columns and there is nothing worth
avoiding, so the definition is the computation. The method every class
inherits except
[TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md).

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md) other
  than
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md).

- x:

  A numeric vector of evaluation points inside the basis interval.

- coef:

  A numeric vector of length `basis@dimension`, or a matrix with that
  many rows and one column per set of coefficients.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector with one value per point, or a matrix with one column
per column of `coef`.

## Details

`coef` is taken in the design's own column order, there being one
obvious one. Several sets of coefficients may be given as the columns of
a matrix, and the result then has one column each.

## See also

[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
for the generic;
[`basis_contract.TensorBasis()`](https://statmodels7.github.io/basis7/reference/basis_contract.TensorBasis.md)
for the case where the design matrix is worth avoiding.
