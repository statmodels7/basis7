# Contract a Tensor Product Basis Against Coefficients

Returns the value of the function the coefficients describe, computed
from the marginal evaluations. A list of factor matrices goes to
[`contract_cp()`](https://statmodels7.github.io/basis7/reference/contract_cp.md),
which forms nothing; an array or a vector is contracted in blocks of
`block` rows, so the peak memory is one block's design matrix and not
the whole one.

## Arguments

- basis:

  A
  [TensorBasis](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  object.

- x:

  A numeric matrix with one column per variable.

- coef:

  An array of dimension `(K_1, ..., K_D)`, a numeric vector of
  `basis@dimension` values in the design's own column order, or a list
  of `D` factor matrices in canonical polyadic form.

- block:

  The number of rows processed at once when `coef` is an array or a
  vector, default `1024`. It bounds the peak memory, which is otherwise
  what forming the design matrix would cost. Ignored for a list.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector with one value per row of `x`.

## The blocked route

Each block of rows is evaluated through
[`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md)
and multiplied by the flattened coefficients, and the block is then
discarded. Measured at four margins of eight functions and 20000 points,
where the full design is 625 MB: 32 MB at the default block and 0.83 s,
against 625 MB and 1.00 s for the design route, agreeing exactly.
Raising `block` past a few thousand buys no speed and costs memory
linearly.

## The flattening

An array is passed through
[`aperm()`](https://rdrr.io/r/base/aperm.html) before flattening,
because the design's columns run with the last margin fastest while an R
array is stored with its first index fastest. A plain vector has no
dimensions to reverse and is taken in the design's column order as it
stands, so the two shapes describe different functions from the same
numbers. See
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md).

A `coef` of the wrong length, or an array whose dimensions are not the
margins', throws with the expected values named.

## See also

[`contract_cp()`](https://statmodels7.github.io/basis7/reference/contract_cp.md)
for the factorized route;
[`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
for the generic and the two shapes of coefficient.
