# Contract an Ordinary Basis Against Coefficients

The design matrix times the coefficients. For a basis of one variable
there is nothing to avoid forming, so the definition is the computation.

## Arguments

- basis:

  An object inheriting from class `basis`.

- x:

  A numeric vector of evaluation points.

- coef:

  A numeric vector, or a matrix of several sets of coefficients.

- ...:

  Unused.

## Value

A numeric vector, or a matrix when `coef` is one.
