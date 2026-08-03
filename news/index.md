# Changelog

## basis7 0.3.1

- The numerical derivative no longer labels the rows of its result after
  the finite-difference stencil each point uses. The stencil is chosen
  per point and recorded in a character vector, whose names travelled
  through the evaluation points and came back as row names for any basis
  whose method propagates the names of `x`.

- Whether a Gram matrix or a penalty is positive definite is now decided
  from its eigenvalues rather than from whether
  [`chol()`](https://rdrr.io/r/base/chol.html) raises. On a matrix with
  an exactly zero eigenvalue the pivot that should be zero comes out
  positive or negative according to rounding, so
  [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  and
  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  gave different answers on different machines about the same matrix.

## basis7 0.3.0

- [`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
  the product of bases, one per variable. Everything follows from the
  marginals because the product separates: a partial derivative
  differentiates one marginal and leaves the others, the integral over
  the box from its lower corner is the product of the marginal
  integrals, and the Gram matrix is the Kronecker product of the
  marginal ones. A tensor of exactly integrated marginals is therefore
  exact at any number of variables, where a quadrature over the box
  would not be.

- A basis now declares how many variables it takes, through
  [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  and the length of its endpoints, so a univariate basis is the case of
  one variable rather than a separate kind of object. Evaluation points
  become a matrix with one column per variable, and the derivative order
  becomes a multi-index. A single non-zero order is refused for a
  product, since it could mean that order in every coordinate or that
  total order.

- [`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
  evaluates a basis against coefficients without necessarily forming the
  design matrix. Coefficients come as a full array, processed in row
  blocks so the peak memory is bounded by the block rather than by the
  sample, or as a list of factor matrices in canonical polyadic form,
  where the cost is linear in the number of variables instead of
  exponential in it. A product of six bases of ten functions has a
  million columns; the factorised route touches neither that matrix nor
  the coefficient array.

- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
  [`print()`](https://rdrr.io/r/base/print.html) and the numerical
  fallbacks understand several variables.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) refuses a
  product and says why.

## basis7 0.2.0

- [`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md),
  the Legendre polynomials by recurrence. They span the same space as
  the raw powers and are chosen over them for conditioning: the Gram
  matrix of raw powers is a Hilbert matrix, and ten of them are already
  close to singular in double precision.

- [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  takes a measure. The default is Lebesgue on the interval; `at` gives
  the empirical measure of a sample, which is the matrix a design matrix
  produces, and `weight` gives a weighted integral. Both are handled in
  the body of the generic, before dispatch, so a method never implements
  them. **Breaking for user-written methods**: a
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  method must now name `at` and `weight` in its signature, because S7
  requires a method’s formals to contain the generic’s.

- `TransformedBasis`, one class for every linear reparametrisation ,
  with three constructors:
  [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
  from the Cholesky factor of the Gram matrix rather than from a grid;
  [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md),
  from the null space of a constraint; and
  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
  the Demmler-Reinsch construction, which diagonalises the empirical
  inner product and the penalty at once and is empirically orthogonal to
  a constant and to the covariate.

  Transforms compose by multiplication rather than by nesting, and a
  transformed basis reports the *parent’s* numerical status, since its
  own methods delegate and multiply.

- [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  factorises the penalty and not the design, so it survives a
  rank-deficient design, which equally spaced knots produce whenever the
  data leave a knot span empty. Verified on a design where the Cholesky
  factor of the design matrix does not exist.

- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  allows for the accuracy of its own reference. Each numerical reference
  is computed at a step and at half of it, and the gap between them
  bounds its uncertainty; the comparison is given that much slack, point
  by point. Without it a spline failed at its own knots, where the third
  derivative jumps and a central difference returns the jump rather than
  the truncation error. A deliberate five per cent error is still caught
  by four orders of magnitude.

- A vignette, `defining-a-basis`.

## basis7 0.1.0

First release.

- The `basis` class and four generics:
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
  [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md).
  Derivative order is an argument, so no order is privileged, and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  returns the inner products a roughness penalty integrates.

- Two families with exact formulas:
  [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
  whose Gram matrix is integrated exactly knot interval by knot
  interval, and
  [`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
  whose derivatives of every order and whose antiderivative come from
  one phase-shift identity, and whose Gram matrix is diagonal in closed
  form over a whole period.

- Numerical fallbacks registered on the base class, so a basis that
  implements only
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
  is complete. Derivatives use one finite-difference stencil built from
  a Vandermonde solve, never a chain of first differences, and switch to
  a one-sided stencil at the interval endpoints.

- [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  is anchored: its value at the lower endpoint is exactly zero, for
  every family. Any antiderivative satisfies a differentiation check, so
  the constant is fixed by the contract instead.

- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  runs six numerical checks and reports a quantity that came from a
  fallback as unchecked rather than as passed, and a property the family
  never claimed as not claimed.

- [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
  answers the same question programmatically.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods; the
  plot draws the basis, any derivative, or the integral.
