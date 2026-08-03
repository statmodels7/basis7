# Package index

## Bases

The families the package ships. Each is complete: all its functions are
kept, so restricting one is a separate, deliberate step.

- [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
  : Construct a B-Spline Basis
- [`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md)
  : Construct a Fourier Basis
- [`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
  : Construct a Legendre Polynomial Basis

## Transforming a basis

Orthonormalisation, identifiability constraints and the Demmler-Reinsch
construction are the same operation with different matrices, so they
share one class.

- [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  : Orthonormalise a Basis
- [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  : Restrict a Basis to a Linear Constraint
- [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  : Demmler-Reinsch Basis
- [`TransformedBasis()`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  : Linearly Transformed Basis

## The interface

What every basis answers. A subclass must implement only the first; the
rest have numerical methods on the base class.

- [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
  : Evaluate a Basis
- [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  : Differentiate a Basis
- [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  : Integrate a Basis
- [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  : Gram Matrix of a Basis
- [`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
  : Column Names of a Basis Matrix
- [`basis()`](https://statmodels7.github.io/basis7/reference/basis.md) :
  Basis Expansion

## Classes

- [`BsplineBasis()`](https://statmodels7.github.io/basis7/reference/BsplineBasis.md)
  : B-Spline Basis
- [`FourierBasis()`](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  : Fourier Basis
- [`PolyBasis()`](https://statmodels7.github.io/basis7/reference/PolyBasis.md)
  : Legendre Polynomial Basis

## Validation

- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  : Validate a Basis
- [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
  : Which of a Basis's Methods Are Numerical

## Internals

Exported for anyone writing a basis of their own: the numerical
machinery the fallbacks are built from, documented rather than hidden.

- [`print.basis`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  : Print a Basis
- [`plot.basis`](https://statmodels7.github.io/basis7/reference/plot.basis.md)
  : Plot a Basis
- [`fd_weights()`](https://statmodels7.github.io/basis7/reference/fd_weights.md)
  : Finite-Difference Weights for an Arbitrary Stencil
- [`fd_offsets()`](https://statmodels7.github.io/basis7/reference/fd_offsets.md)
  : Stencil Offsets for a Derivative Order
- [`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
  : Numerically Differentiate a Matrix-Valued Function
- [`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
  : Gram Matrix by Composite Quadrature
- [`gauss_legendre()`](https://statmodels7.github.io/basis7/reference/gauss_legendre.md)
  : Gauss-Legendre Nodes and Weights
- [`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md)
  : Map a Quadrature Rule onto Intervals
- [`bspline_design()`](https://statmodels7.github.io/basis7/reference/bspline_design.md)
  : Call splines2 for a B-Spline Design Matrix
- [`fourier_trig()`](https://statmodels7.github.io/basis7/reference/fourier_trig.md)
  : The Trigonometric Columns of a Fourier Basis
- [`legendre_table()`](https://statmodels7.github.io/basis7/reference/legendre_table.md)
  : The Legendre Polynomials by Recurrence
- [`legendre_deriv()`](https://statmodels7.github.io/basis7/reference/legendre_deriv.md)
  : Derivatives of the Legendre Polynomials
- [`new_transformed()`](https://statmodels7.github.io/basis7/reference/new_transformed.md)
  : Build a Transformed Basis
- [`empirical_gram()`](https://statmodels7.github.io/basis7/reference/empirical_gram.md)
  : Gram Matrix Against the Empirical Measure
- [`weighted_gram()`](https://statmodels7.github.io/basis7/reference/weighted_gram.md)
  : Gram Matrix Against a Weighted Lebesgue Measure
- [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  : A Finite-Difference Reference, and Where It Can Be Trusted
- [`check_basis_args()`](https://statmodels7.github.io/basis7/reference/check_basis_args.md)
  : Validate the Arguments Every Basis Constructor Takes
- [`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md)
  : Validate Evaluation Points Against a Basis
- [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md)
  : Validate a Derivative Order
- [`name_columns()`](https://statmodels7.github.io/basis7/reference/name_columns.md)
  : Name the Columns of a Basis Matrix
- [`basis_partitions_unity()`](https://statmodels7.github.io/basis7/reference/basis_partitions_unity.md)
  : Does This Basis Sum to One?
- [`rel_close()`](https://statmodels7.github.io/basis7/reference/rel_close.md)
  : Compare Two Matrices Relative to Their Own Magnitude
- [`print_basis_checks()`](https://statmodels7.github.io/basis7/reference/print_basis_checks.md)
  : Print the Outcome of check_basis
- [`is_base_basis_class()`](https://statmodels7.github.io/basis7/reference/is_base_basis_class.md)
  : Is This the Package's Own Base Class?
- [`basis7`](https://statmodels7.github.io/basis7/reference/basis7-package.md)
  [`basis7-package`](https://statmodels7.github.io/basis7/reference/basis7-package.md)
  : basis7: An S7 Framework for Basis Expansions
