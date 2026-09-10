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

Orthonormalization, identifiability constraints and the Demmler-Reinsch
construction are the same operation with different matrices, so they
share one class.

- [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  : Orthonormalize a Basis
- [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  : Restrict a Basis to a Linear Constraint
- [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  : Demmler-Reinsch Basis
- [`TransformedBasis()`](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  : Linearly Transformed Basis

## Smoothers

A smoother carries the four decisions a penalized smooth is made of –
the basis, the penalty, the null space and the reparametrization – as
one object, because the null space is a property of the basis and the
penalty together rather than of either alone. Built at the covariate, it
returns the block and the penalty matrix.

- [`smoother()`](https://statmodels7.github.io/basis7/reference/smoother.md)
  : A Smoother: the Four Decisions of a Penalized Smooth
- [`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
  : A B-Spline Smoother
- [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
  : A Fourier Smoother
- [`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md)
  : A Legendre Smoother
- [`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md)
  : A Cyclic Smoother
- [`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
  : A P-spline Smoother
- [`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md)
  : An Adaptive Smoother
- [`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
  : Build a Smoother at Data
- [`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
  : A Smoother's Block at New Values
- [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md)
  : The Basis a Smoother Builds On
- [`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
  [`smoother_span_apply()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
  : What a Smoother Removes and What It Gives Back
- [`BsplineSmoother()`](https://statmodels7.github.io/basis7/reference/BsplineSmoother.md)
  : The B-Spline Smoother Class
- [`FourierSmoother()`](https://statmodels7.github.io/basis7/reference/FourierSmoother.md)
  : The Fourier Smoother Class
- [`LegendreSmoother()`](https://statmodels7.github.io/basis7/reference/LegendreSmoother.md)
  : The Legendre Smoother Class
- [`CyclicSmoother()`](https://statmodels7.github.io/basis7/reference/CyclicSmoother.md)
  : The Cyclic Smoother Class
- [`PsplineSmoother()`](https://statmodels7.github.io/basis7/reference/PsplineSmoother.md)
  : The P-spline Smoother Class
- [`AdaptiveSmoother()`](https://statmodels7.github.io/basis7/reference/AdaptiveSmoother.md)
  : The Adaptive Smoother Class

## Several variables

A product of bases, and the contraction that computes what a fit needs
from the marginal evaluations without forming the product.

- [`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
  : Construct a Tensor Product Basis
- [`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
  : Evaluate a Basis Against Coefficients
- [`TensorBasis()`](https://statmodels7.github.io/basis7/reference/TensorBasis.md)
  : Tensor Product Basis
- [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  : How Many Variables a Basis Takes

## Differential operators

What a penalty measures. A smoother’s `order` takes one of these in
place of a whole number, and the operator carries its own null space:
the functions a maximal penalty leaves untouched.

- [`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md)
  : The Derivative Operator
- [`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md)
  : The Harmonic Acceleration Operator
- [`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
  : The Oscillator Operator
- [`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md)
  : A Linear Differential Operator From Its Weights
- [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
  : The Null Space of a Differential Operator
- [`operator_null_design()`](https://statmodels7.github.io/basis7/reference/operator_null_design.md)
  : The Null-Space Functions Evaluated
- [`operator_order()`](https://statmodels7.github.io/basis7/reference/operator_order.md)
  : The Order of a Differential Operator
- [`operator_weights()`](https://statmodels7.github.io/basis7/reference/operator_weights.md)
  : The Weights of a Differential Operator
- [`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md)
  : Fill In an Operator's Period From an Interval
- [`LinearOperator()`](https://statmodels7.github.io/basis7/reference/LinearOperator.md)
  : The Linear Differential Operator Class

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
- [`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
  : The Roughness Matrix of a Differential Operator
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
- [`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md)
  : Which Route a Basis's Methods Take

## Internals

Exported for anyone writing a basis of their own: the numerical
machinery the fallbacks are built from, documented rather than hidden.

- [`print.basis`](https://statmodels7.github.io/basis7/reference/print.basis.md)
  : Print a Basis
- [`plot.basis`](https://statmodels7.github.io/basis7/reference/plot.basis.md)
  : Plot a Basis
- [`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
  : Numerically Differentiate a Matrix-Valued Function
- [`numerical_gram()`](https://statmodels7.github.io/basis7/reference/numerical_gram.md)
  : Gram Matrix by Composite Quadrature
- [`numerical_operator_gram()`](https://statmodels7.github.io/basis7/reference/numerical_operator_gram.md)
  : The Roughness Matrix of an Operator by Quadrature
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
- [`chol_pd()`](https://statmodels7.github.io/basis7/reference/chol_pd.md)
  : Cholesky Factorization, With the Rank Decided Before It
- [`empirical_gram()`](https://statmodels7.github.io/basis7/reference/empirical_gram.md)
  : Gram Matrix Against the Empirical Measure
- [`weighted_gram()`](https://statmodels7.github.io/basis7/reference/weighted_gram.md)
  : Gram Matrix Against a Weighted Lebesgue Measure
- [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  : A Finite-Difference Reference, and Where It Can Be Trusted
- [`tensor_design()`](https://statmodels7.github.io/basis7/reference/tensor_design.md)
  : The Row-Wise Kronecker Product of the Marginal Designs
- [`marginal_designs()`](https://statmodels7.github.io/basis7/reference/marginal_designs.md)
  : Evaluate Every Marginal at Its Own Column
- [`khatri_rao()`](https://statmodels7.github.io/basis7/reference/khatri_rao.md)
  : Row-Wise Kronecker Product of Two Matrices
- [`contract_cp()`](https://statmodels7.github.io/basis7/reference/contract_cp.md)
  : Contract a Tensor Product Basis Against Factor Matrices
- [`clamp_to_range()`](https://statmodels7.github.io/basis7/reference/clamp_to_range.md)
  : Reject Points Outside a Range, and Clamp Those On Its Edge
- [`coord()`](https://statmodels7.github.io/basis7/reference/coord.md)
  [`replace_coord()`](https://statmodels7.github.io/basis7/reference/coord.md)
  : One Coordinate of the Evaluation Points
- [`check_basis_args()`](https://statmodels7.github.io/basis7/reference/check_basis_args.md)
  : Validate the Arguments Every Basis Constructor Takes
- [`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md)
  : Validate Evaluation Points Against a Basis
- [`check_order()`](https://statmodels7.github.io/basis7/reference/check_order.md)
  : Validate a Derivative Order
- [`smoother_interval()`](https://statmodels7.github.io/basis7/reference/smoother_interval.md)
  : The Interval a Smoother Expands Over
- [`check_smoother_x()`](https://statmodels7.github.io/basis7/reference/check_smoother_x.md)
  : Check a Smoother's Covariate
- [`check_whole()`](https://statmodels7.github.io/basis7/reference/check_whole.md)
  : Check a Whole-Number Argument
- [`check_interval()`](https://statmodels7.github.io/basis7/reference/check_interval.md)
  : Check a Smoother's Interval Arguments
- [`check_available()`](https://statmodels7.github.io/basis7/reference/check_available.md)
  : Check What a Smoother Asks For Against What Is Built
- [`check_measure()`](https://statmodels7.github.io/basis7/reference/check_measure.md)
  : Check a Smoother's Measure
- [`check_constrain()`](https://statmodels7.github.io/basis7/reference/check_constrain.md)
  : Check a Polynomial Family's Constraint Argument
- [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
  : The Roughness Matrix a Smoother Penalizes With
- [`poly_free()`](https://statmodels7.github.io/basis7/reference/poly_free.md)
  : The Free Columns of a Polynomial Null Space
- [`poly_free_apply()`](https://statmodels7.github.io/basis7/reference/poly_free_apply.md)
  : The Free Columns of a Polynomial Null Space at New Values
- [`free_names()`](https://statmodels7.github.io/basis7/reference/free_names.md)
  : The Names of the Free Columns
- [`shrink_weight()`](https://statmodels7.github.io/basis7/reference/shrink_weight.md)
  : The Weight a Shrunk Null Space Carries
- [`smoother_reparam()`](https://statmodels7.github.io/basis7/reference/smoother_reparam.md)
  : The Coordinates a Smoother's Coefficients Live In
- [`over_penalty()`](https://statmodels7.github.io/basis7/reference/over_penalty.md)
  : One Operation on a Penalty, However Many Components It Has
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
