# Changelog

## basis7 0.15.0

### Repairs

- [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
  reads the roots of the characteristic polynomial on the scale of the
  roots, so its tolerance is relative to their size. With an absolute
  tolerance a long period merged the roots of the harmonic operator:
  from a period of 1e7 the null space read `1, t, t^2` instead of the
  constant, the sine and the cosine, and a short period put a spurious
  rate in the labels (`exp(5.11e-09 t) sin(6283.19 t)` at a period of
  1e-3). The null space is now correct for periods from 1e-3 to 1e9.
- The numerical derivative uses one-sided stencils of `order + 2` nodes
  at the ends of the interval, so they keep the second-order accuracy of
  the central stencil. At even orders the one-sided stencil had one node
  fewer and was first-order accurate: the second derivative of a
  Legendre basis of dimension 5 had a relative error of 5.7e-04 at the
  endpoints, against 1.5e-08 inside; it is now 1.3e-07.
- [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
  removes the constant also when the operator does not have it in its
  null space, as for
  [`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md).
  The block kept the constant as a penalized direction, and beside an
  intercept it had rank 9 of 10. The block now has one column fewer, as
  with the harmonic operator.
- [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
  checks at construction that the constant and the null space of the
  operator leave at least one function to penalize, and names the
  smallest admissible `k`. The error previously came from
  [`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
  and did not name the cause.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a basis
  draws derivatives of order 2 and above; it signalled
  `object 'B' not found`.
- [`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
  accepts an array or a list of factor matrices for a tensor basis built
  with named margins, and a named list of factors. The names made
  [`identical()`](https://rdrr.io/r/base/identical.html) reject
  dimensions that were correct.
- [`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
  of a tensor basis at a matrix of no rows returns `numeric(0)`, and
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  and
  [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  of a B-spline basis at no points return a matrix of no rows, where
  signalled an error.
- [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  accepts an order per variable for a basis of several variables.
- [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  and
  [`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
  signal an error for an operator whose period has not been resolved;
  the matrix was returned with `NA` and `NaN` on the diagonal.
- [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  and
  [`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md)
  signal an error when `panels` or `nodes` reaches a route that computes
  the matrix exactly (a B-spline or Legendre basis, or a Fourier basis
  over a whole period, without `weight`); the two arguments were
  accepted and ignored.
- A basis of several variables evaluated at a plain vector whose length
  is not a multiple of the number of variables signals an error. The
  vector was recycled with a warning, and the last point was made up of
  recycled coordinates.
- A smoother built at a covariate that takes a single value, with
  neither endpoint of the interval given, signals an error that names
  the cause; the error came from .
- The validator of `TensorBasis` checks that every margin takes one
  variable and that `@lower` and `@upper` are the endpoints of the
  margins.
- [`print()`](https://rdrr.io/r/base/print.html) of a smoother built
  with a `penalty` factory states the factory, and names the roughness
  matrix that sets the coordinates `roughness`; it printed the roughness
  matrix as the penalty.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  reports a missing value in an analytic derivative as a failure. The
  comparison returned `NA`, and the table printed `[numerical]` for a
  derivative that has a method. An entry whose reference is missing, and
  whose allowance is therefore infinite, is left out of the comparison,
  as the documentation of
  [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  states.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  compares a numerical Gram matrix with the finer quadrature of 401
  panels of 7 nodes, allowing each entry four times the difference
  between the package’s own fallback and that rule. The comparison was
  skipped, so a Gram method wrong by one percent on a Fourier basis off
  a whole period passed; a correct fallback on a basis with kinks, whose
  quadrature error is 6.6e-4, still passes.
- [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
  of a Fourier basis, with the method called at order 0, gives the
  constant column 1; it gave 0. The generic returns the evaluation at
  order 0 and was not affected.
- [`print()`](https://rdrr.io/r/base/print.html) of a smoother with one
  endpoint of the interval given shows it, as in `[-5, from the data]`;
  it printed `from the data`.
- The message for an order of the wrong length on a basis of one
  variable no longer describes a basis of several variables.
- The internal
  [`gauss_legendre()`](https://statmodels7.github.io/basis7/reference/gauss_legendre.md)
  rejects a number of nodes that is not a positive whole number;
  `nodes = 2.7` was truncated to a two-point rule.
- The internal
  [`chol_pd()`](https://statmodels7.github.io/basis7/reference/chol_pd.md)
  returns `NULL` for an empty matrix or one with a value that is not
  finite, as documented; [`eigen()`](https://rdrr.io/r/base/eigen.html)
  signalled an error first.
- The messages for a `k` too small for the degree state the numbers, and
  the messages that list the operator constructors include
  [`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md).

### Documentation

- The documentation, the README, the vignette and this file are
  rewritten in a plain register, and the statements that did not match
  the code are corrected. Among them: the null space of a derivative
  penalty on a periodic basis is the constant, so a strongly penalized
  Fourier or cyclic smooth contracts to a constant and not to a straight
  line; the one-harmonic operator has a penalty of null dimension 3 on a
  Fourier basis and 1 on a cubic B-spline; the composition of two
  operators has a null space that contains, and can exceed, the union of
  the two; the nodes of
  [`quad_rule()`](https://statmodels7.github.io/basis7/reference/quad_rule.md)
  are ordered by position within the rule and then by interval;
  [`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md)
  takes four times the gap between its two estimates as the bound on its
  error; and the Hilbert matrix of dimension 15 has a condition number
  of 6.1e+20.

## basis7 0.14.0

- [`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
  takes its differences on the coefficients of the basis of Eilers and
  Marx (1996), whose knots are equally spaced beyond the interval as
  well as inside it, and carries the penalty back to the clamped basis
  through the new internal
  [`pspline_map()`](https://statmodels7.github.io/basis7/reference/pspline_map.md).
  The two bases span the same splines, so only the coordinates of the
  penalty change. On the clamped coefficients the null space was only
  approximately the polynomials of degree below `diff` (an of 0.9994 at
  `diff = 2`); it is now exact. The penalty is that of the `bs = "ps"`
  smooths of mgcv: on
  [`MASS::Boston`](https://rdrr.io/pkg/MASS/man/Boston.html),
  `medv ~ s(lstat)` has 6.98906, 8.75744 and 9.41421 effective degrees
  of freedom at `k` of 10, 20 and 40, against 6.98927, 8.75783 and
  9.41421 in mgcv; before, it had 8.29 at `k = 20` against 8.76.
- [`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md)
  takes its differences in the same coordinates, so at equal smoothing
  parameters it gives the penalty of `pspline_smooth(reparam = "none")`,
  and its fit is close to the adaptive P-spline of mgcv’s `bs = "ad"`:
  on [`MASS::mcycle`](https://rdrr.io/pkg/MASS/man/mcycle.html) at
  `k = 40`, `m = 5`, 10.3352 effective degrees of freedom against
  10.3343.
- Every fit with a P-spline or an adaptive smooth changes. `splines` is
  added to Imports for
  [`splines::splineDesign()`](https://rdrr.io/r/splines/splineDesign.html).

## basis7 0.13.1

- [`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md)
  builds with `reparam = "none"` and `reparam = "orthonorm"`, which
  signalled “non-conformable arguments” at every `k` since the family
  was added in 0.10.0. A cyclic basis is itself a transformed basis, and
  [`new_transformed()`](https://statmodels7.github.io/basis7/reference/new_transformed.md)
  flattens a nested transform, so the constrained object carried a
  transform against the parent B-spline basis (15 by 11 instead of 12 by
  11). The congruence of the roughness matrix now uses the local
  transform, computed by the new internal
  [`constraint_null()`](https://statmodels7.github.io/basis7/reference/constraint_null.md),
  which
  [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  also calls. Other families are unaffected: over seventeen
  configurations the results are
  [`identical()`](https://rdrr.io/r/base/identical.html) to those of the
  previous release, except the two that signalled the error.
- The tests build every family in every coordinate system, reapply the
  block at new values, and check that the two periodic families remain
  periodic.

## basis7 0.13.0

- Linear differential operators, , with the roughness matrix . The
  constructors are
  [`deriv_operator()`](https://statmodels7.github.io/basis7/reference/deriv_operator.md),
  [`harmonic_operator()`](https://statmodels7.github.io/basis7/reference/harmonic_operator.md),
  [`oscillator_operator()`](https://statmodels7.github.io/basis7/reference/oscillator_operator.md)
  and
  [`linear_operator()`](https://statmodels7.github.io/basis7/reference/linear_operator.md);
  `*` composes two by multiplying their characteristic polynomials; and
  [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md),
  [`operator_null_design()`](https://statmodels7.github.io/basis7/reference/operator_null_design.md),
  [`operator_order()`](https://statmodels7.github.io/basis7/reference/operator_order.md),
  [`operator_weights()`](https://statmodels7.github.io/basis7/reference/operator_weights.md)
  and
  [`operator_resolve()`](https://statmodels7.github.io/basis7/reference/operator_resolve.md)
  are the accessors.
- An operator is given through the `order` argument of every smoother
  family, a whole number `m` being the shorthand for
  `deriv_operator(m)`. `basis_gram(b, order = op)` computes the matrix
  for a basis alone, through the new generic
  [`basis_operator_gram()`](https://statmodels7.github.io/basis7/reference/basis_operator_gram.md).
- The null space of an operator is read from the roots of its
  characteristic polynomial and not from the rank of the assembled
  matrix. The two differ: the one-harmonic operator has a penalty of
  null dimension 3 on a Fourier basis and 1 on a cubic B-spline at a
  relative tolerance of 1e-10, because a spline represents a sine only
  approximately, while the null space of the operator is
  three-dimensional in both cases.
- [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
  penalizes with the harmonic acceleration operator and shrinks its null
  space by default, where it penalized the squared second derivative and
  kept the null space. On a periodic basis a derivative penalty has the
  constant as its only null direction, so a strongly penalized fit
  contracts to a constant; the harmonic operator leaves the level and
  the fundamental cycle unpenalized, so the fit contracts to a sinusoid.
  Every Fourier smooth changes; `order = 2` restores the previous
  construction exactly.
- At `k = 21`, over 200 observations and eight samples, with the
  smoothing parameter chosen by generalized cross-validation, the
  harmonic penalty had a root mean square error of 0.0596 against 0.0654
  on a truth dominated by its fundamental. On a truth with no periodic
  signal the default `"shrink"` gave 0.0110 against 0.0121 for the
  derivative penalty, and `"keep"` gave 0.0281, which is why the default
  null space of this family is `"shrink"`. It resolves to `"keep"` when
  a `penalty` factory is given, the two being incompatible.
- With `null_space = "keep"` the block of a Fourier smooth changes as
  well: the constraint removes the constant in both constructions, and
  the harmonic operator frees the sine and cosine of the fundamental, so
  nine functions still give eight columns, named `sin1`, `cos1`, `z1` to
  `z6` instead of `z1` to `z8`.
- A periodic family rejects an operator whose null space holds a
  function that is not periodic on the period of the basis, such as in
  `deriv_operator(2) * oscillator_operator(p)`. A plain derivative
  operator restores nothing beyond the constant and is accepted.
- On a Fourier basis over a whole period the roughness matrix of any
  operator is diagonal, with the entry of the pair at frequency equal to
  for the characteristic polynomial.
- On a B-spline basis the roughness matrix of an operator is exact,
  integrated knot interval by knot interval as the derivative one is.
  The equally spaced panels of the base method do not line up with the
  knots and are approximate. An operator of order above the degree of
  the spline is rejected, its leading term being zero there.
- `smoother@order` holds a `LinearOperator` and no longer an integer; a
  whole number is converted by
  [`as_operator()`](https://statmodels7.github.io/basis7/reference/as_operator.md)
  at the constructor, and `operator_order(sm@order)` is the integer.
  Against results recorded before the change, fourteen smoother shapes
  and eight Gram matrices are
  [`identical()`](https://rdrr.io/r/base/identical.html), and the three
  that change are the Fourier default.

## basis7 0.12.0

- [`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md),
  a difference penalty whose weight varies along the covariate. Each
  difference carries a weight , with a B-spline basis over the
  coefficient index, so the penalty is the sum and its smoothing
  parameters are the coefficients of the weight profile.
- The weight functions are a partition of unity, so the components sum
  to and the family at equal smoothing parameters is a P-spline (to
  between 8.9e-16 and 2.7e-15 for `m` from 2 to 12).
- Against a P-spline with one smoothing parameter fitted by mgcv’s REML
  (400 observations, `k = 40`, eight seeds), the median root mean square
  error was 0.0357 against 0.0443 on a truth of variable roughness, at
  17.4 effective degrees of freedom against 23.6, and 0.0315 against
  0.0310 on a truth of constant roughness. Split by region, the gain
  came from the quiet stretches of the function and not from the
  feature.
- The index basis is built over the range of the index, so an affine
  relabeling of the index gives the same profile.
- `reparam = "dr"` is rejected, the Demmler-Reinsch rotation
  diagonalizing one penalty and this family having `m`; the default is
  `"none"`. `null_space = "shrink"` and a `penalty` factory are rejected
  as well.
- The weight profile differs from mgcv’s, whose P-spline basis extends
  its knots beyond the index range; the assembly of the components
  follows mgcv’s constructor and agrees with
  `smoothCon(scale.penalty = FALSE)`.
- A smoother’s penalty may be a list of matrices, which
  [`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
  carries through.
  [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
  may return one matrix or several,
  [`over_penalty()`](https://statmodels7.github.io/basis7/reference/over_penalty.md)
  applies each step to one or to each, and
  [`smoother_reparam()`](https://statmodels7.github.io/basis7/reference/smoother_reparam.md)
  rejects a list under `"dr"`. The number of unpenalized leading columns
  is read from the sum of the components.
- [`print()`](https://rdrr.io/r/base/print.html) of a P-spline smoother
  states `difference of order 2 on the coefficients`, where it printed a
  derivative and a measure.

## basis7 0.11.0

- [`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md),
  the Eilers-Marx smoother: a B-spline basis of `k` functions penalized
  by the sum of squared `diff`-th differences of its coefficients. It
  has no `measure` and no `order`, since nothing is integrated, and
  `diff` belongs to this family alone, because only the coefficients of
  a B-spline form an ordered sequence in which neighbors are comparable.
- [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
  is a generic, so a family may declare a roughness matrix that is not a
  Gram matrix of its derivatives. Its base method is the previous
  function, and over 180 builds of the shipped families the results are
  [`identical()`](https://rdrr.io/r/base/identical.html) before and
  after.
- On a clamped knot sequence the null space of a difference penalty is
  only approximately the polynomials, the Greville abscissae not being
  equally spaced near the ends. The block is constrained against the
  exact polynomials, so the Demmler-Reinsch rotation runs on their
  complement and a strongly penalized fit contracts to a straight line,
  as for
  [`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md).
- The raw P-spline and integrated-derivative penalties differ in scale
  by orders of magnitude, the difference operator carrying no factor of
  the knot spacing; after the Demmler-Reinsch rotation both are the
  identity, so their smoothing parameters are comparable.

## basis7 0.10.1

- The tests of the cyclic smoother no longer name , which this package
  does not depend on; they use a plain function of `n_coef`, and the
  factory of the test for a factory that is stored and never called
  signals an error when called. The tests had passed locally, where
  penalties7 is installed, and failed on the CI platforms.

## basis7 0.10.0

- [`cyclic_smooth()`](https://statmodels7.github.io/basis7/reference/cyclic_smooth.md),
  the local periodic smoother: `k` B-spline functions over one period,
  constrained so that the fit and its first `degree - 1` derivatives
  take the same value at the two ends. It is built with
  [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  on a B-spline basis of dimension `k + degree` over the period itself,
  so the Gram matrix integrates over one period and every quantity is
  exact. On a cubic with `k = 9` over the basis functions agree at the
  two ends to 5.6e-17 in value, 7.1e-15 in the first derivative and
  8.5e-14 in the second.
- The null space of the roughness matrix of a derivative penalty is the
  constant at every order, a non-constant periodic function being never
  a polynomial, so `constrain` is not an argument and a block carries
  `k - 1` columns.
- The periodicity constraints have rank `degree` (checked at every `k`
  from 1 to 12 and every `degree` from 1 to 5), and the construction
  checks the resulting dimension.
- `lower` and `upper` are the ends of the period; left `NULL` they are
  read from the data, which makes the period the observed range.
- A family is added with a class, a
  [`smoother_basis()`](https://statmodels7.github.io/basis7/reference/smoother_basis.md)
  method and a constructor; `modelterms7::s()` reads the cyclic family
  through
  [`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
  with no change outside this package.

## basis7 0.9.0

- The `penalty` argument of
  [`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md),
  [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md)
  and
  [`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md)
  stores a factory of the coefficient count; a non-NULL value was
  rejected before. The smoother never calls it: the layer that builds
  the term calls it at the count settled by the data. Of the toolkit
  this package imports numericals7 only, so
  [`check_penalty()`](https://statmodels7.github.io/basis7/reference/check_penalty.md)
  checks only that the factory is a function with at least one argument.
- A factory does not change the construction: the block, the roughness
  matrix and the unpenalized count are those built without it, the
  reparametrization reading the roughness matrix.
- `penalty` and `null_space = "shrink"` are rejected together and
  accepted separately, the shrinkage being a weight inside the roughness
  matrix that a factory replaces.

## basis7 0.8.3

- The test of reapplication checks that
  [`smoother_apply()`](https://statmodels7.github.io/basis7/reference/smoother_apply.md)
  reapplies the block and does not rebuild it, by mocking the empirical
  Gram matrix and the reparametrization to signal an error.

## basis7 0.8.2

- The four tests that compare a reapplied block with the block it came
  from use a tolerance, through the new test helper
  `expect_reapplied()`. The two routes compute the same product on
  matrices of different shapes, which a BLAS may accumulate differently:
  the results were identical under the reference BLAS and two ulps apart
  under R-devel on the Ubuntu image of the
  101. A rebuilt basis differs by 3.025 on a block whose largest entry
       is 2.1, and still fails the comparison.

## basis7 0.8.1

- [`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
  is exported, for `modelterms7::te()`, which reads the basis and the
  roughness matrix of each margin.

## basis7 0.8.0

- [`fourier_smooth()`](https://statmodels7.github.io/basis7/reference/fourier_smooth.md),
  the periodic smoother. Every column of its block is periodic, so a fit
  built on it is; the construction of `modelterms7::s()` applied to a
  Fourier basis gave a fit with `f(0) - f(1) = 2.2125` on a periodic
  truth.
- A Fourier smooth of nine functions has eight coordinates: the default
  construction removes the constant and the linear function, and a
  periodic basis contains no linear function.
  [`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
  is the generic where a family states what it removes and what it
  restores.
- [`legendre_smooth()`](https://statmodels7.github.io/basis7/reference/legendre_smooth.md),
  the global polynomial smoother, whose null space is the polynomials of
  degree below `order`, as for a B-spline.
- `order` sets the derivative of the penalty and so the function toward
  which a strongly penalized fit contracts; the free columns follow it,
  one at `order = 2` and two at `order = 3`, named `lin` and `poly2`.
  `order` may not exceed `degree`.
- `constrain` sets the directions to which the smooth is made
  orthogonal, which may exceed the null space of the penalty and may not
  fall short of it.
- `measure` is `"lebesgue"`, `"empirical"`, or a weight function.
- `reparam` is `"dr"`, `"none"` or `"orthonorm"`, the last orthonormal
  against the empirical measure.
- `null_space = "shrink"`, at a weight of 0.1, the shrinkage
  construction of mgcv’s `bs = "ts"` smooths in Demmler-Reinsch
  coordinates.
- The default construction is unchanged: the 56 comparisons of 0.7.0 are
  still identical.
- A penalty factory is still rejected.

## basis7 0.7.0

- The smoother class carries the four choices of a penalized smooth: the
  basis, the penalty, the treatment of the null space, and the
  coordinates of the coefficients. The null space is a property of the
  basis and the penalty together, so they are carried as one object.
- [`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
  is the first family. `smoother_build(sm, x)` returns the block `X`,
  the penalty `S`, the number of unpenalized leading columns, the names
  of the coordinates and a blueprint, and
  `smoother_apply(sm, blueprint, newx)` reapplies the recorded
  construction at new values.
- The construction is that of `modelterms7::s()`: against
  `term_build.SmoothTerm()` on seven shapes, 56 comparisons of the
  block, the penalty, the names and the block at new rows are
  [`identical()`](https://rdrr.io/r/base/identical.html).
- Settings that this version did not build (`order`, `measure`,
  `constrain`, `null_space = "shrink"`, `reparam` and a penalty factory)
  are rejected.
- `order` may not exceed `degree`, the roughness matrix being zero above
  it.

## basis7 0.6.0

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a basis
  accepts every argument of
  [`matplot()`](https://rdrr.io/r/graphics/matplot.html). Passing
  `type`, `lty`, `xlab`, `ylab` or `main` signalled
  `formal argument "main" matched by multiple actual arguments`; the
  five are now defaults that a given value replaces.

## basis7 0.5.0

- [`basis_numerical_route()`](https://statmodels7.github.io/basis7/reference/basis_numerical_route.md),
  an exported generic called by
  [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md).
  The predicate read the class on which each method is registered, so a
  Fourier basis whose period is not the width of its interval reported
  an exact Gram matrix computed by quadrature. A family whose route
  depends on its parameters registers a method, whose result replaces
  the owner test;
  [`route_by_owner()`](https://statmodels7.github.io/basis7/reference/route_by_owner.md)
  is the default.
- [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md)
  and
  [`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md)
  of such a Fourier basis report the Gram matrix as numerical.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  no longer compares a numerical Gram matrix with a second quadrature,
  and tests its symmetry and positive semidefiniteness.
  [`print()`](https://rdrr.io/r/base/print.html) names the route on its
  `Numerical:` line.

## basis7 0.4.1

- The finite-difference step comes from
  [`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html).

## basis7 0.4.0

- The finite-difference weights and offsets come from numericals7.
  [`numerical_deriv_matrix()`](https://statmodels7.github.io/basis7/reference/numerical_deriv_matrix.md)
  keeps the cap of the step by the interval and the one-sided stencils
  at the endpoints.

## basis7 0.3.1

- The numerical derivative no longer names the rows of its result after
  the stencil used at each point.
- Whether a Gram matrix or a penalty is positive definite is decided
  from its eigenvalues and not from whether
  [`chol()`](https://rdrr.io/r/base/chol.html) succeeds, which on a
  singular matrix depends on rounding and gave different results on
  different platforms.

## basis7 0.3.0

- [`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
  the product of bases, one per variable. A partial derivative
  differentiates one margin, the integral over the box is the product of
  the marginal integrals, and the Gram matrix is the Kronecker product
  of the marginal ones, so a product of exact margins is exact at any
  number of variables.
- [`basis_nvar()`](https://statmodels7.github.io/basis7/reference/basis_nvar.md)
  gives the number of variables of a basis; evaluation points become a
  matrix with one column per variable, and the derivative order a
  multi-index. A single non-zero order is rejected for a product.
- [`basis_contract()`](https://statmodels7.github.io/basis7/reference/basis_contract.md)
  evaluates a basis against coefficients without forming the full design
  matrix: in row blocks for an array of coefficients, and from the
  margins alone for factor matrices in canonical polyadic form, at a
  cost linear in the number of variables.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md),
  [`print()`](https://rdrr.io/r/base/print.html) and the numerical
  fallbacks handle several variables;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) rejects a
  product.

## basis7 0.2.0

- [`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md),
  the Legendre polynomials by recurrence, better conditioned than the
  raw powers, whose Gram matrix is a Hilbert matrix.
- [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  takes a measure: Lebesgue on the interval by default, the empirical
  measure of a sample with `at`, a weighted integral with `weight`.
  **Breaking for user-written methods**: a
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  method must name `at` and `weight` in its signature, S7 requiring a
  method’s formals to contain the generic’s.
- `TransformedBasis`, one class for the linear reparametrizations , with
  [`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
  [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  and
  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md),
  the Demmler-Reinsch construction. Transforms compose by
  multiplication, and a transformed basis reports the numerical status
  of its parent.
- [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  factorizes the penalty and not the design, so it works on a
  rank-deficient design, which equally spaced knots produce when the
  data leave `degree + 1` or more consecutive knot spans empty.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  allows for the error of its reference: each numerical reference is
  computed at a step and at half of it, and four times the gap between
  them is allowed point by point. A spline no longer fails at its knots,
  and a derivative wrong by five percent still fails.
- A vignette, `defining-a-basis`.

## basis7 0.1.0

First release.

- The `basis` class and four generics:
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
  [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md).
  The derivative order is an argument, and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  returns the inner products a roughness penalty integrates.
- Two families with exact formulas:
  [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
  whose Gram matrix is integrated exactly knot interval by knot
  interval, and
  [`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
  whose derivatives and antiderivative come from one phase-shift
  identity and whose Gram matrix is diagonal in closed form over a whole
  period.
- Numerical methods on the base class, so a basis that implements only
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
  is complete. Derivatives use one finite-difference stencil, never a
  chain of first differences, with a one-sided stencil at the endpoints.
- [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  is anchored at the lower endpoint, where its value is zero for every
  family.
- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  runs six numerical checks, marking a quantity that comes from a
  fallback as `[numerical]` and a property that the family does not have
  as `[not claimed]`.
- [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
  reports which quantities are numerical.
- [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods; the
  plot draws the basis, a derivative or the integral.
