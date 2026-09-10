# basis7 0.10.1

* ⚠️ **The cyclic smoother's tests named `penalties7::`, which this package
  may not.** basis7 sits at the bottom of the dependency graph -- it is the
  reason its validation of a penalty factory is deliberately weak, checking
  that the argument is a function of a count and nothing about what the
  function returns -- so naming that package in a test is the one thing the
  arrangement forbids. Both places now use a plain function of `n_coef`,
  and the factory of the "stored and never called" test RAISES when called,
  so the claim is enforced by the construction rather than counted, which
  is the idiom `test-smoother.R` already used.

  ⚠️ **A local suite could not have caught it and did not**: penalties7 is
  installed here, so the tests passed locally at 941 and failed on all five
  CI platforms with `there is no package called 'penalties7'`. Five of five
  is what says it is deterministic rather than platform arithmetic, and the
  failing step was the test run rather than `setup-r`, which is what
  separates it from an infrastructure outage. Nothing about the family
  itself moved: 944 passing, 0 failures, 0 skips.

# basis7 0.10.0

* **`cyclic_smooth()`, the local periodic smoother.** `k` B-spline functions
  over one period, constrained so that the fit and its first `degree - 1`
  derivatives take the same value at the two ends. It is to
  `fourier_smooth()` what `bspline_smooth()` is to `legendre_smooth()`: a
  local basis where the other is global, so a feature at one point of the
  cycle leaves the rest of it alone.

  Nothing outside this package was touched. A family declares a class, a
  `smoother_basis()` method and a constructor, and `modelterms7::s()` reads
  it through `smoother_build()` like any other -- which is the property the
  separation of 0.7.0 was for, tested rather than asserted.

* **The construction is a constraint and not a fold, and that is what makes
  every quantity exact.** A spline of degree \eqn{d} is periodic exactly
  when its value and its first \eqn{d-1} derivatives agree at the two ends,
  which is \eqn{d} linear conditions, so the periodic splines are the null
  space of those conditions inside an ordinary spline space of dimension
  `k + degree`. `constrain_basis()` already built exactly that.

  The parent therefore lives on the period **itself**, so its Gram matrix
  integrates over one period and the roughness matrix is \eqn{T^\top G T}
  with \eqn{G} the parent's own exact Gram: no widened interval, no
  numerical fallback, and evaluation and derivatives of every order exact
  for free. A folded knot sequence puts its parent on a widened interval,
  and its Gram then integrates over more than one period.

* Measured on a cubic with `k = 9` over \eqn{[0, 1]}. The basis agrees at
  the two ends to 5.6e-17 in value, 7.1e-15 in the first derivative and
  8.5e-14 in the second, where an ordinary B-spline of the same `k`
  disagrees by 4.12; a fitted block agrees to 1.1e-16 at every smoothing
  parameter tried. The roughness matrix agrees with a knot-aligned
  trapezoid of the second derivatives over one period with the gap falling
  by exactly **4.00** at each halving of the step -- the reference's own
  order of convergence, so the residual is the trapezoid's and not the
  matrix's -- while integrating over 99 per cent of the period instead
  moves that matrix by 167.6 against a size of 2820.4.

* ⚠️ The null space of the pair is the **constant at every order**, not a
  space growing with it: a non-constant periodic function is never a
  polynomial. Measured at orders 1, 2 and 3 it is one-dimensional in all
  three and its function is constant to 2.7e-15, where the non-periodic
  family of the same degree has a null space of dimension 3 at order 3. So
  `constrain` is not an argument here, for the reason it is not one on
  `fourier_smooth()`, and a block carries `k - 1` columns.

* ⚠️ **The constraints never lose rank, and the guard that would have been
  written for that is not there.** Measured, the rank is exactly `degree`
  at every `k` from 1 to 12 and every `degree` from 1 to 5, so the bound
  `k >= 3` is the sibling periodic family's floor -- the constant is
  removed, so `k = 2` leaves one column -- and nothing else. What the
  construction does assert is the resulting **dimension**, because
  `constrain_basis()` removes one direction per unit of rank and a
  constraint that lost rank would return a basis wider than the caller
  asked for, silently, and the whole block with it.

* ⚠️ **`lower` and `upper` are the period, and they are a property of the
  problem rather than of the sample.** Left `NULL` they are read from the
  data, which makes the fitted period the observed range and is almost
  never what a periodic model means. The page says so.

# basis7 0.9.0

* **A smoother may carry a penalty of its own.** The `penalty` argument of
  `bspline_smooth()`, `fourier_smooth()` and `legendre_smooth()` has been on
  the constructors and documented on their pages since 0.7.0, and a non-NULL
  value was refused; it is now stored. It is a FACTORY and not a built
  penalty, because how many coefficients a smooth has is settled by the data:
  the constraint, the null space and the reparametrization all move with
  them, and a smoother is a recipe for exactly that reason.

* The validation is deliberately weak, and the dependency graph is why.
  \pkg{basis7} sits at the bottom of it and imports \pkg{numericals7} alone,
  so it cannot name \pkg{penalties7} and cannot ask whether what the function
  returns is a penalty. `check_penalty()` asks only that it be a function of
  one argument. The smoother stores it and NEVER CALLS IT, which a test
  asserts with a factory that raises: whichever layer builds the term calls
  it, at the count the data settle, and checks the result there.

* The construction a factory produces is the construction it would have
  produced without one, asserted by identity on the block, the roughness
  matrix and the unpenalized count. The reparametrization reads the
  ROUGHNESS matrix, which is what orders the coordinates from the smoothest
  to the most wiggly and makes the penalty on them the identity; the factory
  replaces that matrix only for whoever penalizes with it. That ordering is
  what a penalty of another shape is reached for.

* ⚠️ `penalty` and `null_space = "shrink"` are refused TOGETHER, and each is
  accepted alone. The shrinkage is a weight written inside the roughness
  matrix -- a tenth of what a penalized direction carries, which is a ratio
  against that matrix's own eigenvalues -- and a factory replaces the matrix
  with a penalty that has no such eigenvalue to be a tenth of. The free
  columns would keep a weight in a matrix the fit no longer reads, which is
  an argument accepted and ignored.

# basis7 0.8.3

* The reapplication test asserts the ROUTE as well as the numbers. A
  tolerance carries the claim that the two agree and cannot carry the claim
  the test is named for -- that the block is reapplied and not rebuilt --
  since a rebuild that happened to land close would pass it. The route is
  asserted structurally instead: with the empirical Gram and the
  reparametrization mocked to raise, `smoother_apply()` completes, so it
  reaches neither. The mock is proven live in the same test by
  `smoother_build()` raising under it, without which the assertion would
  pass whether or not the bindings were ever replaced.

# basis7 0.8.2

* The four tests comparing a reapplied block against the block it came from
  ask a tolerance rather than an identity, through the new
  `expect_reapplied()` test helper. Both routes end in the same product on
  the same basis object, so they perform the same operations; what differs
  is the shape, the build multiplying a matrix with one row per observation
  and the reapplication one with a row per new point, and a BLAS is free to
  block, vectorize and accumulate a product differently by shape. Measured:
  bit-identical under the reference BLAS here and under R release and
  oldrel-1 on the continuous integration's Ubuntu image, and two ulps apart
  under R-devel on that same image, where it separated three of nine numbers
  of one comparison. The tolerance is chosen so the defect the check exists
  for still fails it. That defect is a basis rebuilt from the new points
  rather than reapplied, which carries different knots, a different
  empirical Gram and a different rotation: measured, it sits 3.025 from a
  block whose largest entry is 2.1, twelve orders above the threshold, while
  the platform's last bits are four orders below it. Injection-checked in
  both directions -- two ulps and a thousand ulps pass, a relative
  perturbation of 1e-11, a rebuilt basis and a block of the wrong width are
  each rejected -- with the control reading the threshold from the helper's
  own default rather than from a copy of it. One of the four sites had gone
  red; the other three carry the same shape and were repaired with it.

# basis7 0.8.1

* `smoother_gram()` is exported. It is the roughness matrix a smoother
  penalizes with, at the order and the measure the smoother carries, and a
  tensor product needs it: what `modelterms7::te()` reads from a margin is
  the basis and this matrix, the constraint and the coordinates being the
  product's own.

# basis7 0.8.0

* **`fourier_smooth()`**, the periodic family, and the reason the smoother
  exists. Handed a Fourier basis, the construction `modelterms7::s()` runs
  today loses the property the basis was chosen for: measured on a periodic
  truth at 300 observations, the fitted curve has **`f(0) - f(1) = 2.2125`**,
  and it does so quietly, being neither an error nor a correct fit. Through
  `fourier_smooth(k = 9, lower = 0, upper = 1)` the same data give
  **2.2e-16**, and the fit is better by a factor of seven, rmse 0.0294
  against 0.2038. Every column of the block is itself periodic, so any fit
  built on it is.

* Nine functions give **eight** coordinates rather than seven. The default
  construction removes the constant and the linear function; a periodic
  basis contains no linear function, so removing one costs a degree of
  freedom and buys nothing. `smoother_span()` is the generic where a family
  says what it removes and what it gives back, and the Fourier method
  answers with the constant alone at every order, which is measured: the
  null function of its Gram matrix has a standard deviation of exactly zero
  at orders 1, 2 and 3.

* **`legendre_smooth()`**, the global polynomial family. Its null space is
  the polynomials of degree below `order`, as a B-spline's is, so it needs
  no `smoother_span()` method of its own -- which is the test of whether
  the seam is in the right place.

* **`order`** is built, and it says what a strongly penalized fit contracts
  toward. Measured at `k = 20`, `degree = 5` and a smoothing parameter of
  1e12: at `order = 1` the fitted values have a standard deviation of
  4.8e-11, which is a constant; at 2 they lie on a straight line to an
  R-squared of 1.0000000; at 3 on a parabola to 1.0000000 while a line
  explains only 0.974. The free columns follow the order rather than a
  fixed count -- one at order 2, **two** at order 3, named `lin` and
  `poly2` -- and `order` may not exceed `degree`, above which the roughness
  matrix is identically zero.

* **`constrain`** is built: the directions the smooth is made orthogonal to,
  which may exceed the penalty's null space and may not fall short of it.
  Measured on `y ~ 1 + x + x^2 + s(x)` at `k = 20`, the largest correlation
  between a column of the block and `x^2` is 0.995 at the default and
  **1.9e-15** at `constrain = 2`, which is exact by construction, and the
  standard error of the quadratic coefficient falls with it by more than an
  order of magnitude, at the cost of one dimension. A constraint below the
  null space is rejected where the two numbers were written: a direction
  neither penalized nor identified is an error several frames down
  otherwise.

* **`measure`** is built: `"lebesgue"`, `"empirical"`, or a weight function.
  It is a different penalty and not a detail -- on a cubic B-spline of
  twelve functions at order 2, the correlation between the Lebesgue Gram
  matrix and the one weighted by a Gaussian of standard deviation 0.25 is
  0.11.

* ⚠️ **The default construction is unchanged, bit for bit.** The pipeline
  was rewritten to consume the new hooks, and the identity gate of 0.7.0 --
  56 comparisons against `term_build.SmoothTerm()` over seven shapes -- is
  still 56 identities and no difference. Passing the roughness matrix and
  the constraint to `dr_basis()` explicitly gives what it builds for itself
  at order 2, and the first free column is written as
  `(x - mean(x)) / sd(x)` rather than as the general regression it is a
  case of, because a change of arithmetic there would move fits that are
  not being asked to move.

* **`reparam`** is built, in all three coordinate systems. `"dr"` is the
  Demmler-Reinsch rotation, where the penalty is the identity; `"none"`
  leaves the constrained basis as it stands and the penalty is the
  congruence of the roughness matrix, which is not diagonal; `"orthonorm"`
  rotates so that \eqn{X'X = I} over the observed covariate, against the
  **empirical** measure, `orthonorm_basis()` being the \eqn{L^2} one. The
  three describe the same space, and that is how it is checked: an
  unpenalized fit cannot tell them apart, the fitted values of the three
  agreeing to 1.8e-15 at `k = 12` over 300 observations. It is the
  reparametrized part that is orthonormal -- with a free column prepended
  the whole block is not, and `null_space = "drop"` gives \eqn{X'X = I} for
  the block itself, at 3.1e-15.

* ⚠️ **The block is evaluated from the basis object `smoother_apply()`
  evaluates, so the two are the same arithmetic and not the same formula.**
  Built as a local matrix product instead, the orthonormal route reapplied
  to something that was not identical to the block: `new_transformed()`
  flattens a nested transform, so a later evaluation computes
  \eqn{B(T_1T_2)} where building it in two steps computes \eqn{(BT_1)T_2},
  and the two agree in exact arithmetic and not in the last bit. The
  penalty is read off that same object's transform for the same reason.

* **`null_space = "shrink"`** is built, at a weight of 0.1. It is
  \pkg{mgcv}'s rule translated rather than a number chosen here: measured
  on `mgcv::s(bs = "ts")`, the shrinkage construction leaves the positive
  eigenvalues of the penalty exactly as they were, 7708.76 down to 20.82,
  and replaces each zero with 2.082, a tenth of the smallest positive one.
  In Demmler-Reinsch coordinates every penalized direction has eigenvalue
  exactly 1, so the rule is that single number. ⚠️ Both 0.1 and 1 let the
  term leave the model -- the fitted values reach a standard deviation
  under 1e-6 at a large smoothing parameter, against 0.31 when the null
  space is kept -- and they differ in rate: on a genuinely linear truth at
  a smoothing parameter of 100, the error against that truth is 0.0723 at
  0.1 and 0.4049 at 1, so the heavier weight destroys a real linear trend
  at a smoothing parameter chosen to smooth the wiggles.

* ⚠️ **A penalty factory is the one setting still rejected**, that being a
  lot of its own.

* ⚠️ **The guard that a constraint must leave something to smooth is
  reachable only through `constrain`**, which was found by writing its
  test: with the default, `k` is at least `degree + 1` and `order` is at
  most `degree`, so `k` always exceeds the number of directions removed.
  The guard is kept because a caller can write a large `constrain`, and the
  test reaches it that way.

# basis7 0.7.0

* A **smoother** is a new class: the four decisions a penalized smooth is
  made of, carried as one object. Which functions span the space, what
  counts as roughness, which directions the penalty leaves alone and what
  becomes of them, and which coordinates the coefficients live in. They are
  independent of one another and a basis determines none of them but the
  first, which is why they travel together rather than as separate
  arguments at a call site: the null space is a property of the basis and
  the penalty TOGETHER, so a basis and a penalty chosen separately can be
  an illegal pair. The second-derivative Gram matrix of a cubic B-spline
  has a two-dimensional null space, that of a Fourier basis is
  one-dimensional at every order because the basis contains no linear
  function, and the Gram matrix of a B-spline of degree m at an order above
  m is identically zero, which penalizes nothing.

* `bspline_smooth()` is the first family, and `smoother_build(sm, x)`
  resolves it at a covariate and returns the pair `(X, S)` of
  `solve(X'X + lambda * S, X'y)`, with the count of leading columns the
  penalty does not cover, the names of the coordinates, and a blueprint.
  `smoother_apply(sm, blueprint, newx)` reapplies the recorded construction
  at new values rather than rebuilding it: a rebuild is a basis over
  another interval and a rotation of another Gram matrix, hence a different
  function of the covariate, which the suite pins by measuring that the two
  disagree.

* It is a recipe rather than a built object because two of the four
  decisions need the data: the Demmler-Reinsch rotation diagonalizes the
  pencil of the empirical Gram matrix against the penalty, and the default
  interval is read from the covariate.

* ⚠️ **The construction is the one `modelterms7::s()` has always run, and
  what says so is an identity rather than a tolerance.** Against
  `term_build.SmoothTerm()` on seven shapes -- with and without `by`,
  factor and numeric, sparse and dense, the null space kept and dropped,
  and a non-default `k` and `degree` -- the block, the penalty object, the
  coefficient names and the block at new rows are `identical()`, 56
  comparisons and no difference. The gate was written against a reference
  captured from the installed packages before any of this code existed, and
  it was injection-checked: five defects make it red in the places they
  belong, and the sharpest, the lower endpoint of the interval moved by ONE
  ULP, turns 28 comparisons red. ⚠️ A sixth injection changed nothing and
  is recorded because it looked like a blind spot and is not: adding one
  more `.Machine$double.eps` to the PAD is half an ulp of the endpoint, so
  the subtraction rounds to the same double and the perturbation is not an
  injection at all.

* ⚠️ **A setting this version does not build is rejected rather than
  ignored.** `order`, `measure`, `constrain`, `null_space = "shrink"`,
  `reparam` and a penalty factory are on the constructor because they are
  part of its interface, and each reaches arithmetic that is not written
  yet. An argument accepted and ignored would report a fit of a model the
  caller did not ask for. What is built is `order = 2`, the Lebesgue
  measure, the null space of the basis and the penalty as the constraint,
  `null_space` in `"keep"` and `"drop"`, the Demmler-Reinsch coordinates
  and the quadratic roughness penalty.

* `order` may not exceed `degree`, and the constructor says so where the
  two numbers were written: above the degree the roughness matrix is
  identically zero, so the penalty would leave every direction free.

# basis7 0.6.0

* `plot()` on a basis accepts every argument `matplot()` does. The method
  named `type`, `lty`, `xlab`, `ylab` and `main` in its call, so passing any
  of the five through `...` matched the same formal twice and R threw
  `formal argument "main" matched by multiple actual arguments` before
  anything was drawn. `main` and `xlab` are the two a reader reaches for
  first, and the error named neither the plot method nor where the argument
  came from.

  The five are defaults now, replaced by a value the caller gives, so
  `plot(b, main = "my title")` retitles the panel and `plot(b, type = "p",
  pch = 16)` draws points. Nothing else changes: with none of the five given
  the device output is byte-identical to what it was, which the suite
  asserts.

# basis7 0.5.0

* `basis_numerical_route()` is a new exported generic, and it is what
  `basis_is_numerical()` now asks. The predicate answered by reading which
  class each method is registered on, which says where a method came from
  and not what it does: a method registered on a concrete class that then
  calls the fallback was reported as exact. `basis_gram.FourierBasis()`
  does exactly that whenever the period is not the interval width, so a
  Fourier basis built with `omega = 0.7` reported all three routes exact
  while its Gram matrix came from composite Gauss-Legendre.

* A family whose route depends on its own parameters registers a method and
  is believed over the owner test, which stays the default method. The
  generic is exported because a basis written outside the package has the
  same need; take the default through
  `basis_numerical_route(S7::super(basis, basis7::basis))` and set what
  your own branching decides. `route_by_owner()` is that default under a
  name, for the package's own override, whose formal `basis` shadows the
  class of the same name.

* The correction reaches the two wrapper classes without either of them
  changing: `orthonorm_basis()` of such a Fourier basis reports its Gram
  matrix numerical, since a transformed basis delegates to its parent, and
  so does a `tensor_basis()` carrying one, taking any over its margins.

* `check_basis()` therefore stops comparing that Gram matrix against a
  finer quadrature, which is the treatment its `deriv` and `integral`
  branches already gave a numerical quantity, and holds it to symmetry and
  positive semidefiniteness instead. What was measured is that the two
  agreed to 2.6e-14, both being `numerical_gram()` at different settings.
  `print()` names the route on its `Numerical:` line.

# basis7 0.4.1

* The finite-difference step comes from `numericals7::fd_step()` as the
  offsets and the weights already did, instead of restating its formula.

# basis7 0.4.0

* The finite-difference weights and offsets move to numericals7, where the
  toolkit's one stencil library now lives; this package's Vandermonde
  construction was the most general of the three the toolkit carried and is
  the one that survived. `numerical_deriv_matrix()` keeps its own policy --
  the interval-aware step cap and the switch to one-sided stencils at the
  endpoints -- and speaks to the shared weights for everything else.

# basis7 0.3.1

* The numerical derivative no longer labels the rows of its result after the
  finite-difference stencil each point uses. The stencil is chosen per point
  and recorded in a character vector, whose names traveled through the
  evaluation points and came back as row names for any basis whose method
  propagates the names of `x`.

* Whether a Gram matrix or a penalty is positive definite is now decided from
  its eigenvalues rather than from whether `chol()` raises. On a matrix with an
  exactly zero eigenvalue the pivot that should be zero comes out positive or
  negative according to rounding, so `orthonorm_basis()` and `dr_basis()` gave
  different answers on different machines about the same matrix.

# basis7 0.3.0

* `tensor_basis()`, the product of bases, one per variable. Everything follows
  from the marginals because the product separates: a partial derivative
  differentiates one marginal and leaves the others, the integral over the box
  from its lower corner is the product of the marginal integrals, and the Gram
  matrix is the Kronecker product of the marginal ones. A tensor of exactly
  integrated marginals is therefore exact at any number of variables, where a
  quadrature over the box would not be.

* A basis now declares how many variables it takes, through `basis_nvar()` and
  the length of its endpoints, so a univariate basis is the case of one
  variable rather than a separate kind of object. Evaluation points become a
  matrix with one column per variable, and the derivative order becomes a
  multi-index. A single non-zero order is rejected for a product, since it
  could mean that order in every coordinate or that total order.

* `basis_contract()` evaluates a basis against coefficients without
  necessarily forming the design matrix. Coefficients come as a full array,
  processed in row blocks so the peak memory is bounded by the block rather
  than by the sample, or as a list of factor matrices in canonical polyadic
  form, where the cost is linear in the number of variables instead of
  exponential in it. A product of six bases of ten functions has a million
  columns; the factorized route touches neither that matrix nor the
  coefficient array.

* `check_basis()`, `print()` and the numerical fallbacks understand several
  variables. `plot()` rejects a product with an explanatory message.

# basis7 0.2.0

* `poly_basis()`, the Legendre polynomials by recurrence. They span the same
  space as the raw powers and are chosen over them for conditioning: the Gram
  matrix of raw powers is a Hilbert matrix, and ten of them are already close
  to singular in double precision.

* `basis_gram()` takes a measure. The default is Lebesgue on the interval;
  `at` gives the empirical measure of a sample, which is the matrix a design
  matrix produces, and `weight` gives a weighted integral. Both are handled in
  the body of the generic, before dispatch, so a method never implements them.
  **Breaking for user-written methods**: a `basis_gram()` method must now name
  `at` and `weight` in its signature, because S7 requires a method's formals to
  contain the generic's.

* `TransformedBasis`, one class for every linear reparametrization
  \eqn{B \mapsto BT}, with three constructors: `orthonorm_basis()`, from the
  Cholesky factor of the Gram matrix rather than from a grid;
  `constrain_basis()`, from the null space of a constraint; and `dr_basis()`,
  the Demmler-Reinsch construction, which diagonalizes the empirical inner
  product and the penalty at once and is empirically orthogonal to a constant
  and to the covariate.

  Transforms compose by multiplication rather than by nesting, and a
  transformed basis reports the *parent's* numerical status, since its own
  methods delegate and multiply.

* `dr_basis()` factorizes the penalty and not the design, so it survives a
  rank-deficient design, which equally spaced knots produce whenever the data
  leave a knot span empty. Verified on a design where the Cholesky factor of
  the design matrix does not exist.

* `check_basis()` allows for the accuracy of its own reference. Each numerical
  reference is computed at a step and at half of it, and the gap between them
  bounds its uncertainty; the comparison is given that much slack, point by
  point. Without it a spline failed at its own knots, where the third
  derivative jumps and a central difference returns the jump rather than the
  truncation error. A deliberate five per cent error is still caught by four
  orders of magnitude.

* A vignette, `defining-a-basis`.

# basis7 0.1.0

First release.

* The `basis` class and four generics: `basis_eval()`, `basis_deriv()`,
  `basis_int()` and `basis_gram()`. Derivative order is an argument, so no
  order is privileged, and `basis_gram()` returns the inner products a
  roughness penalty integrates.

* Two families with exact formulas: `bspline_basis()`, whose Gram matrix is
  integrated exactly knot interval by knot interval, and `fourier_basis()`,
  whose derivatives of every order and whose antiderivative come from one
  phase-shift identity, and whose Gram matrix is diagonal in closed form over
  a whole period.

* Numerical fallbacks registered on the base class, so a basis that implements
  only `basis_eval()` is complete. Derivatives use one finite-difference
  stencil built from a Vandermonde solve, never a chain of first differences,
  and switch to a one-sided stencil at the interval endpoints.

* `basis_int()` is anchored: its value at the lower endpoint is exactly zero,
  for every family. Any antiderivative satisfies a differentiation check, so
  the constant is fixed by the contract instead.

* `check_basis()` runs six numerical checks and reports a quantity that came
  from a fallback as unchecked rather than as passed, and a property the
  family never claimed as not claimed.

* `basis_is_numerical()` answers the same question programmatically.

* `print()` and `plot()` methods; the plot draws the basis, any derivative, or
  the integral.
