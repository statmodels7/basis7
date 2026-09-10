#' @include basis_class.R
NULL


#' Evaluate a Basis
#'
#' @description
#' Returns the design matrix of a basis at the given points: one row per
#' evaluation point, one column per basis function, entry \eqn{(i, j)} equal
#' to \eqn{\varphi_j(x_i)}. This is the matrix a regression on the basis is
#' fitted against, and the one every other generic in the package is defined
#' in terms of.
#'
#' @details
#' # The design matrix
#'
#' An expansion with coefficients \eqn{\beta} is evaluated as
#'
#' \deqn{f(x) = \sum_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\,\beta,}
#'
#' so `basis_eval(b, x) %*% beta` is the fitted function at `x`, and the
#' matrix is the design block a linear model on the basis uses. Its columns
#' carry the names [basis_colnames()] declares, and every other matrix the
#' basis produces carries the same ones in the same order.
#'
#' # The one generic a basis must implement
#'
#' [basis_deriv()], [basis_int()] and [basis_gram()] all have numerical
#' methods registered on the abstract [basis] class, computed from this one,
#' so a subclass supplying its evaluation alone answers all four. Registering
#' a closed form for any of the three later takes over through dispatch;
#' [basis_is_numerical()] reports which are still on the fallback.
#'
#' # What the generic does before dispatching
#'
#' The generic body validates `x` through [check_eval_points()] and passes
#' the validated version on, so every method, including one written outside
#' the package, gets the same guarantees without writing them: a point outside
#' the interval throws, a point that is an endpoint up to a relative `1e-8`
#' arrives clamped exactly onto that endpoint, and a basis of several
#' variables receives a matrix of [basis_nvar()] columns whatever shape the
#' caller passed.
#'
#' `NA` is neither checked nor clamped and flows through arithmetic, so a
#' missing evaluation point gives a row of `NA`.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points inside the basis interval: a numeric vector for
#'   a basis of one variable, or a matrix of [basis_nvar()] columns for a
#'   basis of several, where a plain vector is taken by row. `NA` is allowed
#'   and produces a missing row. A point outside the interval throws.
#' @param ... Passed to methods. No shipped family reads anything from it.
#'
#' @return A numeric matrix of `n` rows and `basis@dimension` columns, `n`
#'   being `length(x)` for one variable and `nrow(x)` for several, with column
#'   names [basis_colnames()].
#'
#' @seealso [basis_deriv()] for its derivatives, [basis_int()] for its
#'   anchored integral, [basis_gram()] for its inner products, and
#'   [basis_contract()] to get \eqn{B(x)\beta} without forming \eqn{B(x)}
#'   for a product basis.
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' round(basis_eval(b, c(0, 0.5, 1)), 4)
#'
#' # A B-spline basis is a partition of unity: every row sums to one.
#' rowSums(basis_eval(b, seq(0, 1, length.out = 5)))
#'
#' # Evaluating an expansion is one matrix product.
#' beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
#' drop(basis_eval(b, c(0.3, 0.7)) %*% beta)
#'
#' # Outside the interval it throws; an endpoint up to rounding is clamped.
#' try(basis_eval(b, 1.5))
#' all.equal(basis_eval(b, 1 + 5e-9), basis_eval(b, 1))
#'
#' # A basis of several variables takes one column per variable.
#' tb <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#' dim(basis_eval(tb, cbind(c(0.1, 0.5), c(0.2, 0.6))))
#'
#' @export
basis_eval <- S7::new_generic("basis_eval", "basis", function(basis, x, ...) {
  x <- check_eval_points(basis, x)
  S7::S7_dispatch()
})


#' Differentiate a Basis
#'
#' @description
#' Returns the `order`-th derivative of every basis function at the given
#' points, as a matrix of the same shape [basis_eval()] returns. Since an
#' expansion is linear in its coefficients, this matrix times \eqn{\beta} is
#' the `order`-th derivative of the fitted function, so a derivative of a fit
#' needs no refitting.
#'
#' @details
#' # Order is an argument
#'
#' Derivative order is an argument, not a family of generics, because it is
#' unbounded: a Fourier basis is differentiable to any order, and a spline of
#' degree \eqn{k} has \eqn{k} non-trivial derivatives and zeros above that.
#' An order beyond what the family carries returns the zero matrix, which is
#' the value of the derivative; nothing is thrown, and the zeros are exact.
#'
#' `order = 0` short-circuits to [basis_eval()] in the generic body, so a loop
#' over orders needs no special case and pays nothing for the zero.
#'
#' # A basis of several variables takes a multi-index
#'
#' For a product basis `order` has one entry per variable and names a mixed
#' partial: `c(2, 0)` is \eqn{\partial^2/\partial x_1^2} and `c(1, 1)` is
#' \eqn{\partial^2/\partial x_1 \partial x_2}. A single non-zero number is
#' refused, having two readings; `0` alone is accepted, meaning no derivative
#' under either. See [check_order()].
#'
#' # The numerical fallback
#'
#' A subclass registering no method gets the one on the abstract [basis]
#' class, which applies **one** stencil of the order asked for to
#' [basis_eval()], never a composition of lower-order differences. The
#' offsets, weights and step come from [numericals7::fd_offsets()],
#' [numericals7::fd_weights()] and [numericals7::fd_step()], and the stencil
#' is shifted to one side near an endpoint so that no node leaves the
#' interval. [basis_is_numerical()] says whether this is the route in use, and
#' [check_basis()] measures the agreement.
#'
#' For a basis of several variables the fallback differentiates one coordinate
#' at a time, so a mixed partial such as `c(1, 1)` throws there: a stencil in
#' the plane has the product of two errors, and the one family that needs
#' mixed partials, [tensor_basis()], computes them exactly from its margins.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points inside the basis interval: a numeric vector for
#'   a basis of one variable, or a matrix of [basis_nvar()] columns for a
#'   basis of several.
#' @param order The derivative order. A single non-negative whole number for a
#'   basis of one variable, default `1`; a vector of [basis_nvar()] such
#'   numbers for a basis of several, or a single `0`. A negative, fractional
#'   or missing order throws.
#' @param ... Passed to methods. No shipped family and no fallback reads
#'   anything from it.
#'
#' @return A numeric matrix of `n` rows and `basis@dimension` columns, with
#'   column names [basis_colnames()]. All zero above the smoothness of the
#'   family.
#'
#' @seealso [basis_eval()], which is `order = 0`; [basis_int()] for the
#'   opposite direction; [basis_is_numerical()] to learn whether a family
#'   answers this from a formula or a stencil.
#'
#' @examples
#' b <- bspline_basis(dimension = 6, degree = 3)
#' round(basis_deriv(b, c(0.25, 0.5, 0.75), order = 1), 3)
#'
#' # A cubic spline has three derivatives and then exact zeros.
#' all(basis_deriv(b, 0.5, order = 4) == 0)
#'
#' # Order 0 is the evaluation itself.
#' identical(basis_deriv(b, 0.5, order = 0), basis_eval(b, 0.5))
#'
#' # The derivative of a fit is the derivative of the basis times the same
#' # coefficients, so no refitting is involved.
#' beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
#' h <- 1e-5
#' fp <- (drop(basis_eval(b, 0.4 + h) %*% beta) -
#'        drop(basis_eval(b, 0.4 - h) %*% beta)) / (2 * h)
#' c(exact = drop(basis_deriv(b, 0.4, order = 1) %*% beta), difference = fp)
#'
#' # On a product basis the order is a multi-index, one entry per variable.
#' tb <- tensor_basis(bspline_basis(dimension = 4), poly_basis(dimension = 3))
#' round(basis_deriv(tb, cbind(0.5, 0.5), order = c(1, 1)), 4)
#' try(basis_deriv(tb, cbind(0.5, 0.5), order = 1))
#'
#' @export
basis_deriv <- S7::new_generic(
  "basis_deriv", "basis",
  function(basis, x, order = 1L, ...) {
    x <- check_eval_points(basis, x)
    order <- check_order(order, basis_nvar(basis))
    if (all(order == 0L)) {
      return(basis_eval(basis, x, ...))
    }
    S7::S7_dispatch()
  }
)


#' Integrate a Basis
#'
#' @description
#' Returns the definite integral of every basis function from the lower
#' endpoint of the basis interval up to each evaluation point, as a matrix of
#' the same shape [basis_eval()] returns. Times a coefficient vector it gives
#' the antiderivative of the fitted function that vanishes at the lower
#' endpoint.
#'
#' @details
#' # The anchored integral
#'
#' Column \eqn{j} of the result is
#'
#' \deqn{\int_{a}^{x} \varphi_j(t)\,\mathrm{d}t,}
#'
#' with \eqn{a} the lower endpoint of the basis interval, so the integral of
#' an expansion is the expansion against the same coefficients:
#' \eqn{\int_a^x \sum_j \beta_j \varphi_j = \sum_j \beta_j I_j(x)}.
#' The integral over the whole interval is the row at `basis@upper`, and the
#' integral over \eqn{[u, v]} is the difference of two rows.
#'
#' # Why the anchor is fixed
#'
#' The value at `basis@lower` is exactly zero, for every basis and every
#' column, and that is part of the contract every implementation owes. Any
#' antiderivative satisfies the differentiation check, so with the constant of
#' integration left free two bases could disagree while both being right, and
#' a sum of them would be wrong with nothing to report it.
#'
#' # A basis of several variables
#'
#' The integral is taken over the box from the lower corner to the point, one
#' iterated integral per variable, so on two variables column \eqn{j} is
#' \eqn{\int_{a_1}^{x_1}\int_{a_2}^{x_2} \varphi_j}.
#'
#' # The numerical fallback
#'
#' A subclass registering no method gets the one on the abstract [basis]
#' class: composite Gauss-Legendre from the lower endpoint to each point,
#' `nodes = 12` per panel by default. Exact for a polynomial integrand of
#' degree up to `2 * nodes - 1`, and accurate to the panel width elsewhere.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points inside the basis interval, the upper limits of
#'   the integrals: a numeric vector for a basis of one variable, or a matrix
#'   of [basis_nvar()] columns for a basis of several.
#' @param ... Passed to methods. The numerical fallback reads `nodes` from it,
#'   the number of Gauss-Legendre nodes per panel, default `12`.
#'
#' @return A numeric matrix of `n` rows and `basis@dimension` columns, with
#'   column names [basis_colnames()]. The row at `basis@lower` is exactly
#'   zero.
#'
#' @seealso [basis_eval()] and [basis_deriv()] for the other direction, and
#'   [basis_gram()] for integrals of products of basis functions.
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#'
#' # Zero at the lower endpoint, by the convention above.
#' basis_int(b, b@lower)
#' round(basis_int(b, c(0.25, 0.5, 1)), 4)
#'
#' # Differentiating the integral returns the basis.
#' bs <- bspline_basis(dimension = 6)
#' x <- c(0.2, 0.55, 0.9)
#' h <- 1e-6
#' fd <- (basis_int(bs, x + h) - basis_int(bs, x - h)) / (2 * h)
#' max(abs(fd - basis_eval(bs, x)))
#'
#' # The integral of a fitted curve, and the integral over a subinterval.
#' beta <- c(0.2, 1.1, -0.4, 0.8, 0.1, -0.6)
#' drop(basis_int(bs, 1) %*% beta)
#' drop((basis_int(bs, 0.75) - basis_int(bs, 0.25)) %*% beta)
#'
#' @export
basis_int <- S7::new_generic("basis_int", "basis", function(basis, x, ...) {
  x <- check_eval_points(basis, x)
  S7::S7_dispatch()
})


#' Gram Matrix of a Basis
#'
#' @description
#' Returns the matrix of inner products of the `order`-th derivatives of the
#' basis functions,
#' \deqn{G_{ab} = \int_{a}^{b} \varphi_a^{(d)}(t)\,\varphi_b^{(d)}(t)\,
#'   \mathrm{d}t,}
#' which is the matrix of a roughness penalty: \eqn{\beta^\top G_2 \beta} is
#' exactly \eqn{\int (f'')^2} for the function \eqn{f} the coefficients
#' describe. The three shipped families compute it in closed form, so no
#' quadrature error enters a penalized fit.
#'
#' @details
#' # What it is and what it is not
#'
#' The Gram matrix is an inner product of basis functions, so it depends on
#' the basis and the interval and on nothing else. It says which combinations
#' of coefficients are wiggly; how hard to shrink them is a separate decision
#' belonging to whatever fits the model.
#'
#' It is symmetric and positive semidefinite by construction. At `order = 0`
#' it is positive definite for a basis of linearly independent functions. At
#' any `order >= 1` it is singular, the constant differentiating to zero, and
#' its null space has dimension `order` for a family that contains the
#' polynomials of that degree.
#'
#' # Three measures
#'
#' The default is Lebesgue measure on the basis interval, computed from a
#' closed form by every shipped family.
#'
#' `at` replaces it with the empirical measure of the points given,
#' \eqn{B^{(d)\top} B^{(d)} / n}. That is the matrix a design matrix
#' produces, and the one to diagonalize against when the construction should
#' depend on where the data lie; [dr_basis()] uses it.
#'
#' `weight` takes a weighted Lebesgue measure \eqn{\int B B^\top w}, by
#' composite Gauss-Legendre over 50 panels of 12 nodes. A weight is an
#' arbitrary function, so no family has a closed form and the quadrature is
#' always run. Measured against the closed forms with \eqn{w \equiv 1}: 5e-15
#' for Fourier and Legendre at order 0, and 2e-11 at order 2; for a cubic
#' B-spline 5e-12 at order 0 and 2e-6 relative at order 2, where the second
#' derivative has kinks at knots that the panel breaks do not line up with.
#' Give at most one of `at` and `weight`; both together throws.
#'
#' # Where the arguments are handled
#'
#' `at` and `weight` are dealt with in the body of the generic, before
#' dispatch, so a method never sees either and returns the plain Lebesgue
#' matrix alone. A method must still carry both names in its signature, S7
#' requiring a method's formals to contain the generic's.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param order What the inner products are taken of. A single non-negative
#'   whole number, default `0`, giving the inner products of the basis
#'   functions themselves; `2` is the usual roughness penalty; one entry per
#'   variable for a basis of several. A [LinearOperator] instead gives
#'   \eqn{\int (Lb)(Lb)^\top}, the roughness matrix of that operator, and
#'   routes to [basis_operator_gram()].
#' @param at An optional numeric vector of points, or a matrix of
#'   [basis_nvar()] columns. When given, the inner products are taken against
#'   the empirical measure of those points and divided by their number.
#'   Missing values are dropped before evaluation; no usable point left
#'   throws.
#' @param weight An optional function of one numeric vector returning one
#'   non-negative value per point, a density to weight the integral by.
#'   Refused for a basis of several variables. A function returning the wrong
#'   length, an `NA` or a negative value throws.
#' @param ... Passed to methods, and to [weighted_gram()] when `weight` is
#'   given, where `panels` and `nodes` control the quadrature.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins. Positive semidefinite, and
#'   singular for any `order >= 1`.
#'
#' @seealso [basis_deriv()], whose columns it takes the inner products of;
#'   [orthonorm_basis()], which makes the `order = 0` matrix the identity;
#'   [dr_basis()], which diagonalizes one against the other.
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#'
#' # Diagonal at order 0, the Fourier functions being orthogonal.
#' round(basis_gram(b), 6)
#' round(basis_gram(b, order = 1), 4)
#'
#' # The order-2 matrix is the penalty matrix: beta' G beta is the integrated
#' # squared second derivative, which a fine trapezoid rule confirms to 6e-11.
#' p <- poly_basis(dimension = 5)
#' beta <- c(0.2, 1.1, -0.4, 0.8, 0.1)
#' tt <- seq(0, 1, length.out = 200001)
#' fpp <- drop(basis_deriv(p, tt, order = 2) %*% beta)^2
#' c(quadratic_form = drop(t(beta) %*% basis_gram(p, order = 2) %*% beta),
#'   quadrature = sum(head(fpp, -1) + tail(fpp, -1)) / 2 * (tt[2] - tt[1]))
#'
#' # Singular from order 1 up: order d leaves exactly d zero eigenvalues,
#' # the polynomials of degree below d differentiating away.
#' bs <- bspline_basis(dimension = 8)
#' vapply(0:3, function(d) {
#'   ev <- eigen(basis_gram(bs, order = d), only.values = TRUE)$values
#'   sum(abs(ev) < 1e-8 * max(abs(ev)))
#' }, integer(1))
#'
#' # Against the empirical measure of a sample instead of the interval.
#' set.seed(1)
#' round(basis_gram(b, at = runif(2000)), 3)
#'
#' # And against a weighted Lebesgue measure.
#' round(basis_gram(b, weight = function(x) dbeta(x, 2, 5)), 4)
#'
#' @export
basis_gram <- S7::new_generic(
  "basis_gram", "basis",
  function(basis, order = 0L, at = NULL, weight = NULL, ...) {
    # AN OPERATOR IS A SPELLING OF 'order', not a second argument beside it:
    # the question both answer is what the penalty differentiates, and a
    # whole number is the shorthand for deriv_operator() of it.
    if (is_operator(order)) {
      if (!is.null(at) && !is.null(weight)) {
        stop("Give at most one of 'at' and 'weight'.", call. = FALSE)
      }
      return(basis_operator_gram(basis, order, at = at, weight = weight, ...))
    }
    order <- check_order(order, basis_nvar(basis))
    if (!is.null(at) && !is.null(weight)) {
      stop("Give at most one of 'at' and 'weight'.", call. = FALSE)
    }
    if (!is.null(at)) return(empirical_gram(basis, order, at))
    if (!is.null(weight)) return(weighted_gram(basis, order, weight, ...))
    S7::S7_dispatch()
  }
)


#' The Roughness Matrix of a Differential Operator
#'
#' @description
#' The Gram matrix of \eqn{Lb}, that is
#' \deqn{R = \int_a^b (Lb)(Lb)^\top \, \mathrm{d}\mu,}
#' the matrix for which \eqn{\lVert Lx \rVert^2 = c^\top R c} when
#' \eqn{x = b^\top c}. It is what [basis_gram()] returns when its `order` is
#' a [LinearOperator], and the generic exists so that a family with a closed
#' form can declare one.
#'
#' @details
#' The base method integrates numerically, evaluating \eqn{Lb} through
#' [operator_eval()] at Gauss-Legendre nodes, at the covariate values for the
#' empirical measure, or against a weight function. A [FourierBasis] over a
#' full period overrides it with an exact diagonal form; see
#' [basis_operator_gram.FourierBasis()].
#'
#' @param basis A [basis].
#' @param op A [LinearOperator], with its period resolved.
#' @param at The covariate values, for the empirical measure, or `NULL`.
#' @param weight A density to integrate against, or `NULL`.
#' @param ... Passed to methods, and on to the quadrature (`panels`,
#'   `nodes`).
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns.
#'
#' @seealso [basis_gram()], the generic that routes to it;
#'   [operator_eval()] for \eqn{Lb}.
#'
#' @examples
#' b <- fourier_basis(lower = 0, upper = 365, dimension = 7)
#' round(diag(basis_gram(b, order = harmonic_operator(365))), 6)
#' @export
basis_operator_gram <- S7::new_generic(
  "basis_operator_gram", "basis",
  function(basis, op, at = NULL, weight = NULL, ...) S7::S7_dispatch()
)


#' Gram Matrix Against the Empirical Measure
#'
#' @description
#' Computes \eqn{B^{(d)\top} B^{(d)} / n} at the given points: the inner
#' products a design matrix produces, in place of those of the functions on
#' their interval. Called from the body of [basis_gram()] when `at` is
#' supplied, so no method ever sees this case.
#'
#' @details
#' Missing points are dropped **before** the basis is evaluated. A basis is
#' entitled to refuse a vector that is entirely missing, and its refusal would
#' name the wrong thing here, so `at` with no usable point throws
#' `'at' has no usable points.` instead. For a basis of several variables `at`
#' is coerced to a matrix and rows with any missing entry are dropped whole.
#'
#' The result is symmetrized as `(G + t(G))/2` before it is returned, the two
#' triangles of a crossproduct differing in their last bits, and given
#' [basis_colnames()] on both margins.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param order The derivative order, already validated by [check_order()].
#' @param at A numeric vector of points, or a matrix of [basis_nvar()]
#'   columns. Not range-checked here; [basis_deriv()] does that and throws for
#'   a point outside the interval.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins.
#'
#' @seealso [basis_gram()], its only caller, and [weighted_gram()] for the
#'   other alternative measure.
#'
#' @keywords internal
empirical_gram <- function(basis, order, at) {
  if (!is.numeric(at)) stop("'at' must be numeric.", call. = FALSE)
  # Dropped before evaluating rather than after: a basis is entitled to refuse
  # a vector that is entirely missing, and the refusal would name the wrong
  # thing here.
  if (basis_nvar(basis) > 1L) {
    at <- as.matrix(at)
    at <- at[stats::complete.cases(at), , drop = FALSE]
    if (!nrow(at)) stop("'at' has no usable points.", call. = FALSE)
  } else {
    at <- at[!is.na(at)]
    if (!length(at)) stop("'at' has no usable points.", call. = FALSE)
  }
  b <- basis_deriv(basis, at, order = order)
  g <- crossprod(b) / nrow(b)
  g <- (g + t(g)) / 2
  nm <- basis_colnames(basis)
  dimnames(g) <- list(nm, nm)
  g
}


#' Gram Matrix Against a Weighted Lebesgue Measure
#'
#' @description
#' Computes \eqn{\int_a^b B^{(d)}(t)\, B^{(d)}(t)^\top w(t)\,\mathrm{d}t}
#' by composite Gauss-Legendre. A weight is an arbitrary function, so no
#' family has a closed form for it and the quadrature is always run. Called
#' from the body of [basis_gram()] when `weight` is supplied.
#'
#' @details
#' The interval is cut into `panels` equal pieces and an `nodes`-point
#' Gauss-Legendre rule is placed on each, which integrates a polynomial of
#' degree up to `2 * nodes - 1` exactly on every panel. The weight is folded
#' into the quadrature weights and the matrix formed as a crossproduct of
#' \eqn{\sqrt{w_i}\,B^{(d)}(t_i)}, which keeps the result positive
#' semidefinite whatever the weight does.
#'
#' Measured against the closed forms at \eqn{w \equiv 1} with the defaults,
#' worst absolute entry: 5e-15 at order 0 and 2e-11 at order 2 for Fourier and
#' Legendre; 5e-12 at order 0 and 1.4e-3 at order 2 for a cubic B-spline over
#' eight knots, 2e-6 of the matrix's own scale, the second derivative there
#' having kinks the panel breaks do not line up with. Raise `panels` when a
#' family's derivative is not smooth.
#'
#' A basis of several variables is refused: the rule above is
#' one-dimensional.
#'
#' @param basis A basis object of one variable, of any class inheriting from
#'   [basis]. More than one variable throws.
#' @param order The derivative order, already validated by [check_order()].
#' @param weight A function of one numeric vector returning one non-negative
#'   value per point. A wrong length, an `NA` or a negative value throws.
#' @param panels The number of equal subintervals, default `50`.
#' @param nodes The number of Gauss-Legendre nodes per subinterval, default
#'   `12`.
#' @param ... Unused.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins.
#'
#' @seealso [basis_gram()], its only caller, and [empirical_gram()] for the
#'   other alternative measure.
#'
#' @keywords internal
weighted_gram <- function(basis, order, weight, panels = 50L, nodes = 12L, ...) {
  if (!is.function(weight)) {
    stop("'weight' must be a function of one numeric vector.", call. = FALSE)
  }
  if (basis_nvar(basis) > 1L) {
    stop("'weight' is not supported for a basis of several variables.",
      call. = FALSE
    )
  }
  breaks <- seq(basis@lower, basis@upper, length.out = panels + 1L)
  r <- quad_rule(breaks, nodes)
  w <- weight(r$nodes)
  if (length(w) != length(r$nodes) || anyNA(w) || any(w < 0)) {
    stop(
      "'weight' must return one non-negative value per point.",
      call. = FALSE
    )
  }
  b <- basis_deriv(basis, r$nodes, order = order)
  g <- crossprod(sqrt(r$weights * w) * b)
  g <- (g + t(g)) / 2
  nm <- basis_colnames(basis)
  dimnames(g) <- list(nm, nm)
  g
}


#' Column Names of a Basis Matrix
#'
#' @description
#' Returns the names every matrix the basis produces carries, so that the
#' evaluation, the derivatives of every order, the anchored integral and the
#' Gram matrix of one basis agree on their columns. The default numbers the
#' functions after the family; a family whose functions have identities of
#' their own overrides it.
#'
#' @details
#' The default method takes the first two characters of `@basis_name` and
#' appends `1` to `@dimension`, giving `bs1 ... bs6` for a B-spline.
#' [fourier_basis()] and [poly_basis()] override it with names carrying
#' meaning: `const`, `sin1`, `cos1`, `sin2` for the first, `P0`, `P1`, `P2`
#' for the second.
#'
#' A wrapper numbers its own columns under a short prefix, so an
#' orthonormalized basis reads `on1 ... on5`. A [tensor_basis()] pastes its
#' margins' names, one term per pair, as `bs1.const`, `bs1.sin1`,
#' `bs2.const`, so a coefficient's name says which marginal function it
#' belongs to in each variable.
#'
#' Overriding it is how a subclass gives its columns meaning. The method must
#' return exactly `basis@dimension` strings; [name_columns()] sets them
#' without checking, so a shorter vector is recycled by R and silently
#' mislabels.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param ... Passed to methods. No shipped family reads anything from it.
#'
#' @return A character vector of length `basis@dimension`.
#'
#' @seealso [basis_eval()], whose columns these name, and [print.basis()] for
#'   the object's summary.
#'
#' @examples
#' # The default numbers the functions after the family name.
#' basis_colnames(bspline_basis(dimension = 4))
#'
#' # Fourier and Legendre name theirs instead.
#' basis_colnames(fourier_basis(dimension = 5))
#' basis_colnames(poly_basis(dimension = 4))
#'
#' # Every matrix the basis produces carries the same names.
#' b <- bspline_basis(dimension = 4)
#' identical(colnames(basis_eval(b, 0.5)), basis_colnames(b))
#' identical(colnames(basis_gram(b, order = 2)), basis_colnames(b))
#'
#' @export
basis_colnames <- S7::new_generic("basis_colnames", "basis", function(basis, ...) {
  S7::S7_dispatch()
})

#' @name basis_colnames.basis
#' @title Default Column Names
#'
#' @description
#' Numbers the basis functions after the family, taking the first two
#' characters of `@basis_name` and appending `1` to `@dimension`: `bs1 ... bs6`
#' for a B-spline, `on1 ... on5` for an orthonormalized basis. The method
#' every class inherits unless it registers one of its own, as
#' [fourier_basis()], [poly_basis()] and [tensor_basis()] do.
#'
#' @details
#' Two characters is enough to tell the families apart at a glance in a
#' coefficient table without making the names long. Nothing depends on the
#' names being distinct across bases, and a model combining two B-spline
#' blocks will see `bs1` twice unless whatever assembles the design
#' disambiguates them.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A character vector of length `basis@dimension`.
#'
#' @seealso [basis_colnames()] for the generic.
#'
#' @keywords internal
S7::method(basis_colnames, basis) <- function(basis, ...) {
  paste0(substr(basis@basis_name, 1L, 2L), seq_len(basis@dimension))
}


#' Validate a Derivative Order
#'
#' @description
#' Checks that `order` is a non-negative whole number, or a vector of them of
#' length `nvar`, and returns it as an integer vector of length `nvar`. Called
#' from the bodies of [basis_deriv()] and [basis_gram()], so both report the
#' same errors in the same words.
#'
#' @details
#' Anything failing the first test throws
#' `'order' must be a non-negative integer.`; that covers a negative value, a
#' fraction, an infinity, an `NA` and a non-numeric.
#'
#' A vector of length `nvar` is returned as it stands, and a single `0` is
#' repeated to that length. For a basis of several variables a single
#' **non-zero** order throws with a longer message, because a scalar has two
#' readings there: that order in every coordinate, or that total order.
#' Choosing one silently would fit a different model from the one the caller
#' wrote. Zero is exempt, meaning no derivative under either reading.
#'
#' @param order The value supplied by the caller: a single non-negative whole
#'   number, or a vector of `nvar` of them.
#' @param nvar The number of variables the basis takes, from [basis_nvar()].
#'   Default `1`.
#'
#' @return `order` as an integer vector of length `nvar`.
#'
#' @seealso [basis_deriv()] and [basis_gram()], its two callers.
#'
#' @keywords internal
check_order <- function(order, nvar = 1L) {
  ok <- is.numeric(order) && length(order) >= 1L && all(is.finite(order)) &&
    all(order >= 0) && all(order == round(order))
  if (!ok) {
    stop("'order' must be a non-negative integer.", call. = FALSE)
  }
  order <- as.integer(order)

  if (length(order) == nvar) return(order)
  if (length(order) == 1L && (nvar == 1L || order == 0L)) {
    return(rep(order, nvar))
  }
  stop(sprintf(
    paste0(
      "'order' must have one entry per variable (%d), or be 0. A single ",
      "non-zero order is ambiguous for a basis of several variables: it ",
      "could mean that order in each coordinate, or that total order."
    ),
    nvar
  ), call. = FALSE)
}
