#' @include numerical_fallbacks.R
NULL


#' B-Spline Basis
#'
#' @description
#' The S7 class of B-spline bases, the objects [bspline_basis()] returns. It
#' adds no property to [basis] and exists as the class the spline methods
#' dispatch on. A B-spline basis is piecewise polynomial, each function
#' supported on a few knot intervals, and its functions sum to one.
#'
#' @details
#' # The Cox-de Boor recurrence
#'
#' On a knot sequence \eqn{t_1 \le \cdots \le t_{d+m+1}} the functions are
#' defined from the indicators upward:
#'
#' \deqn{B_{j,0}(x) = \mathbf{1}\{t_j \le x < t_{j+1}\},}
#'
#' \deqn{B_{j,q}(x) = \frac{x - t_j}{t_{j+q} - t_j} B_{j,q-1}(x)
#'   + \frac{t_{j+q+1} - x}{t_{j+q+1} - t_{j+1}} B_{j+1,q-1}(x),
#'   \qquad q = 1, \dots, m,}
#'
#' a term with a zero denominator being taken as zero.
#'
#' # The two properties that follow
#'
#' \eqn{B_{j,m}} vanishes outside \eqn{[t_j, t_{j+m+1}]}, so at most
#' \eqn{m + 1} columns are non-zero in any row and the design matrix is
#' banded: at `degree = 3` exactly four of them, whatever the dimension. And
#' \eqn{\sum_j B_{j,m}(x) = 1} on the interval, measured to 2.2e-16, so the
#' basis carries its own constant and is collinear with an intercept in the
#' same design.
#'
#' # Where the numbers come from
#'
#' Evaluation, derivatives and integrals come from
#' [splines2::bSpline()], through the single wrapper [bspline_design()], which
#' computes all three from the recurrence and takes no differences. The Gram
#' matrix is integrated here instead, knot interval by knot interval, and is
#' exact: see [basis_gram.BsplineBasis()].
#'
#' # Its `basis_params`
#'
#' Three entries: `degree`, the interior `knots` as a numeric vector, and
#' `boundary_knots`, which is `c(lower, upper)`.
#'
#' @inheritParams basis
#'
#' @return An object of class `BsplineBasis`, inheriting from [basis], with
#'   the same five properties and `basis_params` holding `degree`, `knots` and
#'   `boundary_knots`. Call [bspline_basis()] instead of this class directly;
#'   it places the knots and checks the dimension against the degree.
#'
#' @references
#' de Boor, C. (2001). *A Practical Guide to Splines*, revised edition.
#' Springer.
#'
#' @seealso [bspline_basis()], the constructor;
#'   [basis_eval.BsplineBasis()], [basis_deriv.BsplineBasis()],
#'   [basis_int.BsplineBasis()] and [basis_gram.BsplineBasis()] for the
#'   methods registered on it.
#'
#' @examples
#' b <- bspline_basis(dimension = 5)
#' S7::S7_inherits(b, BsplineBasis)
#' b@basis_params
#'
#' # Local support: at most degree + 1 columns are non-zero in a row.
#' E <- basis_eval(bspline_basis(dimension = 12), seq(0.05, 0.95, by = 0.1))
#' rowSums(E != 0)
#'
#' # And the rows sum to one.
#' max(abs(rowSums(E) - 1))
#'
#' @export
BsplineBasis <- S7::new_class("BsplineBasis", parent = basis)


#' Construct a B-Spline Basis
#'
#' @description
#' Returns a basis of `dimension` B-splines of the given degree on
#' \eqn{[\ell, u]}, with the interior knots placed at equal spacing. This is
#' the general-purpose choice: the functions have local support, so a
#' coefficient moves the fitted curve only near its own knots, and the
#' conditioning does not deteriorate as the dimension grows.
#'
#' @details
#' # Dimension, degree and knots
#'
#' A basis of \eqn{K} functions of degree \eqn{m} has \eqn{K - m - 1}
#' interior knots, so \eqn{K \ge m + 1}; a smaller `dimension` throws, naming
#' both numbers. At equality there is no interior knot and the basis is the
#' polynomials of degree \eqn{m} on the whole interval, which it spans
#' exactly. The knots are `seq(lower, upper, length.out = K - m + 1)` with the
#' endpoints dropped, so they are equally spaced; a quantile placement is not
#' offered, and a caller wanting one can build the class directly.
#'
#' `degree = 0` gives indicator functions of the knot intervals, a step basis,
#' and `degree = 1` the piecewise linear hat functions.
#'
#' # The basis is complete
#'
#' All `dimension` functions are kept, so the rows of [basis_eval()] sum to
#' one and the basis spans the constant. Beside an intercept the design is
#' therefore rank deficient by one. Dropping a function is a linear
#' transformation of the basis: use [constrain_basis()], which keeps the
#' object a basis, and leave the choice of constraint to whatever owns the
#' meaning of the term.
#'
#' # Where each quantity comes from
#'
#' [basis_eval()], [basis_deriv()] and [basis_int()] call
#' [splines2::bSpline()], which evaluates the recurrence and its exact
#' derivative and integral. The Gram matrix is integrated in this package,
#' knot interval by knot interval with a rule sized from the degree, and is
#' exact to rounding. Nothing about a B-spline basis is differenced, and
#' [basis_is_numerical()] reports all three `FALSE`.
#'
#' @param lower,upper The endpoints of the interval, each a single finite
#'   number with `lower < upper`. Default \eqn{[0, 1]}. They are the boundary
#'   knots, and evaluating outside them throws.
#' @param dimension The number of basis functions, a single whole number of at
#'   least `degree + 1`, default `5`. It is the number of columns, so it fixes
#'   the flexibility of the fit; the interior knot count follows as
#'   `dimension - degree - 1`.
#' @param degree The degree of the piecewise polynomials, a single
#'   non-negative whole number, default `3` for cubic splines. `0` gives
#'   indicators of the knot intervals and `1` piecewise linear functions. A
#'   spline of degree \eqn{m} has \eqn{m} non-trivial derivatives; above that
#'   [basis_deriv()] returns exact zeros.
#'
#' @return An object of class [BsplineBasis], with `basis_name` `"bspline"`,
#'   `basis_params` holding `degree`, `knots` and `boundary_knots`, and column
#'   names `bs1`, `bs2`, and so on.
#'
#' @references
#' de Boor, C. (2001). *A Practical Guide to Splines*, revised edition.
#' Springer.
#'
#' @seealso [fourier_basis()] for a periodic basis and [poly_basis()] for a
#'   global polynomial one; [constrain_basis()] to remove the constant;
#'   [dr_basis()] to rotate this basis into the Demmler-Reinsch form a
#'   penalized fit uses; [basis_gram()] for its roughness penalty.
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' b
#'
#' # Local support: each function is non-zero on a few knot intervals only.
#' round(basis_eval(b, c(0.1, 0.5, 0.9)), 3)
#'
#' # The basis carries its own constant, so its rows sum to one.
#' max(abs(rowSums(basis_eval(b, seq(0, 1, length.out = 25))) - 1))
#'
#' # The second-derivative Gram matrix, the matrix of a roughness penalty.
#' round(basis_gram(b, order = 2), 2)
#'
#' # At dimension = degree + 1 there is no interior knot, and the basis is
#' # the polynomials of that degree: a cubic is fitted exactly.
#' p <- bspline_basis(dimension = 4, degree = 3)
#' length(p@basis_params$knots)
#' x <- seq(0, 1, length.out = 40)
#' max(abs(lm.fit(basis_eval(p, x), x^3)$fitted.values - x^3))
#'
#' # Degree 0 gives indicators of the knot intervals.
#' basis_eval(bspline_basis(dimension = 4, degree = 0), c(0.1, 0.3, 0.6, 0.9))
#'
#' # Too few functions for the degree is refused, with both numbers named.
#' try(bspline_basis(dimension = 3, degree = 3))
#'
#' @export
bspline_basis <- function(lower = 0, upper = 1, dimension = 5, degree = 3) {
  dimension <- check_basis_args(lower, upper, dimension)

  if (!is.numeric(degree) || length(degree) != 1L || !is.finite(degree) ||
    degree < 0 || degree != round(degree)) {
    stop("'degree' must be a single non-negative integer.", call. = FALSE)
  }
  degree <- as.integer(degree)

  if (dimension < degree + 1L) {
    stop(sprintf(
      paste0(
        "'dimension' (%d) is too small for 'degree' (%d): a B-spline basis ",
        "of degree m needs at least m + 1 functions."
      ),
      dimension, degree
    ), call. = FALSE)
  }

  n_interior <- dimension - degree - 1L
  interior <- if (n_interior > 0L) {
    seq(lower, upper, length.out = n_interior + 2L)[-c(1L, n_interior + 2L)]
  } else {
    numeric(0)
  }

  BsplineBasis(
    basis_name = "bspline",
    dimension = dimension,
    lower = lower,
    upper = upper,
    basis_params = list(
      degree = degree,
      knots = interior,
      boundary_knots = c(lower, upper)
    )
  )
}


#' Evaluate a B-Spline Basis
#'
#' @name basis_eval.BsplineBasis
#'
#' @description
#' Returns the B-spline design matrix at the given points, evaluated by
#' [splines2::bSpline()] from the Cox-de Boor recurrence. At most
#' `degree + 1` entries of any row are non-zero, and the row sums are one to
#' 2.2e-16.
#'
#' @details
#' The call goes through [bspline_design()] with `intercept = TRUE`, so all
#' `dimension` functions are returned and none is dropped for
#' identifiability. The result is stripped of the ten attributes
#' \pkg{splines2} attaches and of its `BSpline` class, leaving a plain
#' matrix, so a consumer never has to know where the numbers came from.
#'
#' @param basis A [BsplineBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names `bs1`, `bs2`, and so on.
#'
#' @seealso [bspline_design()], the one call into \pkg{splines2};
#'   [basis_eval()] for the generic.
#'
#' @keywords internal
S7::method(basis_eval, BsplineBasis) <- function(basis, x, ...) {
  name_columns(bspline_design(basis, x), basis)
}


#' Derivatives of a B-Spline Basis
#'
#' @name basis_deriv.BsplineBasis
#'
#' @description
#' Returns the `order`-th derivative of every B-spline, exactly, from the
#' derivative form of the Cox-de Boor recurrence in [splines2::bSpline()]. An
#' order above the degree short-circuits to the zero matrix, which is the
#' value of that derivative; nothing is thrown.
#'
#' @details
#' A spline of degree \eqn{m} is piecewise polynomial of that degree, so it
#' has \eqn{m} non-trivial derivatives and the rest vanish: a cubic gives
#' three, and `order = 4` is exactly zero everywhere. The short-circuit
#' happens here because \pkg{splines2} rejects a `derivs` above the degree
#' instead of returning zeros.
#'
#' The derivative of order \eqn{m} is a step function, discontinuous at each
#' interior knot, and the value returned at a knot is the one the recurrence
#' gives there. That matters for a Gram matrix, which is why
#' [basis_gram.BsplineBasis()] integrates knot interval by knot interval and
#' never across one.
#'
#' @param basis A [BsplineBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param order The derivative order, a single non-negative whole number,
#'   default `1`. Above `basis@basis_params$degree` the result is exactly
#'   zero.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names `bs1`, `bs2`, and so on.
#'
#' @seealso [bspline_design()], the one call into \pkg{splines2};
#'   [basis_deriv()] for the generic.
#'
#' @keywords internal
S7::method(basis_deriv, BsplineBasis) <- function(basis, x, order = 1L, ...) {
  if (order > basis@basis_params$degree) {
    out <- matrix(0, length(x), basis@dimension)
    out[is.na(x), ] <- NA_real_
    return(name_columns(out, basis))
  }
  name_columns(bspline_design(basis, x, derivs = order), basis)
}


#' Integral of a B-Spline Basis
#'
#' @name basis_int.BsplineBasis
#'
#' @description
#' Returns \eqn{\int_{\ell}^{x} B_j(t)\,\mathrm{d}t} for every function,
#' exactly, from [splines2::bSpline()] with `integral = TRUE`. That function
#' anchors its integral at the lower boundary knot, which is the same
#' convention [basis_int()] states, so no correction is applied here.
#'
#' @details
#' The integral of a spline of degree \eqn{m} is a spline of degree
#' \eqn{m + 1} on the same knots, and \pkg{splines2} evaluates it from the
#' recurrence, with no quadrature anywhere.
#'
#' The row at `basis@upper` holds the area under each function. Because the
#' basis is a partition of unity, those areas sum to the width of the
#' interval, which is a cheap check that the two conventions agree.
#'
#' @param basis A [BsplineBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, exactly zero in the row at `basis@lower`.
#'
#' @seealso [basis_int()] for the generic and the anchoring convention;
#'   [bspline_design()], the one call into \pkg{splines2}.
#'
#' @keywords internal
S7::method(basis_int, BsplineBasis) <- function(basis, x, ...) {
  name_columns(bspline_design(basis, x, integral = TRUE), basis)
}


#' Gram Matrix of a B-Spline Basis
#'
#' @name basis_gram.BsplineBasis
#'
#' @description
#' Returns the inner products of the `order`-th derivatives exactly, by
#' integrating over one knot interval at a time with a Gauss-Legendre rule
#' sized so that it reproduces the integrand exactly. This is the matrix a
#' roughness penalty on a spline is built from, so its exactness is worth the
#' small amount of work.
#'
#' @details
#' # Why it is exact
#'
#' On one knot interval the order-\eqn{d} derivative of a spline of degree
#' \eqn{m} is a polynomial of degree \eqn{m - d}, so the integrand
#' \eqn{B_a^{(d)} B_b^{(d)}} has degree \eqn{2(m - d)}. A Gauss-Legendre
#' rule with \eqn{m - d + 1} nodes is exact to degree \eqn{2(m - d) + 1},
#' which is one higher, so the only error left is floating point. Against the
#' same knot-aligned construction run at 20 nodes instead, the worst entry
#' agrees to 3.1e-16, 2.1e-14, 8.0e-13 and 7.3e-12 at orders 0 to 3 on a
#' cubic basis of six functions.
#'
#' The breaks are the boundary knots and the interior knots, so no panel
#' straddles a knot, and the exactness claim rests on that: at `order = m` the
#' derivative is a step function, and a rule spanning a knot would integrate
#' the wrong thing. Measured on the same basis, the general
#' `weight` route of [basis_gram()], whose 50 equal panels do not line up with
#' the knots, is out by 1.4e-3 at order 2 and by 65 at order 3.
#'
#' # Above the degree
#'
#' An `order` above `degree` returns the zero matrix, every derivative having
#' vanished. Below it the matrix is singular with an `order`-dimensional null
#' space, the polynomials of lower degree differentiating away.
#'
#' @param basis A [BsplineBasis] object.
#' @param order The derivative order, a single non-negative whole number,
#'   default `0`. `2` is the usual roughness penalty for a cubic.
#' @param at,weight Handled in the body of [basis_gram()] before dispatch, so
#'   they never arrive here. Named only because S7 requires a method's formals
#'   to contain the generic's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with column names `bs1`, `bs2`, and so on. Banded, entry \eqn{(a, b)}
#'   being zero whenever the two supports do not overlap.
#'
#' @seealso [basis_gram()] for the generic and the alternative measures;
#'   [quad_rule()], which builds the composite rule.
#'
#' @keywords internal
S7::method(basis_gram, BsplineBasis) <- function(basis, order = 0L, at = NULL,
                                                 weight = NULL, ...) {
  nm <- basis_colnames(basis)
  degree <- basis@basis_params$degree

  if (order > degree) {
    return(matrix(0, basis@dimension, basis@dimension, dimnames = list(nm, nm)))
  }

  breaks <- unique(c(basis@lower, basis@basis_params$knots, basis@upper))
  r <- quad_rule(breaks, n = degree - order + 1L)
  b <- basis_deriv(basis, r$nodes, order = order)
  g <- crossprod(sqrt(r$weights) * b)
  g <- (g + t(g)) / 2
  dimnames(g) <- list(nm, nm)
  g
}


#' The Roughness Matrix of an Operator on a B-Spline Basis
#'
#' @name basis_operator_gram.BsplineBasis
#'
#' @description
#' Exact for the length measure, integrating knot interval by knot interval
#' as [basis_gram.BsplineBasis()] does, and with one guard in front of it: an
#' operator of order above the degree of the spline is rejected rather than
#' integrated.
#'
#' @details
#' # Why it is exact
#'
#' \eqn{Lb} is a linear combination of derivatives of a spline of degree
#' \eqn{p}, so it is piecewise polynomial of degree at most \eqn{p} with the
#' same breakpoints, and \eqn{(Lb)(Lb)^\top} is piecewise of degree at most
#' \eqn{2p}. Gauss-Legendre with \eqn{p + 1} nodes on each knot interval is
#' exact to degree \eqn{2p + 1}, so nothing is approximated. The panels of
#' the base method's rule do not line up with the knots, where a derivative
#' of the spline jumps, and it is that misalignment rather than the node
#' count that costs the accuracy: measured at `dimension = 10`, `degree = 3`,
#' the base rule at 50 panels differs from the derivative route by 5e-8
#' relative, and this one by 3e-15.
#'
#' # The guard
#'
#' The \eqn{m}-th derivative of a spline of degree \eqn{p} is identically
#' zero for \eqn{m > p}, so the leading term of \eqn{Lb} vanishes and what
#' would be integrated is the operator with its highest derivative deleted.
#' That is a different penalty with a different null space, and nothing about
#' the result would say so. [basis_gram()] returns an exact zero matrix in
#' the same situation for a plain derivative, where the answer is at least
#' unmistakable; here it is not, so the refusal names the degree to raise.
#'
#' A spline does not contain the null space of a periodic operator exactly,
#' which is admissible and is not this guard's business: see the section on
#' [bspline_smooth()]'s page.
#'
#' @param basis A [basis].
#' @param op A [LinearOperator], with its period resolved.
#' @param at The covariate values, for the empirical measure, or `NULL`.
#' @param weight A density to integrate against, or `NULL`.
#' @param ... Passed on to the quadrature (`panels`, `nodes`).
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns.
#'
#' @seealso [basis_operator_gram()] for the generic.
#'
#' @examples
#' b <- bspline_basis(lower = 0, upper = 1, dimension = 10, degree = 3)
#' dim(basis_gram(b, order = harmonic_operator(1)))
#'
#' # an operator of order 5 on a cubic spline is refused, not truncated
#' try(basis_gram(b, order = harmonic_operator(1, harmonics = 2)))
#' @keywords internal
S7::method(basis_operator_gram, BsplineBasis) <- function(basis, op, at = NULL,
                                                          weight = NULL, ...) {
  check_operator(op)
  degree <- basis@basis_params$degree
  m <- operator_order(op)
  if (m > degree) {
    stop(sprintf(paste0(
      "the operator has order %d and the spline has degree %d, so D^%d of",
      " every\n  basis function is zero and the penalty would be the",
      " operator with its leading\n  term deleted. Raise 'degree' to at",
      " least %d."
    ), m, degree, m, m), call. = FALSE)
  }
  if (!is.null(at) || !is.null(weight)) {
    return(numerical_operator_gram(basis, op, at = at, weight = weight, ...))
  }
  breaks <- unique(c(basis@lower, basis@basis_params$knots, basis@upper))
  r <- quad_rule(breaks, n = degree + 1L)
  lb <- operator_eval(basis, r$nodes, op)
  g <- crossprod(sqrt(r$weights) * lb)
  g <- (g + t(g)) / 2
  nm <- basis_colnames(basis)
  dimnames(g) <- list(nm, nm)
  g
}


#' Call splines2 for a B-Spline Design Matrix
#'
#' @description
#' The single point at which this package talks to \pkg{splines2}. It
#' assembles the knot arguments from `basis@basis_params`, calls
#' [splines2::bSpline()] once, and returns a plain matrix, so the dependency
#' stays behind the S7 interface and no caller has to know its argument names
#' or its return class.
#'
#' @details
#' `intercept = TRUE` is passed always, so all `dimension` functions come
#' back; \pkg{splines2} would otherwise drop the first.
#'
#' The returned object is of class `BSpline` and carries ten attributes,
#' among them `x`, `knots`, `degree` and `intercept`. Rebuilding it as
#' `matrix(as.numeric(out), ...)` strips every one, which matters because
#' those attributes would survive arithmetic and reappear on a matrix that no
#' longer describes them. The dimensions are taken from `length(x)` and
#' `basis@dimension`, so a mismatch with what \pkg{splines2} returned surfaces
#' here.
#'
#' `derivs` and `integral` are mutually exclusive in practice, each caller
#' setting at most one.
#'
#' @param basis A [BsplineBasis] object.
#' @param x A numeric vector of evaluation points, already validated by the
#'   generic.
#' @param derivs The derivative order, default `0`. Must not exceed
#'   `basis@basis_params$degree`; [basis_deriv.BsplineBasis()] short-circuits
#'   above that and never calls here.
#' @param integral `TRUE` to return the integral anchored at the lower
#'   boundary knot instead of the functions. Default `FALSE`.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, of class `matrix` alone, with no dimnames and none of the
#'   attributes \pkg{splines2} attaches. Callers add the column names through
#'   [name_columns()].
#'
#' @seealso [splines2::bSpline()], the function called;
#'   [basis_eval.BsplineBasis()], [basis_deriv.BsplineBasis()] and
#'   [basis_int.BsplineBasis()], its three callers.
#'
#' @keywords internal
bspline_design <- function(basis, x, derivs = 0L, integral = FALSE) {
  p <- basis@basis_params
  out <- splines2::bSpline(
    x,
    knots = p$knots,
    Boundary.knots = p$boundary_knots,
    degree = p$degree,
    intercept = TRUE,
    derivs = derivs,
    integral = integral
  )
  matrix(as.numeric(out), nrow = length(x), ncol = basis@dimension)
}
