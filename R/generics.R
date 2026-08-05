#' @include basis_class.R
NULL


#' Evaluate a Basis
#'
#' @description
#' Returns the design matrix of the basis at the given points: one row per
#' evaluation point, one column per basis function.
#'
#' @details
#' This is the only generic a basis must implement. Everything else in the
#' package has a numerical method registered on the \code{\link{basis}} class
#' and is therefore available from this one alone.
#'
#' The generic validates the evaluation points before dispatching, so every
#' method, including one written outside the package, refuses a point outside
#' the basis interval and receives points that are endpoints up to rounding
#' already clamped onto them.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns, with the column names the basis declares.
#'
#' @seealso \code{\link{basis_deriv}}, \code{\link{basis_int}},
#'   \code{\link{basis_gram}}
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' round(basis_eval(b, c(0, 0.5, 1)), 4)
#'
#' # a B-spline basis is a partition of unity
#' rowSums(basis_eval(b, seq(0, 1, length.out = 5)))
#'
#' @export
basis_eval <- S7::new_generic("basis_eval", "basis", function(basis, x, ...) {
  x <- check_eval_points(basis, x)
  S7::S7_dispatch()
})


#' Differentiate a Basis
#'
#' @description
#' Returns the \code{order}-th derivative of every basis function at the given
#' points, as a matrix of the same shape as \code{\link{basis_eval}}.
#'
#' @details
#' Derivative order is an argument rather than a family of generics because it
#' is unbounded: a Fourier basis is differentiable to any order, and a spline
#' of degree \eqn{k} has \eqn{k} non-trivial derivatives followed by zeros.
#' An order beyond what the family supports returns the zero matrix rather than
#' raising, which is the value of the derivative and not an omission.
#'
#' \code{order = 0} returns \code{\link{basis_eval}}, so that a loop over
#' orders needs no special case.
#'
#' A subclass that registers no method for this generic gets the numerical one
#' of the \code{\link{basis}} class, which applies a single central difference
#' stencil to \code{\link{basis_eval}}.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param order The derivative order, a non-negative integer.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#'
#' @seealso \code{\link{basis_eval}}, \code{\link{basis_is_numerical}}
#'
#' @examples
#' b <- bspline_basis(dimension = 6, degree = 3)
#' round(basis_deriv(b, c(0.25, 0.5, 0.75), order = 1), 3)
#'
#' # a cubic spline has no fourth derivative
#' all(basis_deriv(b, 0.5, order = 4) == 0)
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
#' endpoint of the basis interval up to each evaluation point.
#'
#' @details
#' The convention is fixed and is part of the contract: the value at
#' \code{basis@lower} is exactly zero, for every basis and every column. Any
#' antiderivative would satisfy the differentiation check, so without a fixed
#' constant of integration two bases could disagree while both being right,
#' and a sum of them would be wrong in a way nothing would report.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#'
#' @seealso \code{\link{basis_eval}}, \code{\link{basis_deriv}}
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#' basis_int(b, b@lower) # exactly zero, by the convention above
#' round(basis_int(b, c(0.25, 0.5, 1)), 4)
#'
#' @export
basis_int <- S7::new_generic("basis_int", "basis", function(basis, x, ...) {
  x <- check_eval_points(basis, x)
  S7::S7_dispatch()
})


#' Gram Matrix of a Basis
#'
#' @description
#' Returns the matrix of inner products of the \code{order}-th derivatives of
#' the basis functions,
#' \deqn{G_{ab} = \int_{\ell}^{u} B_a^{(d)}(t)\, B_b^{(d)}(t)\, \mathrm{d}t.}
#'
#' @details
#' The Gram matrix is what a roughness penalty integrates: the quadratic form
#' \eqn{\beta^\top G_2 \beta} is \eqn{\int (f'')^2}, for \eqn{f} the function
#' the coefficients describe. It is basis mathematics, an inner product, and
#' says nothing about which combination a model should be shrunk by.
#'
#' It is symmetric and positive semidefinite by construction, and singular
#' whenever the order-\eqn{d} derivatives are linearly dependent, which for
#' \eqn{d \ge 1} they always are: constants differentiate to zero.
#'
#' The inner product is taken against a measure, and which measure matters.
#' The default is Lebesgue on the basis interval, which is what a roughness
#' penalty integrates. Supplying \code{at} takes the empirical measure of those
#' points instead, \eqn{B^\top B / n}, which is the matrix a design matrix
#' actually produces and the one a basis is diagonalized against when the
#' construction is meant to depend on where the data lie. Supplying
#' \code{weight} takes a weighted Lebesgue measure.
#'
#' Both alternatives are handled in the body of the generic, before dispatch,
#' so a method never sees them and never has to implement them; it always
#' returns the plain Lebesgue matrix. A method must still name the arguments in
#' its signature, because S7 requires a method's formals to contain the
#' generic's.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order, a non-negative integer. Zero, the
#'   default, gives the inner products of the basis functions themselves.
#' @param at An optional numeric vector of points. When given, the inner
#'   products are taken against the empirical measure of those points rather
#'   than against Lebesgue measure.
#' @param weight An optional function of one numeric vector, a density to
#'   weight the integral by.
#' @param ... Passed to methods.
#'
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
#'
#' @seealso \code{\link{basis_deriv}}
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#' round(basis_gram(b), 6) # diagonal, by orthogonality
#' round(basis_gram(b, order = 1), 4)
#'
#' # against the empirical measure of a sample instead
#' set.seed(1)
#' round(basis_gram(b, at = runif(2000)), 3)
#'
#' @export
basis_gram <- S7::new_generic(
  "basis_gram", "basis",
  function(basis, order = 0L, at = NULL, weight = NULL, ...) {
    order <- check_order(order, basis_nvar(basis))
    if (!is.null(at) && !is.null(weight)) {
      stop("Give at most one of 'at' and 'weight'.", call. = FALSE)
    }
    if (!is.null(at)) return(empirical_gram(basis, order, at))
    if (!is.null(weight)) return(weighted_gram(basis, order, weight, ...))
    S7::S7_dispatch()
  }
)


#' Gram Matrix Against the Empirical Measure
#'
#' @description
#' \eqn{B^{(d)\top} B^{(d)} / n} at the given points: the inner products a
#' design matrix produces, rather than those of the functions on their
#' interval.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order.
#' @param at A numeric vector of points.
#'
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
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
#' \eqn{\int B^{(d)} B^{(d)\top} w(t)\, \mathrm{d}t}, by composite
#' Gauss-Legendre. A weight is an arbitrary function, so no family has a closed
#' form for it and the quadrature is always used.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order.
#' @param weight A function of one numeric vector.
#' @param panels The number of subintervals.
#' @param nodes The number of quadrature nodes per subinterval.
#' @param ... Unused.
#'
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
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
#' The names every matrix the basis produces carries. The default numbers the
#' functions after the family name; a subclass whose functions have their own
#' identities overrides it.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param ... Passed to methods.
#'
#' @return A character vector of length \code{basis@dimension}.
#'
#' @examples
#' basis_colnames(bspline_basis(dimension = 4))
#' basis_colnames(fourier_basis(dimension = 5))
#'
#' @export
basis_colnames <- S7::new_generic("basis_colnames", "basis", function(basis, ...) {
  S7::S7_dispatch()
})

#' @name basis_colnames.basis
#' @title Default Column Names
#' @description
#' Numbers the basis functions after the family name, as in \code{bs1},
#' \code{bs2}, and so on.
#' @param basis An object inheriting from class \code{basis}.
#' @param ... Unused.
#' @return A character vector of length \code{basis@dimension}.
#' @keywords internal
S7::method(basis_colnames, basis) <- function(basis, ...) {
  paste0(substr(basis@basis_name, 1L, 2L), seq_len(basis@dimension))
}


#' Validate a Derivative Order
#'
#' @description
#' A single non-negative integer for a basis of one variable, and one per
#' variable for a basis of several.
#'
#' @details
#' A multi-index is required rather than recycled from a scalar, because a
#' scalar has two plausible readings for a product basis -- that order in every
#' coordinate, or that total order -- and guessing between them would be a
#' silent choice. The single exception is zero, which means no derivative under
#' either reading.
#'
#' @param order The value supplied by the caller.
#' @param nvar The number of variables the basis takes.
#'
#' @return \code{order}, as an integer vector of length \code{nvar}.
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
