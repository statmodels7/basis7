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
    order <- check_order(order)
    if (order == 0L) {
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
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order, a non-negative integer. Zero, the
#'   default, gives the inner products of the basis functions themselves.
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
#' @export
basis_gram <- S7::new_generic(
  "basis_gram", "basis",
  function(basis, order = 0L, ...) {
    order <- check_order(order)
    S7::S7_dispatch()
  }
)


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
#' @param order The value supplied by the caller.
#'
#' @return \code{order}, as an integer.
#'
#' @keywords internal
check_order <- function(order) {
  if (!is.numeric(order) || length(order) != 1L || !is.finite(order) ||
    order < 0 || order != round(order)) {
    stop("'order' must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(order)
}
