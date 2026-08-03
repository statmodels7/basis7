#' @include numerical_fallbacks.R
NULL


#' B-Spline Basis
#'
#' @description
#' The S7 class of B-spline bases. Constructed by
#' \code{\link{bspline_basis}}.
#'
#' @details
#' Evaluation, derivatives and integrals come from \pkg{splines2}, which
#' computes all three from the recurrence rather than by differencing. The
#' Gram matrix is integrated exactly, interval by interval: on each knot
#' interval the integrand is a polynomial of known degree, so a Gauss-Legendre
#' rule sized from that degree leaves no quadrature error.
#'
#' @inheritParams basis
#'
#' @return An object of class \code{BsplineBasis}. Use
#'   \code{\link{bspline_basis}} rather than calling the class directly, so
#'   that the knots are placed and the arguments checked.
#'
#' @seealso \code{\link{bspline_basis}}
#'
#' @examples
#' b <- bspline_basis(dimension = 5)
#' S7::S7_inherits(b, BsplineBasis)
#'
#' @export
BsplineBasis <- S7::new_class("BsplineBasis", parent = basis)


#' Construct a B-Spline Basis
#'
#' @description
#' A basis of B-splines of the given degree on \eqn{[\ell, u]}, with interior
#' knots placed at equal spacing.
#'
#' @details
#' The basis is complete: all \code{dimension} functions are kept, so the rows
#' of \code{\link{basis_eval}} sum to one. Dropping a function for
#' identifiability is a linear transformation of the basis and a decision for
#' the layer that knows what the term means.
#'
#' A basis of \eqn{K} functions with degree \eqn{m} has \eqn{K - m - 1}
#' interior knots, so \eqn{K \ge m + 1} is required, with equality giving the
#' polynomials of degree \eqn{m} on the whole interval and no interior knot at
#' all.
#'
#' @param lower,upper The endpoints of the interval.
#' @param dimension The number of basis functions, at least
#'   \code{degree + 1}.
#' @param degree The degree of the piecewise polynomials. Three, the default,
#'   gives cubic splines.
#'
#' @return An object of class \code{\link{BsplineBasis}}.
#'
#' @references
#' de Boor, C. (2001). \emph{A Practical Guide to Splines}. Springer.
#'
#' @seealso \code{\link{fourier_basis}}, \code{\link{check_basis}}
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' b
#'
#' # local support: each function is non-zero on a few knot intervals only
#' round(basis_eval(b, c(0.1, 0.5, 0.9)), 3)
#'
#' # the second-derivative Gram matrix, which a roughness penalty integrates
#' round(basis_gram(b, order = 2), 2)
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
#' @description
#' The design matrix of the B-spline recurrence, from
#' \code{\link[splines2]{bSpline}}.
#' @param basis A \code{\link{BsplineBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_eval, BsplineBasis) <- function(basis, x, ...) {
  name_columns(bspline_design(basis, x), basis)
}


#' Derivatives of a B-Spline Basis
#'
#' @name basis_deriv.BsplineBasis
#' @description
#' Exact derivatives from the B-spline recurrence. An order above the degree
#' gives the zero matrix, which is the value of that derivative.
#' @param basis A \code{\link{BsplineBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
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
#' @description
#' The exact integral from the lower boundary knot, from
#' \code{\link[splines2]{bSpline}}, which follows the same convention as this
#' package: the value at the lower endpoint is zero.
#' @param basis A \code{\link{BsplineBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_int, BsplineBasis) <- function(basis, x, ...) {
  name_columns(bspline_design(basis, x, integral = TRUE), basis)
}


#' Gram Matrix of a B-Spline Basis
#'
#' @name basis_gram.BsplineBasis
#' @description
#' Exact inner products, integrated knot interval by knot interval.
#' @details
#' On each knot interval the order-\eqn{d} derivative of a spline of degree
#' \eqn{m} is a polynomial of degree \eqn{m - d}, so the integrand of the Gram
#' matrix has degree \eqn{2(m - d)}. A Gauss-Legendre rule with \eqn{m - d + 1}
#' nodes integrates degree \eqn{2(m - d) + 1} exactly, so the result carries no
#' quadrature error, only floating-point error.
#' @param basis A \code{\link{BsplineBasis}} object.
#' @param order The derivative order.
#' @param ... Unused.
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
#' @keywords internal
S7::method(basis_gram, BsplineBasis) <- function(basis, order = 0L, ...) {
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


#' Call splines2 for a B-Spline Design Matrix
#'
#' @description
#' The single point at which this package talks to \pkg{splines2}, so that the
#' knot arguments are assembled once and the dependency stays behind the S7
#' interface.
#'
#' @param basis A \code{\link{BsplineBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param derivs The derivative order, or zero.
#' @param integral Whether to return the integral instead.
#'
#' @return A numeric matrix, stripped of the attributes \pkg{splines2} attaches.
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
