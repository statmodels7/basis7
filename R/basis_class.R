#' Basis Expansion
#'
#' @description
#' The abstract parent of every basis in the package. A basis is a finite
#' collection of functions on an interval, and the object carries what is
#' needed to evaluate that collection, differentiate it, integrate it, and
#' compute its inner products.
#'
#' @details
#' Concrete bases are subclasses. Each implements at least
#' \code{\link{basis_eval}}; \code{\link{basis_deriv}},
#' \code{\link{basis_int}} and \code{\link{basis_gram}} fall back to numerical
#' methods registered on this class, so a subclass that implements only its
#' evaluation is immediately complete. Every closed form registered later takes
#' over through dispatch, with no change to calling code.
#'
#' Bases are complete: a B-spline basis carries all its functions and its rows
#' sum to one. Restricting a basis, whether for identifiability or to separate
#' a linear from a nonlinear part, is a linear transformation of it and belongs
#' to the layer that owns that decision.
#'
#' @param basis_name A short name for the family, used when printing.
#' @param dimension The number of functions in the basis, a positive integer.
#' @param lower,upper The endpoints of the interval the basis lives on.
#' @param basis_params A named list of whatever else the subclass needs.
#'
#' @return An object inheriting from class \code{basis}.
#'
#' @seealso \code{\link{basis_eval}}, \code{\link{check_basis}},
#'   \code{\link{bspline_basis}}, \code{\link{fourier_basis}}
#'
#' @examples
#' # `basis` is abstract; construct a concrete subclass
#' b <- bspline_basis(dimension = 6)
#' b@dimension
#' c(b@lower, b@upper)
#'
#' @export
basis <- S7::new_class(
  "basis",
  abstract = TRUE,
  properties = list(
    basis_name = S7::class_character,
    dimension = S7::class_integer,
    lower = S7::class_numeric,
    upper = S7::class_numeric,
    basis_params = S7::class_list
  ),
  validator = function(self) {
    if (length(self@basis_name) != 1L) {
      return("@basis_name must be a single string")
    }
    if (length(self@dimension) != 1L || is.na(self@dimension) ||
      self@dimension < 1L) {
      return("@dimension must be a single positive integer")
    }
    if (length(self@lower) != 1L || length(self@upper) != 1L ||
      !is.finite(self@lower) || !is.finite(self@upper)) {
      return("@lower and @upper must be single finite numbers")
    }
    if (self@lower >= self@upper) {
      return("@lower must be strictly less than @upper")
    }
    NULL
  }
)


#' Validate the Arguments Every Basis Constructor Takes
#'
#' @description
#' Checks the interval and the number of functions, and returns the dimension
#' as an integer. Called by every constructor before anything family-specific.
#'
#' @param lower,upper The endpoints of the interval.
#' @param dimension The number of basis functions.
#'
#' @return \code{dimension}, as an integer.
#'
#' @keywords internal
check_basis_args <- function(lower, upper, dimension) {
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) ||
    !is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("'lower' and 'upper' must be single finite numbers.", call. = FALSE)
  }
  if (lower >= upper) {
    stop("'lower' must be strictly less than 'upper'.", call. = FALSE)
  }
  if (!is.numeric(dimension) || length(dimension) != 1L ||
    !is.finite(dimension) || dimension < 1 || dimension != round(dimension)) {
    stop("'dimension' must be a single positive integer.", call. = FALSE)
  }
  as.integer(dimension)
}


#' Validate Evaluation Points Against a Basis
#'
#' @description
#' Checks that \code{x} is numeric and lies inside the basis interval, and
#' returns it unchanged. Missing values are allowed and travel through to a
#' missing row.
#'
#' @details
#' A basis is defined on its interval and nowhere else, so a point outside it
#' is refused rather than extrapolated. Which extrapolation rule a model wants,
#' if any, is a decision for the layer that knows what the covariate means; a
#' silent answer here would take that decision away from it.
#'
#' The comparison uses a tolerance relative to the width of the interval, so
#' that a point which is an endpoint up to rounding is accepted and then
#' clamped onto the endpoint exactly.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points.
#'
#' @return \code{x}, with near-endpoint values clamped onto the endpoints.
#'
#' @keywords internal
check_eval_points <- function(basis, x) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  lo <- basis@lower
  hi <- basis@upper
  tol <- 1e-8 * (hi - lo)

  bad <- !is.na(x) & (x < lo - tol | x > hi + tol)
  if (any(bad)) {
    stop(sprintf(
      "%d of %d evaluation points fall outside the basis interval [%s, %s].",
      sum(bad), length(x), format(lo), format(hi)
    ), call. = FALSE)
  }
  x[!is.na(x) & x < lo] <- lo
  x[!is.na(x) & x > hi] <- hi
  x
}


#' Name the Columns of a Basis Matrix
#'
#' @description
#' Gives a basis matrix the column names the basis declares, so that
#' evaluation, derivatives and integrals of the same basis always agree on
#' them.
#'
#' @param m A numeric matrix with as many columns as the basis has functions.
#' @param basis An object inheriting from class \code{basis}.
#'
#' @return \code{m}, with column names set.
#'
#' @keywords internal
name_columns <- function(m, basis) {
  colnames(m) <- basis_colnames(basis)
  m
}
