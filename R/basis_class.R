#' Basis Expansion
#'
#' @description
#' The abstract parent of every basis in the package. A basis is a finite
#' collection of functions on an interval, and the object carries what is
#' needed to evaluate that collection, differentiate it, integrate it, and
#' compute its inner products.
#'
#' @details
#' A basis of dimension \eqn{d} on \eqn{[a, b]} is the collection
#' \eqn{\varphi_1, \dots, \varphi_d}, and the object exists so that a
#' function may be written as a linear combination of them,
#'
#' \deqn{f(x) = \sum_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\beta,
#'   \qquad B(x)_{ij} = \varphi_j(x_i),}
#'
#' with \eqn{B(x)} the \eqn{n \times d} design matrix
#' \code{\link{basis_eval}} returns. Fitting \eqn{f} is then a linear
#' problem in \eqn{\beta} whatever the family, which is what makes
#' derivatives, the anchored integral and the Gram matrix properties of the
#' basis rather than of the fit.
#'
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
#' A basis lives on an interval, or, when it is a product of several, on a box:
#' \code{lower} and \code{upper} then have one entry per variable and
#' \code{\link{basis_nvar}} reports how many. Everything else is unchanged, and
#' a univariate basis is the case of one variable rather than a separate kind
#' of object.
#'
#' @param basis_name A short name for the family, used when printing.
#' @param dimension The number of functions in the basis, a positive integer.
#' @param lower,upper The endpoints of the interval the basis lives on, or one
#'   endpoint per variable for a basis of several.
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
    if (length(self@lower) < 1L || length(self@lower) != length(self@upper) ||
      !all(is.finite(self@lower)) || !all(is.finite(self@upper))) {
      return("@lower and @upper must be finite and of the same length")
    }
    if (any(self@lower >= self@upper)) {
      return("every @lower must be strictly less than its @upper")
    }
    NULL
  }
)


#' How Many Variables a Basis Takes
#'
#' @description
#' One for an ordinary basis, and the number of marginals for a product of
#' several.
#'
#' @details
#' The number is read off the endpoints, which carry one entry per variable, so
#' a basis declares its input dimension by construction rather than by saying
#' so separately and possibly disagreeing.
#'
#' @param basis An object inheriting from class \code{basis}.
#'
#' @return A positive integer.
#'
#' @examples
#' basis_nvar(bspline_basis(dimension = 5))
#' basis_nvar(tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3)))
#'
#' @seealso \code{\link{basis_colnames}}, \code{\link{basis_is_numerical}}
#' @export
basis_nvar <- function(basis) {
  length(basis@lower)
}


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
#' is rejected rather than extrapolated. Which extrapolation rule a model wants,
#' if any, is a decision for the layer that owns the meaning of the covariate; a
#' silent answer here would take that decision away from it.
#'
#' The comparison uses a tolerance relative to the width of the interval, so
#' that a point which is an endpoint up to rounding is accepted and then
#' clamped onto the endpoint exactly.
#'
#' For a basis of several variables the points are a matrix with one column per
#' variable, and each column is checked against its own endpoints. A basis of
#' one variable keeps taking, and returning, a plain vector.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points, or a matrix with one column
#'   per variable.
#'
#' @return \code{x}, with near-endpoint values clamped onto the endpoints: a
#'   vector for a basis of one variable and a matrix otherwise.
#'
#' @keywords internal
check_eval_points <- function(basis, x) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  d <- basis_nvar(basis)

  if (d == 1L) {
    if (is.matrix(x)) {
      if (ncol(x) != 1L) {
        stop("'x' must be a vector for a basis of one variable.", call. = FALSE)
      }
      x <- as.numeric(x)
    }
    return(clamp_to_range(x, basis@lower, basis@upper, "the basis interval"))
  }

  if (!is.matrix(x)) x <- matrix(x, ncol = d, byrow = TRUE)
  if (ncol(x) != d) {
    stop(sprintf(
      "'x' must have %d columns, one per variable of the basis.", d
    ), call. = FALSE)
  }
  for (j in seq_len(d)) {
    x[, j] <- clamp_to_range(
      x[, j], basis@lower[j], basis@upper[j],
      sprintf("the range of variable %d", j)
    )
  }
  x
}


#' Reject Points Outside a Range, and Clamp Those On Its Edge
#'
#' @param z A numeric vector.
#' @param lo,hi The endpoints.
#' @param what A phrase naming the range, used in the error message.
#'
#' @return \code{z}, clamped.
#'
#' @keywords internal
clamp_to_range <- function(z, lo, hi, what) {
  tol <- 1e-8 * (hi - lo)
  bad <- !is.na(z) & (z < lo - tol | z > hi + tol)
  if (any(bad)) {
    stop(sprintf(
      "%d of %d evaluation points fall outside %s [%s, %s].",
      sum(bad), length(z), what, format(lo), format(hi)
    ), call. = FALSE)
  }
  z[!is.na(z) & z < lo] <- lo
  z[!is.na(z) & z > hi] <- hi
  z
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
