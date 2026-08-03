#' @include numerical_fallbacks.R
NULL


#' Validate a Basis
#'
#' @description
#' Runs the numerical checks that a basis must pass, and prints the outcome of
#' each. It is meant above all for a basis written outside the package, where a
#' hand-derived derivative or a misplaced constant of integration is the
#' likeliest mistake.
#'
#' @details
#' The checks are:
#' \enumerate{
#'   \item \strong{shape}: every generic returns a matrix of the declared size,
#'     with the declared column names, for a vector and for a single point;
#'   \item \strong{derivatives}: each analytic order agrees with one numerical
#'     differentiation of the order below it;
#'   \item \strong{integral}: it differentiates back to the basis, and is
#'     exactly zero at the lower endpoint;
#'   \item \strong{partition of unity}: the rows sum to one, for the families
#'     that have that property;
#'   \item \strong{Gram}: symmetric, positive semidefinite, and equal to an
#'     independent quadrature;
#'   \item \strong{missing values}: a missing evaluation point gives a missing
#'     row and nothing else.
#' }
#'
#' An order whose value comes from the numerical fallback is reported as
#' \code{[numerical]} rather than as passed. Checking such a value against a
#' numerical reference would be the same arithmetic twice, agreeing however
#' wrong the basis is, and a validator that reports agreement in that case is
#' worse than one that reports nothing: it says a thing was verified when it
#' was not.
#'
#' The Gram check is run against a rule the basis does not itself use, so that
#' a basis whose own method is quadrature is still compared with something
#' independent.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param n The number of points at which to test.
#' @param orders The derivative orders to check.
#' @param tol The relative tolerance for the derivative and integral checks.
#' @param verbose Whether to print the outcome.
#'
#' @return A named logical vector, invisibly, with \code{NA} for a check that
#'   was not run. The attribute \code{"numerical"} records which generics fell
#'   back.
#'
#' @seealso \code{\link{basis_is_numerical}}
#'
#' @examples
#' invisible(check_basis(bspline_basis(dimension = 6)))
#' invisible(check_basis(fourier_basis(dimension = 5)))
#'
#' @export
check_basis <- function(basis, n = 41L, orders = 1:2, tol = 1e-6,
                        verbose = TRUE) {
  num <- basis_is_numerical(basis)
  k <- basis@dimension
  nm <- basis_colnames(basis)

  # Points strictly inside the interval: the numerical reference differences
  # through neighbouring points, and at an endpoint it would have to switch to
  # a one-sided stencil, whose error is larger and would be read as a failure
  # of the basis rather than of the reference. The endpoints are covered by
  # the shape and integral checks, which need no reference.
  pad <- 0.05 * (basis@upper - basis@lower)
  x <- seq(basis@lower + pad, basis@upper - pad, length.out = n)

  res <- c(
    shape = NA, deriv = NA, integral = NA, partition = NA,
    gram = NA, missing = NA
  )

  ## 1. shapes and names
  b <- basis_eval(basis, x)
  res[["shape"]] <- is.matrix(b) && nrow(b) == n && ncol(b) == k &&
    identical(colnames(b), nm) &&
    identical(dim(basis_eval(basis, x[1])), c(1L, k)) &&
    identical(dim(basis_int(basis, x)), c(n, k)) &&
    identical(dim(basis_deriv(basis, x, order = 1L)), c(n, k))

  ## 2. derivatives: order k against one differentiation of order k - 1
  if (!num[["basis_deriv"]]) {
    ok <- TRUE
    for (d in orders) {
      ana <- basis_deriv(basis, x, order = d)
      ref <- numerical_deriv_matrix(
        function(z) basis_deriv(basis, z, order = d - 1L),
        x, 1L, basis@lower, basis@upper
      )
      ok <- ok && rel_close(ana, ref, tol)
    }
    res[["deriv"]] <- ok
  }

  ## 3. the integral: differentiates back, and is zero at the lower endpoint
  if (!num[["basis_int"]]) {
    back <- numerical_deriv_matrix(
      function(z) basis_int(basis, z), x, 1L, basis@lower, basis@upper
    )
    at_lower <- basis_int(basis, basis@lower)
    res[["integral"]] <- rel_close(back, basis_eval(basis, x), tol) &&
      all(at_lower == 0)
  }

  ## 4. partition of unity, where the family has it
  if (basis_partitions_unity(basis)) {
    res[["partition"]] <- max(abs(rowSums(b) - 1)) < 1e-10
  }

  ## 5. the Gram matrix
  g <- basis_gram(basis, order = 0L)
  sym <- max(abs(g - t(g))) < 1e-10
  psd <- min(eigen(g, symmetric = TRUE, only.values = TRUE)$values) > -1e-8
  if (num[["basis_gram"]]) {
    res[["gram"]] <- sym && psd
  } else {
    # A different rule from any the package uses for its own answer, so the
    # comparison is not the same arithmetic twice.
    ref <- numerical_gram(basis, 0L, panels = 401L, nodes = 7L)
    res[["gram"]] <- sym && psd && rel_close(g, ref, 1e-6)
  }

  ## 6. missing values travel through
  xm <- c(x[1:3], NA_real_)
  bm <- basis_eval(basis, xm)
  res[["missing"]] <- all(is.na(bm[4L, ])) && !anyNA(bm[1:3, ])

  if (verbose) print_basis_checks(basis, res, num)
  attr(res, "numerical") <- num
  invisible(res)
}


#' Does This Basis Sum to One?
#'
#' @description
#' Whether the family is a partition of unity, which
#' \code{\link{check_basis}} tests only for the families that claim it.
#'
#' @param basis An object inheriting from class \code{basis}.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
basis_partitions_unity <- function(basis) {
  S7::S7_inherits(basis, BsplineBasis)
}


#' Compare Two Matrices Relative to Their Own Magnitude
#'
#' @description
#' Whether two matrices agree to a relative tolerance, with the denominator
#' taken from the values themselves rather than floored at one.
#'
#' @details
#' Flooring the denominator at one would flatten a disagreement between two
#' small numbers into apparent agreement, which is exactly the region a basis
#' spends most of its time in: a B-spline is zero on most of its interval.
#' Values whose scale is at the level of rounding error are excluded instead of
#' being compared, since neither side carries information there.
#'
#' @param a,b Numeric matrices of the same shape.
#' @param tol The relative tolerance.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
rel_close <- function(a, b, tol) {
  scale <- pmax(abs(a), abs(b))
  big <- scale > 1e-8 * max(scale, na.rm = TRUE)
  if (!any(big)) return(TRUE)
  max(abs(a[big] - b[big]) / scale[big]) < tol
}


#' Print the Outcome of check_basis
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param res The logical vector of results.
#' @param num The logical vector of numerical fallbacks.
#'
#' @return \code{NULL}, invisibly.
#'
#' @keywords internal
print_basis_checks <- function(basis, res, num) {
  labels <- c(
    shape = "shapes and column names",
    deriv = "derivatives against finite differences",
    integral = "integral differentiates back, zero at lower",
    partition = "partition of unity",
    gram = "Gram symmetric, PSD, matches quadrature",
    missing = "missing values give missing rows"
  )
  # A check can be skipped for two different reasons, and saying which matters:
  # a value that came from the numerical fallback was not verified, whereas a
  # property the family never claimed was not applicable in the first place.
  skipped <- c(
    shape = "[numerical]", deriv = "[numerical]", integral = "[numerical]",
    partition = "[not claimed]", gram = "[numerical]", missing = "[numerical]"
  )
  cat("check_basis: ", basis@basis_name, " (", basis@dimension,
    " functions)\n",
    sep = ""
  )
  for (i in names(res)) {
    tag <- if (is.na(res[[i]])) {
      skipped[[i]]
    } else if (res[[i]]) {
      "[PASSED]"
    } else {
      "[FAILED]"
    }
    cat(sprintf("  %-11s %-46s %s\n", i, labels[[i]], tag))
  }
  if (any(num)) {
    cat("  computed numerically: ", paste(names(num)[num], collapse = ", "),
      "\n",
      sep = ""
    )
  }
  invisible(NULL)
}
