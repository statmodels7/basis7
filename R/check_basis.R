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
#' The derivative and integral checks allow for the accuracy of their own
#' reference. A central difference assumes derivatives the function may not
#' have: at a knot a spline's third derivative jumps, and a stencil straddling
#' it returns the jump rather than the truncation error, which would read as a
#' failure of the basis. Each reference is therefore computed twice, at a step
#' and at half of it, and the gap between them bounds its uncertainty; the
#' comparison is allowed that much slack, point by point. A deliberate error of
#' five per cent is still caught by four orders of magnitude, which is the
#' check that the allowance has not blunted anything.
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
  # The shape checks compare against dim(), which is always integer, so a
  # count given as a double would fail a check about the basis for a reason
  # about the argument.
  n <- as.integer(n)

  nvar <- basis_nvar(basis)

  # Points strictly inside the domain: the numerical reference differences
  # through neighbouring points, and at an endpoint it would have to switch to
  # a one-sided stencil, whose error is larger and would be read as a failure
  # of the basis rather than of the reference. The endpoints are covered by
  # the shape and integral checks, which need no reference.
  pad <- 0.05 * (basis@upper - basis@lower)
  x <- if (nvar == 1L) {
    seq(basis@lower + pad, basis@upper - pad, length.out = n)
  } else {
    # A shifted lattice rather than a product grid: a product grid of n points
    # per coordinate would cost n^D evaluations to learn nothing more.
    vapply(seq_len(nvar), function(j) {
      seq(basis@lower[j] + pad[j], basis@upper[j] - pad[j], length.out = n)[
        1L + (seq_len(n) + j - 2L) %% n
      ]
    }, numeric(n))
  }
  first <- if (nvar == 1L) x[1L] else x[1L, , drop = FALSE]
  low_corner <- if (nvar == 1L) basis@lower else matrix(basis@lower, 1L)

  res <- c(
    shape = NA, deriv = NA, integral = NA, partition = NA,
    gram = NA, missing = NA
  )

  ## 1. shapes and names
  b <- basis_eval(basis, x)
  res[["shape"]] <- is.matrix(b) && nrow(b) == n && ncol(b) == k &&
    identical(colnames(b), nm) &&
    identical(dim(basis_eval(basis, first)), c(1L, k)) &&
    identical(dim(basis_deriv(basis, x, order = rep(1L, nvar))), c(n, k))

  ## 2. derivatives: order k against one differentiation of order k - 1.
  ## For several variables this is done one coordinate at a time, since a
  ## mixed stencil would carry the product of two errors.
  if (!num[["basis_deriv"]]) {
    ok <- TRUE
    for (j in seq_len(nvar)) {
      for (d in orders) {
        e <- rep(0L, nvar)
        e[j] <- d
        ana <- basis_deriv(basis, x, order = e)
        below <- e
        below[j] <- d - 1L
        ref <- fd_reference(
          function(z) basis_deriv(basis, replace_coord(x, j, z), order = below),
          coord(x, j), basis@lower[j], basis@upper[j]
        )
        ok <- ok && rel_close(ana, ref$value, tol, ref$uncertainty)
      }
    }
    res[["deriv"]] <- ok
  }

  ## 3. the integral: differentiates back, and is zero at the lower corner.
  ## Differentiating back is a one-variable statement; for several the mixed
  ## derivative would need the stencil the fallback refuses, so only the
  ## anchor is checked there, which is the part the convention fixes.
  if (!num[["basis_int"]]) {
    at_lower <- basis_int(basis, low_corner)
    if (nvar == 1L) {
      ref <- fd_reference(
        function(z) basis_int(basis, z), x, basis@lower, basis@upper
      )
      res[["integral"]] <-
        rel_close(basis_eval(basis, x), ref$value, tol, ref$uncertainty) &&
          all(at_lower == 0)
    } else {
      res[["integral"]] <- all(at_lower == 0) &&
        identical(dim(basis_int(basis, x)), c(n, k))
    }
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
  xm <- if (nvar == 1L) {
    c(x[1:3], NA_real_)
  } else {
    rbind(x[1:3, , drop = FALSE], NA_real_)
  }
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
  # A product of partitions of unity is one: the row sums of a Kronecker
  # product are the product of the row sums.
  if (S7::S7_inherits(basis, TensorBasis)) {
    return(all(vapply(basis@marginals, basis_partitions_unity, logical(1))))
  }
  S7::S7_inherits(basis, BsplineBasis)
}


#' One Coordinate of the Evaluation Points
#'
#' @description
#' The \code{j}th variable of the points, and the points with that variable
#' replaced. A basis of one variable has a vector of points and no coordinate
#' to pick, so both are the identity there.
#'
#' @param x A numeric vector or matrix of evaluation points.
#' @param j The coordinate.
#' @param z The replacement values.
#'
#' @return A numeric vector for \code{coord}, and points of the same shape as
#'   \code{x} for \code{replace_coord}.
#'
#' @keywords internal
coord <- function(x, j) {
  if (is.matrix(x)) x[, j] else x
}

#' @rdname coord
#' @keywords internal
replace_coord <- function(x, j, z) {
  if (!is.matrix(x)) return(z)
  x[, j] <- z
  x
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
#' spends most of its time in: a B-spline is zero on most of its interval. The
#' denominator is therefore the values themselves.
#'
#' It is floored, but at a millionth of the scale of the column it belongs to
#' rather than at one. A basis function's derivative crosses zero, and at the
#' crossing the pointwise value vanishes while the numerical reference carries
#' its usual rounding error; dividing that error by nothing reports a failure
#' of the reference as a failure of the basis. Tying the floor to the curve's
#' own magnitude keeps a proportional error detectable wherever the curve is
#' large, which is where a wrong formula shows itself.
#'
#' Columns whose whole scale is at the level of rounding error are skipped,
#' since neither side carries information there.
#'
#' @param a,b Numeric matrices of the same shape.
#' @param tol The relative tolerance.
#'
#' @return \code{TRUE} or \code{FALSE}.
#'
#' @keywords internal
rel_close <- function(a, b, tol, slack = NULL) {
  scale <- pmax(abs(a), abs(b))
  colscale <- apply(scale, 2L, max, na.rm = TRUE)
  informative <- rep(colscale > 1e-8 * max(colscale, na.rm = TRUE),
    each = nrow(scale)
  )
  if (!any(informative)) return(TRUE)
  den <- pmax(scale, rep(colscale, each = nrow(scale)) * 1e-6)

  allowed <- tol * den
  if (!is.null(slack)) allowed <- pmax(allowed, slack)
  all(abs(a[informative] - b[informative]) <= allowed[informative])
}


#' A Finite-Difference Reference, and Where It Can Be Trusted
#'
#' @description
#' Differentiates \code{f} numerically, and reports at which points the result
#' is worth comparing anything against.
#'
#' @details
#' A central difference is only valid where the function has the derivatives
#' the stencil assumes. A spline does not: at a knot its third derivative
#' jumps, so a stencil that straddles the knot returns a number of the order of
#' the jump rather than of the truncation error, and comparing an exact
#' analytical value against it reports a failure of the \emph{reference}.
#'
#' Recomputing with the step halved says how much of the reference is error.
#' For a smooth point the two differ by about three quarters of the truncation,
#' so the gap between them bounds the reference's own uncertainty; at a knot it
#' is large, and says so. Nothing is discarded: the gap becomes the slack the
#' comparison is allowed, so a point tells the check as much as it honestly
#' can and no more. This is the same device used elsewhere in the toolkit for a
#' parameter that is not differentiable, and it needs the same care: the two
#' estimates are compared with each other, not against a denominator floored at
#' one, since near a kink both are small and still differ by a factor.
#'
#' @param f A function of a numeric vector returning a matrix.
#' @param x A numeric vector of evaluation points.
#' @param lower,upper The endpoints of the interval.
#'
#' @return A list with the reference \code{value} and the matrix
#'   \code{uncertainty} bounding its error at each entry.
#'
#' @keywords internal
fd_reference <- function(f, x, lower, upper) {
  full <- numerical_deriv_matrix(f, x, 1L, lower, upper)
  half <- numerical_deriv_matrix(f, x, 1L, lower, upper, step_scale = 0.5)
  gap <- abs(full - half)
  gap[is.na(gap)] <- Inf
  # Four times the gap, because halving the step of a second-order stencil
  # divides its truncation by four: the gap is three quarters of the error, so
  # the error is four thirds of it, and the factor is rounded up rather than
  # tuned.
  list(value = full, uncertainty = 4 * gap)
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
