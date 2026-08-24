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
#' The identities verified, for a basis \eqn{\varphi_1, \dots, \varphi_d} on
#' \eqn{[a, b]}, are
#'
#' \deqn{\varphi_j^{(k)}(x) = \frac{\mathrm{d}}{\mathrm{d}x}\,
#'     \varphi_j^{(k-1)}(x), \qquad
#'   \frac{\mathrm{d}}{\mathrm{d}x} \int_a^x \varphi_j(t)\,\mathrm{d}t
#'     = \varphi_j(x), \qquad
#'   \int_a^a \varphi_j = 0,}
#'
#' together with \eqn{\sum_j \varphi_j(x) = 1} where the family has that
#' property and
#' \eqn{G_{jl} = \int \varphi_j \varphi_l \,\mathrm{d}\mu = G_{lj}} with
#' \eqn{G \succeq 0} for the Gram matrix.
#'
#' The checks are:
#'
#' 1. **shape**: every generic returns a matrix of the declared size,
#'    with the declared column names, for a vector and for a single point;
#' 2. **derivatives**: each analytic order agrees with one numerical
#'    differentiation of the order below it;
#' 3. **integral**: it differentiates back to the basis, and is
#'    exactly zero at the lower endpoint;
#' 4. **partition of unity**: the rows sum to one, for the families
#'    that have that property;
#' 5. **Gram**: symmetric, positive semidefinite, and equal to an
#'    independent quadrature;
#' 6. **missing values**: a missing evaluation point gives a missing
#'    row and nothing else.
#'
#' # What is skipped, and what is weakened
#'
#' Where a quantity comes from the numerical fallback, checking it against a
#' numerical reference would be the same arithmetic twice, agreeing however
#' wrong the basis is. The `deriv` and `integral` checks are therefore not run
#' at all in that case: they report `NA` and print `[numerical]`.
#'
#' The `gram` check is **weakened, not skipped**. Symmetry and positive
#' semidefiniteness are properties of the matrix whatever produced it, so they
#' are still tested and the row still prints `[PASSED]`; what is dropped is the
#' comparison against an independent quadrature. The `shape` and `missing`
#' checks read nothing about the fallback and always run.
#'
#' The independent quadrature is 401 panels of 7 nodes, a rule no basis in the
#' package uses for its own answer, so a family whose own Gram matrix is a
#' quadrature is still compared against different arithmetic.
#'
#' `partition` is `NA` and prints `[not claimed]` for a family that is not a
#' partition of unity. That is a different thing from `[numerical]`, and the
#' printout distinguishes them: a property the family never claimed was not
#' applicable, where a numerical value was simply not verified.
#'
#' # How the reference's own error is allowed for
#'
#' A central difference assumes derivatives the function may not have. At a
#' knot a spline's third derivative jumps, and a stencil straddling it returns
#' a number of the order of the jump, which compared against an exact value
#' reads as a failure of the basis. Each reference is therefore computed
#' twice, at a step and at half of it, and the gap between them bounds its own
#' uncertainty; the comparison is allowed that much slack point by point, so
#' each point contributes exactly the accuracy its reference supports. See
#' [fd_reference()].
#'
#' A deliberate error of five per cent is still caught by four orders of
#' magnitude, which is the check that the allowance has not blunted anything.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param n The number of points to test at, default `41`. They are equally
#'   spaced across the interval with five per cent trimmed off each end,
#'   because a one-sided stencil at an endpoint carries a larger error that
#'   would read as a failure of the basis. Passed through `as.integer()`.
#' @param orders The derivative orders to check, default `1:2`. Each is
#'   compared against one differentiation of the order below it, never against
#'   a chain of lower-order differences. Raising it tests more and costs one
#'   pass per order.
#' @param tol The relative tolerance for the derivative and integral checks,
#'   default `1e-6`. Relative to the values themselves, with the denominator
#'   taken from them; see [rel_close()]. The Gram comparison uses
#'   `1e-6` regardless.
#' @param verbose Whether to print the table, default `TRUE`. The result is
#'   returned invisibly either way.
#'
#' @return Invisibly, a named logical vector of length 6, in the order
#'   `shape`, `deriv`, `integral`, `partition`, `gram`, `missing`, with `NA`
#'   for a check that was not run. It carries the attribute `"numerical"`, the
#'   named logical vector [basis_is_numerical()] returns.
#'
#' @seealso [basis_is_numerical()] for the attribute and which route each
#'   quantity takes; `vignette("defining-a-basis")`, which uses this function
#'   to develop one.
#'
#' @examples
#' # Every shipped family passes all six checks.
#' invisible(check_basis(bspline_basis(dimension = 6)))
#' invisible(check_basis(fourier_basis(dimension = 5)))
#'
#' # The result is a logical vector, and the attribute says which quantities
#' # were computed from a formula.
#' r <- check_basis(poly_basis(dimension = 5), verbose = FALSE)
#' r
#' attr(r, "numerical")
#'
#' # A basis defined from its evaluation alone: the two checks that would
#' # compare a difference against a difference are not run, and the Gram check
#' # keeps its symmetry and definiteness half.
#' Bumps <- S7::new_class("Bumps", parent = basis)
#' S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
#'   out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
#'                           "-")^2 / 0.12^2)
#'   colnames(out) <- basis_colnames(basis)
#'   out
#' }
#' invisible(check_basis(Bumps(basis_name = "bumps", dimension = 4L,
#'                             lower = 0, upper = 1)))
#'
#' # A derivative five per cent wrong is caught.
#' Wrong <- S7::new_class("Wrong", parent = BsplineBasis)
#' S7::method(basis_deriv, Wrong) <- function(basis, x, order = 1L, ...) {
#'   1.05 * S7::method(basis_deriv, BsplineBasis)(basis, x, order = order)
#' }
#' b <- bspline_basis(dimension = 6)
#' invisible(check_basis(Wrong(basis_name = "wrong", dimension = b@dimension,
#'                             lower = 0, upper = 1,
#'                             basis_params = b@basis_params)))
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
  # through neighboring points, and at an endpoint it would have to switch to
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
#' Reports whether the family is a partition of unity, so that
#' [check_basis()] tests the row sums only where the property is claimed. A
#' family that is not one would fail a check it never promised to pass.
#'
#' @details
#' `TRUE` for a [BsplineBasis], which sums to one by construction, and for a
#' [TensorBasis] all of whose margins do: the row sums of a Kronecker product
#' are the products of the row sums.
#'
#' `FALSE` for everything else, including a [TransformedBasis] over a
#' B-spline. That is deliberate: an orthonormalization or a constraint takes
#' linear combinations of the columns, and the sum of the
#' new columns is generally not one. It is also conservative, so a
#' transformation that happens to preserve the property is untested rather
#' than wrongly failed.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#'
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [check_basis()], its only caller.
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
#' A matched pair for reading and writing one variable of a set of evaluation
#' points. `coord(x, j)` returns the `j`th variable; `replace_coord(x, j, z)`
#' returns the points with that variable replaced by `z` and the others left
#' alone. [check_basis()] uses them to sweep one coordinate at a time when the
#' basis takes several.
#'
#' @details
#' Both are the identity for a basis of one variable, where `x` is a plain
#' vector with no coordinate to pick: `coord()` returns `x` and
#' `replace_coord()` returns `z`, so the caller writes one loop over
#' `seq_len(basis_nvar(basis))` with no branch for the univariate case.
#'
#' Neither validates anything. `j` outside the columns of `x` gives R's own
#' subscript error, and a `z` of the wrong length is recycled by R's usual
#' rules.
#'
#' @param x A numeric vector of evaluation points, or a matrix of one column
#'   per variable.
#' @param j The coordinate to read or write, a column index. Ignored when `x`
#'   is a vector.
#' @param z The replacement values, for `replace_coord()` only: a numeric
#'   vector as long as `x` has rows.
#'
#' @return `coord()` returns a numeric vector: column `j` of `x`, or `x`
#'   itself when it is not a matrix. `replace_coord()` returns points of the
#'   same shape as `x`: the matrix with column `j` overwritten, or `z` itself
#'   when `x` is not a matrix.
#'
#' @seealso [check_basis()], their only caller, and [basis_nvar()] for the
#'   count they loop over.
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
#' Reports whether two matrices agree to a relative tolerance, with the
#' denominator taken from the values themselves and floored at a millionth of
#' the column's own scale. The one comparison [check_basis()] makes, so that
#' its derivative, integral and Gram checks all read agreement the same way.
#'
#' @details
#' # Why the denominator is not floored at one
#'
#' Flooring at one would flatten a disagreement between two small numbers into
#' apparent agreement, and small numbers are most of what a basis produces: a
#' B-spline is exactly zero on most of its interval. The denominator is
#' `pmax(abs(a), abs(b))`.
#'
#' # Why it is floored at all
#'
#' A basis function's derivative crosses zero, and at the crossing the
#' pointwise value vanishes while the numerical reference carries its usual
#' rounding error. Dividing that error by nothing reports a failure of the
#' reference as a failure of the basis. The floor is `1e-6` times the largest
#' absolute value in the same column, so a proportional error stays detectable
#' wherever the curve is large, which is where a wrong formula shows itself.
#'
#' # Columns with nothing to say
#'
#' A column whose whole scale is below `1e-8` of the largest column's is
#' skipped: neither side carries information there. If every column is
#' skipped the answer is `TRUE`.
#'
#' @param a,b Numeric matrices of the same shape. Neither is privileged; the
#'   comparison is symmetric.
#' @param tol The relative tolerance, a single positive number.
#' @param slack An optional numeric matrix of the same shape as `a`, an
#'   absolute allowance added entry by entry through `pmax(tol * den, slack)`.
#'   [check_basis()] passes [fd_reference()]'s `uncertainty` here, so a point
#'   whose finite-difference reference is unreliable is allowed the error that
#'   reference has. `NULL`, the default, allows none.
#'
#' @return A single `TRUE` or `FALSE`: `TRUE` when every informative entry
#'   agrees within its own allowance.
#'
#' @seealso [check_basis()], its caller, and [fd_reference()], which supplies
#'   `slack`.
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
#' Differentiates `f` once numerically and returns both the estimate and a
#' bound on its own error, entry by entry. [check_basis()] uses the second to
#' decide how much slack the comparison at each point deserves, so a point
#' where the reference is unreliable does not report the basis as wrong.
#'
#' @details
#' A central difference is only valid where the function has the derivatives
#' the stencil assumes. A spline does not: at a knot its third derivative
#' jumps, so a stencil that straddles the knot returns a number of the order of
#' the jump and not of the truncation error, and comparing an exact
#' analytical value against it reports a failure of the *reference*.
#'
#' Recomputing with the step halved says how much of the reference is error.
#' For a smooth point the two differ by about three quarters of the truncation,
#' so the gap between them bounds the reference's own uncertainty; at a knot it
#' is large. Nothing is discarded: the gap becomes the slack allowed to the
#' comparison, so each point contributes exactly the accuracy its reference
#' supports. This is the same device used elsewhere in the toolkit for a
#' parameter that is not differentiable, and it needs the same care: the two
#' estimates are compared with each other, not against a denominator floored at
#' one, since near a kink both are small and still differ by a factor.
#'
#' @param f A function of one numeric vector returning a numeric matrix with
#'   one row per point. [check_basis()] passes a closure over
#'   [basis_deriv()] or [basis_int()].
#' @param x A numeric vector of evaluation points.
#' @param lower,upper The endpoints of the interval, so that the stencil can
#'   be shifted to one side near an end instead of leaving the domain.
#'
#' @return A list of two matrices of the same shape as `f(x)`: `value`, the
#'   first derivative estimated at the full step, and `uncertainty`, four
#'   times the gap between that estimate and the one at half the step. An
#'   entry where either estimate is `NA` gets an infinite uncertainty, so the
#'   comparison there is unconstrained.
#'
#' @seealso [check_basis()], its only caller; [rel_close()], which consumes
#'   `uncertainty`; [numerical_deriv_matrix()], which computes each estimate.
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
#' @description
#' Prints the six-row table [check_basis()] shows when `verbose` is `TRUE`: a
#' header naming the family and its dimension, one line per check with a
#' description and a verdict, and a closing line listing the quantities that
#' came from the numerical fallback.
#'
#' @details
#' A verdict is `[PASSED]`, `[FAILED]`, or one of two things for an `NA`. A
#' check not run because the quantity is a finite difference prints
#' `[numerical]`; the partition-of-unity check on a family that is not one
#' prints `[not claimed]`. Saying which matters: a value that was not verified
#' is a gap, where a property never promised is not.
#'
#' Only `deriv`, `integral` and `partition` can be `NA` today. The other three
#' rows carry a `[numerical]` label in the table that nothing currently
#' reaches, `shape`, `gram` and `missing` all being assigned on every path.
#'
#' @param basis A basis object, of any class inheriting from [basis]. Read for
#'   its `@basis_name` and `@dimension` only.
#' @param res The named logical vector of results, in the order `shape`,
#'   `deriv`, `integral`, `partition`, `gram`, `missing`.
#' @param num The named logical vector [basis_is_numerical()] returns, printed
#'   as the closing line when any entry is `TRUE`.
#'
#' @return `NULL`, invisibly. Called for the printed output.
#'
#' @seealso [check_basis()], its only caller.
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
