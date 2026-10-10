#' @include generics.R
NULL




#' Numerically Differentiate a Matrix-Valued Function
#'
#' @description
#' Estimates the `order`-th derivative of `f` at each point of `x` by a single
#' finite-difference stencil, symmetric where the interval leaves room for it
#' and one-sided where it does not. One stencil of the order wanted, never a
#' composition of lower-order differences, so the error is the truncation of
#' that one formula. It computes the numerical derivatives of every basis; the
#' numerical integral and Gram matrix use Gauss-Legendre quadrature through
#' [quad_rule()] instead.
#'
#' @details
#' # Choice of the stencil
#'
#' A basis is evaluated at its endpoints as readily as anywhere else, and a
#' symmetric stencil centered on an endpoint would need points outside the
#' interval, where the basis signals an error. Each point therefore gets the
#' central stencil when both sides have room, and otherwise the one-sided
#' stencil that points inward. The one-sided stencils have `order + 2` nodes,
#' so that they keep the second-order accuracy of the central one, which at an
#' even order has one node fewer; the central stencil is padded with zero
#' weights to the same length. The weights come from
#' [numericals7::fd_weights()], and the central offsets from
#' [numericals7::fd_offsets()].
#'
#' The one-sided stencils carry a larger error constant. Measured against
#' exact Legendre derivatives of dimension 5 on \eqn{[0, 1]}, relative to the
#' scale of the result, the first derivative agrees to 5.1e-10 at the two
#' endpoints and to 3.9e-10 in the interior, the second to 1.3e-07 and
#' 1.5e-08.
#'
#' # The step
#'
#' The step is that of [numericals7::fd_step()],
#' \eqn{\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)}, which balances
#' truncation against rounding for order \eqn{d}. The rule is not repeated in
#' this package, so that it has a single definition in \pkg{numericals7}.
#'
#' It is then capped at `0.4 * (upper - lower) / (2 * reach)`, `reach` being
#' the half-width of the central stencil, so that every stencil fits inside
#' the interval: `0.2` of the width at orders 1 and 2, and `0.1` at orders 3
#' and 4.
#'
#' # What it costs in accuracy
#'
#' Measured against exact Legendre derivatives of dimension 5 at 21 interior
#' points of \eqn{[0, 1]}, relative to the scale of the result: 3.9e-10 at
#' order 1, 1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08 at order 4,
#' and about ten times more at the endpoints from order 2 on. A closed form,
#' where one exists, is exact.
#'
#' @param f A function of one numeric vector returning a numeric matrix with
#'   one row per element. Called `max(2 * reach + 1, order + 2)` times, from
#'   three times at order 1 to six at order 4.
#' @param x A numeric vector of evaluation points. `NA` entries are given the
#'   central stencil and propagate to an `NA` row.
#' @param order The derivative order, a single positive whole number.
#' @param lower,upper The endpoints of the interval on which `f` is defined,
#'   used to choose each point's stencil and to cap the step.
#' @param step_scale A multiplier on the step, default `1`. [fd_reference()]
#'   passes `0.5` to measure the uncertainty of the reference by the change in
#'   the result.
#'
#' @return A numeric matrix with `length(x)` rows and as many columns as `f`
#'   returns, with no dimnames; callers add them through [name_columns()].
#'
#' @seealso [numericals7::fd_weights()], [numericals7::fd_offsets()] and
#'   [numericals7::fd_step()], which supply the weights, the central offsets
#'   and the step;
#'   [basis_deriv.basis()], its main caller.
#'
#' @keywords internal
numerical_deriv_matrix <- function(f, x, order, lower, upper, step_scale = 1) {
  off <- numericals7::fd_offsets(order)
  reach <- off$reach

  # the step is numericals7's, like the nodes and the weights below: written
  # out here it would be a second copy of a rule that has one home
  h <- numericals7::fd_step(x, order, accuracy = 2L) * step_scale
  h <- pmin(h, 0.4 * (upper - lower) / (2 * reach))

  # Which stencil each point can afford: symmetric when both sides have room,
  # otherwise the one-sided stencil that points into the interval.
  room_left <- (x - lower) >= reach * h
  room_right <- (upper - x) >= reach * h
  kind <- ifelse(room_left & room_right, "c", ifelse(room_right, "f", "b"))
  kind[is.na(x)] <- "c"

  # The one-sided stencils have order + 2 nodes, so that they keep the
  # second-order accuracy of the central one: with 2 * reach + 1 nodes an even
  # order is only first-order accurate at the ends. The central stencil is
  # padded with zero weights to the same length.
  nnode <- max(2L * reach + 1L, order + 2L)
  pad <- nnode - (2L * reach + 1L)
  offs <- list(
    c = c(off$central, rep(0L, pad)),
    f = seq.int(0L, nnode - 1L),
    b = seq.int(-(nnode - 1L), 0L)
  )
  w <- list(
    c = c(numericals7::fd_weights(off$central, order), rep(0, pad)),
    f = numericals7::fd_weights(offs$f, order),
    b = numericals7::fd_weights(offs$b, order)
  )

  n <- length(x)
  out <- NULL
  for (k in seq_len(nnode)) {
    # USE.NAMES = FALSE, because `kind` is a character vector and vapply would
    # otherwise name the offsets after the stencil each point uses. Those names
    # travel into `x + s * h` and out again as the row names of whatever the
    # basis returns, so an evaluation method built on outer() would label its
    # rows "c" and "f".
    s <- vapply(kind, function(z) offs[[z]][k], numeric(1), USE.NAMES = FALSE)
    wk <- vapply(kind, function(z) w[[z]][k], numeric(1), USE.NAMES = FALSE)
    val <- f(x + s * h)
    if (is.null(out)) out <- matrix(0, n, ncol(val))
    out <- out + wk * val
  }
  out / h^order
}


#' Gauss-Legendre Nodes and Weights
#'
#' @description
#' Returns the `n`-point Gauss-Legendre rule on \eqn{[-1, 1]}: the nodes and
#' the weights that integrate every polynomial of degree up to \eqn{2n - 1}
#' exactly. At `n = 5` the rule reproduces \eqn{\int t^9} as `0` and first
#' departs at \eqn{t^{10}}, returning 0.17889 against 0.18182.
#'
#' @details
#' The construction is Golub-Welsch: the nodes are the eigenvalues of the
#' symmetric tridiagonal Jacobi matrix of the Legendre recurrence, with
#' off-diagonal \eqn{i/\sqrt{4i^2 - 1}}, and the weights are twice the square
#' of the first component of each eigenvector. The weights sum to 2, the
#' length of the interval.
#'
#' They are computed at call time, so every node count is available. The
#' exact spline rules need this, one rule per knot interval sized from the
#' degree and the derivative order, which a fixed table could not cover.
#'
#' @param n The number of nodes, a single positive whole number. `1` returns
#'   the midpoint rule directly. Any other value signals an error, which
#'   reaches the user through the `nodes` argument of the quadrature
#'   routines.
#'
#' @return A list of two numeric vectors of length `n`: `nodes`, in increasing
#'   order, and `weights`, positive and summing to 2.
#'
#' @references
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature
#' rules. *Mathematics of Computation* **23**, 221-230.
#'
#' @seealso [quad_rule()], which maps this rule onto a sequence of intervals.
#'
#' @keywords internal
gauss_legendre <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
    n != round(n)) {
    stop("the number of quadrature nodes must be a positive whole number.",
      call. = FALSE
    )
  }
  n <- as.integer(n)
  if (n == 1L) return(list(nodes = 0, weights = 2))
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- b
  jacobi[cbind(i + 1L, i)] <- b
  e <- eigen(jacobi, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = 2 * e$vectors[1L, ord]^2)
}


#' Map a Quadrature Rule onto Intervals
#'
#' @description
#' Places an `n`-point Gauss-Legendre rule on each interval between
#' consecutive breakpoints and returns the pooled nodes and weights, so that
#' `sum(weights * f(nodes))` is the integral over the whole span. The
#' composite rule every quadrature in the package is built from.
#'
#' @details
#' Each interval gets the same `n` nodes, affinely mapped from
#' \eqn{[-1, 1]}, and its weights scaled by half its width. The result is
#' exact for any function that is a polynomial of degree at most
#' \eqn{2n - 1} **on each interval separately**, which is why the callers
#' choose their breaks with care: [basis_gram.BsplineBasis()] uses the knots,
#' so no interval straddles the point where a spline's derivative jumps.
#'
#' The function does not validate its arguments. `breaks` must be
#' increasing and of length at least two. The nodes are ordered by their
#' position within the rule and then by interval: the first node of every
#' interval comes first, then the second node of every interval, and so on.
#'
#' @param breaks A numeric vector of at least two increasing breakpoints. The
#'   first and last are the ends of the span.
#' @param n The number of nodes per interval, a single positive whole number.
#'
#' @return A list of two numeric vectors of length
#'   `n * (length(breaks) - 1)`: `nodes` and `weights`, in the order described
#'   above.
#'
#' @seealso [gauss_legendre()], which supplies the rule on \eqn{[-1, 1]};
#'   [numerical_gram()] and [basis_int.basis()], which consume it.
#'
#' @keywords internal
quad_rule <- function(breaks, n) {
  gl <- gauss_legendre(n)
  lo <- breaks[-length(breaks)]
  hi <- breaks[-1L]
  half <- (hi - lo) / 2
  mid <- (hi + lo) / 2
  list(
    nodes = as.numeric(outer(half, gl$nodes) + mid),
    weights = as.numeric(outer(half, gl$weights))
  )
}


#' Numerical Derivatives of a Basis
#'
#' @name basis_deriv.basis
#' @title Numerical Derivatives of a Basis
#'
#' @description
#' The derivative method every basis inherits from the abstract [basis] class:
#' one finite-difference stencil of the requested order, applied to
#' [basis_eval()]. [basis_deriv()] is therefore available for a subclass that
#' supplies its evaluation alone, and a closed form registered later takes
#' over through dispatch with no change to calling code.
#'
#' @details
#' # Accuracy
#'
#' One stencil of the requested order is used, never a chain of first
#' differences. See [numerical_deriv_matrix()] for the stencil, the step, the
#' endpoint rule and the measured accuracy.
#'
#' # A basis of several variables
#'
#' The stencil differentiates along one coordinate at a time, replacing that
#' column of the points and holding the others. A mixed partial such as
#' `c(1, 1)` therefore signals an error: a stencil in the plane carries the
#' product of two errors, and [tensor_basis()], the family that needs mixed
#' partials, computes them exactly from its margins.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points inside the basis interval: a numeric vector for
#'   a basis of one variable, or a matrix of [basis_nvar()] columns for a
#'   basis of several.
#' @param order The derivative order, already checked by the generic: a
#'   single non-negative whole number, or one per variable. More than one
#'   non-zero entry signals an error.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names [basis_colnames()].
#'
#' @seealso [numerical_deriv_matrix()], which does the work;
#'   [basis_is_numerical()], which reports whether the derivatives of an
#'   object come from here; [basis_deriv()] for the generic.
#'
#' @keywords internal
S7::method(basis_deriv, basis) <- function(basis, x, order = 1L, ...) {
  d <- basis_nvar(basis)
  if (d == 1L) {
    out <- numerical_deriv_matrix(
      function(z) basis_eval(basis, z),
      x, order, basis@lower, basis@upper
    )
    return(name_columns(out, basis))
  }

  # For several variables the stencil differentiates along one coordinate at a
  # time. A mixed partial needs a stencil in the plane, whose error is the
  # product of two, and no basis in the package needs one: a tensor product
  # computes its own exactly. Refused rather than approximated badly.
  active <- which(order > 0L)
  if (length(active) > 1L) {
    stop(
      "The numerical fallback differentiates one variable at a time; a mixed ",
      "partial derivative has to be supplied by the basis.",
      call. = FALSE
    )
  }
  j <- active[1L]
  out <- numerical_deriv_matrix(
    function(z) {
      xj <- x
      xj[, j] <- z
      basis_eval(basis, xj)
    },
    x[, j], order[j], basis@lower[j], basis@upper[j]
  )
  name_columns(out, basis)
}


#' Numerical Integral of a Basis
#'
#' @name basis_int.basis
#' @title Numerical Integral of a Basis
#'
#' @description
#' The integration method every basis inherits from the abstract [basis]
#' class: composite Gauss-Legendre from the lower endpoint, accumulated over
#' the sorted evaluation points, so the whole set costs one pass instead of
#' one quadrature per point.
#'
#' @details
#' # Accumulation
#'
#' The points are sorted and made unique, the rule is placed on each segment
#' between consecutive ones, and the segment integrals are cumulated. A point
#' equal to `basis@lower` gives an empty first segment, so its row is exactly
#' zero and the anchoring convention of [basis_int()] holds by construction
#' with no cancellation behind it. Duplicated points are computed once and
#' matched back.
#'
#' # Accuracy
#'
#' The `nodes`-point rule integrates a polynomial of degree up to
#' `2 * nodes - 1` exactly on each segment, so on a polynomial family of
#' lower degree the result is exact to rounding. On a family that is not
#' polynomial the error is that of the rule on each segment, which shrinks
#' with the spacing of the evaluation points; a single distant point is
#' integrated by one rule over the whole span.
#'
#' # One variable only
#'
#' A basis of several variables signals an error. The integral there is over
#' a box, one iterated integral per variable, and a family that needs it
#' supplies its own, as [basis_int.TensorBasis()] does from its margins.
#'
#' @param basis A basis object of one variable, of any class inheriting from
#'   [basis]. More than one variable signals an error.
#' @param x A numeric vector of evaluation points inside the basis interval,
#'   the upper limits of the integrals. Need not be sorted or unique. `NA`
#'   gives an `NA` row.
#' @param nodes The number of Gauss-Legendre nodes per segment, default `12`,
#'   so a polynomial integrand of degree up to 23 is integrated exactly.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names [basis_colnames()], exactly zero in the row at
#'   `basis@lower`.
#'
#' @seealso [quad_rule()], which places the rule; [basis_int()] for the
#'   generic and the anchoring convention.
#'
#' @keywords internal
S7::method(basis_int, basis) <- function(basis, x, nodes = 12L, ...) {
  if (basis_nvar(basis) > 1L) {
    stop(
      "The numerical fallback integrates one variable; a multiple integral ",
      "over a box has to be supplied by the basis, as a tensor product does.",
      call. = FALSE
    )
  }
  k <- basis@dimension
  out <- matrix(NA_real_, length(x), k)
  ok <- !is.na(x)

  if (any(ok)) {
    u <- sort(unique(x[ok]))
    # Segment i runs from edges[i] to edges[i + 1]; the first starts at the
    # lower endpoint. A point equal to the lower endpoint gives an empty first
    # segment, so its integral is exactly zero rather than nearly so.
    edges <- c(basis@lower, u)
    seg <- matrix(0, length(u), k)
    for (i in seq_along(u)) {
      if (edges[i] < edges[i + 1L]) {
        r <- quad_rule(edges[c(i, i + 1L)], nodes)
        seg[i, ] <- colSums(r$weights * basis_eval(basis, r$nodes))
      }
    }
    acc <- matrix(apply(seg, 2L, cumsum), nrow = length(u))
    out[ok, ] <- acc[match(x[ok], u), , drop = FALSE]
  }
  name_columns(out, basis)
}


#' Numerical Gram Matrix of a Basis
#'
#' @name basis_gram.basis
#' @title Numerical Gram Matrix of a Basis
#'
#' @description
#' The inner-product method every basis inherits from the abstract [basis]
#' class: composite Gauss-Legendre over `panels` equal subintervals of the
#' interval, applied to the requested derivatives. A one-line wrapper over
#' [numerical_gram()], which does the work and which [basis_gram.FourierBasis()]
#' also calls when the period is not the width of the interval.
#'
#' @details
#' Its accuracy is bounded by that of the derivatives it integrates, which
#' come from [basis_deriv()]: a family with a closed-form derivative gets the
#' accuracy of the quadrature, and one with the finite-difference derivative
#' carries that error into the matrix.
#'
#' All three shipped families and both wrappers register their own method, so
#' this one is reached only by a basis defined outside the package.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param order The derivative order, a single non-negative whole number,
#'   default `0`, or one per variable.
#' @param at,weight Handled in the body of [basis_gram()] before dispatch, so
#'   they never arrive here. Named only because S7 requires a method's formals
#'   to contain the generic's.
#' @param panels The number of equal subintervals, default `50`. For a basis
#'   of several variables the rule is a product and each coordinate gets
#'   `max(2, ceiling(panels^(1/d)))` panels, so 8 apiece at `panels = 50` on
#'   two variables.
#' @param nodes The number of Gauss-Legendre nodes per subinterval, default
#'   `12`.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins.
#'
#' @seealso [numerical_gram()], which it calls; [basis_gram()] for the generic
#'   and the alternative measures.
#'
#' @keywords internal
S7::method(basis_gram, basis) <- function(basis, order = 0L, at = NULL,
                                          weight = NULL, panels = 50L,
                                          nodes = 12L, ...) {
  numerical_gram(basis, order, panels, nodes)
}


#' Gram Matrix by Composite Quadrature
#'
#' @description
#' Computes the inner products of the order-`order` derivatives by
#' Gauss-Legendre on equal subintervals. Shared by [basis_gram.basis()], by
#' [basis_gram.FourierBasis()] when the period is not the interval width, and
#' by [check_basis()] as the independent reference it compares a family's own
#' Gram matrix against.
#'
#' @details
#' # One variable
#'
#' The interval is cut into `panels` equal pieces with an `nodes`-point rule
#' on each, and the matrix is formed as a crossproduct of
#' \eqn{\sqrt{w_i}\,B^{(d)}(t_i)}, which keeps it positive semidefinite
#' whatever the integrand does. It is then symmetrized as `(G + t(G))/2`.
#'
#' # Several variables
#'
#' The rule is a product over the box: the nodes are the product grid of the
#' marginal rules and the weights their products. Each coordinate gets
#' `max(2, ceiling(panels^(1/d)))` panels, so the total node count stays near
#' `panels * nodes^d` and does not grow as `panels^d`. At the defaults on two
#' variables that is 8 panels of 12 nodes per coordinate, 9216 points.
#'
#' # Accuracy
#'
#' Equal panels are not aligned with the knots of a spline, so a family
#' whose derivative has kinks is integrated less accurately than a smooth
#' one; on a polynomial family the result is exact to rounding. A larger
#' `panels` reduces the error where the breaks cannot be aligned.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param order The derivative order, default `0`. Passed through
#'   [check_order()], so a single non-zero order on a basis of several
#'   variables signals an error.
#' @param panels The number of equal subintervals, default `50`, divided among
#'   the coordinates as above for a basis of several variables.
#' @param nodes The number of Gauss-Legendre nodes per subinterval, default
#'   `12`, exact for a polynomial integrand of degree up to 23 on each panel.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins.
#'
#' @seealso [quad_rule()], which builds the composite rule;
#'   [basis_gram()] for the generic.
#'
#' @keywords internal
numerical_gram <- function(basis, order = 0L, panels = 50L, nodes = 12L) {
  d <- basis_nvar(basis)
  order <- check_order(order, d)

  if (d == 1L) {
    r <- quad_rule(seq(basis@lower, basis@upper, length.out = panels + 1L), nodes)
    pts <- r$nodes
    w <- r$weights
  } else {
    # A product rule over the box: the nodes are the lattice of the marginal
    # rules and the weights their products. Kept coarse per coordinate, since
    # the count grows as the power of the number of variables.
    per <- max(2L, as.integer(ceiling(panels^(1 / d))))
    rules <- lapply(seq_len(d), function(j) {
      quad_rule(seq(basis@lower[j], basis@upper[j], length.out = per + 1L), nodes)
    })
    grid <- expand.grid(lapply(rules, function(r) seq_along(r$nodes)))
    pts <- vapply(seq_len(d), function(j) rules[[j]]$nodes[grid[[j]]], numeric(nrow(grid)))
    w <- Reduce(`*`, lapply(seq_len(d), function(j) rules[[j]]$weights[grid[[j]]]))
  }

  b <- basis_deriv(basis, pts, order = order)
  g <- crossprod(sqrt(w) * b)
  g <- (g + t(g)) / 2
  nm <- basis_colnames(basis)
  dimnames(g) <- list(nm, nm)
  g
}


#' @name basis_operator_gram
#' @keywords internal
S7::method(basis_operator_gram, basis) <- function(basis, op, at = NULL,
                                                   weight = NULL, ...) {
  numerical_operator_gram(basis, op, at = at, weight = weight, ...)
}


#' The Roughness Matrix of an Operator by Quadrature
#'
#' @description
#' Integrates \eqn{(Lb)(Lb)^\top} numerically: over Gauss-Legendre panels for
#' the length measure, over the covariate values for the empirical measure,
#' and against a density where one is given. It is the base method of
#' [basis_operator_gram()] and the route a family with no closed form takes.
#'
#' @details
#' The nodes are those of the rule used by [numerical_gram()], and the only
#' difference is what is evaluated at them: \eqn{Lb} from [operator_eval()]
#' instead of one derivative. A piecewise polynomial basis therefore gets an
#' approximation here where its own derivative Gram matrix is exact, with the
#' accuracy of the quadrature.
#'
#' @param basis A [basis] of one variable.
#' @param op A [LinearOperator], with its period resolved.
#' @param at The covariate values for the empirical measure, or `NULL`.
#' @param weight A density to integrate against, or `NULL`.
#' @param panels,nodes The quadrature: how many equal panels the interval is
#'   split into and how many Gauss-Legendre nodes each carries.
#' @param ... Ignored.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns.
#'
#' @seealso [basis_operator_gram()], the generic.
#'
#' @keywords internal
numerical_operator_gram <- function(basis, op, at = NULL, weight = NULL,
                                    panels = 50L, nodes = 12L, ...) {
  check_operator(op)
  if (basis_nvar(basis) > 1L) {
    stop(paste0(
      "a differential operator is not defined for a basis of several",
      " variables:\n  which variable it differentiates is not said. Give",
      " a whole 'order' instead."
    ), call. = FALSE)
  }
  nm <- basis_colnames(basis)

  if (!is.null(at)) {
    if (!is.numeric(at)) stop("'at' must be numeric.", call. = FALSE)
    at <- at[!is.na(at)]
    if (!length(at)) stop("'at' has no usable points.", call. = FALSE)
    lb <- operator_eval(basis, at, op)
    g <- crossprod(lb) / nrow(lb)
  } else {
    breaks <- seq(basis@lower, basis@upper, length.out = panels + 1L)
    r <- quad_rule(breaks, nodes)
    w <- r$weights
    if (!is.null(weight)) {
      if (!is.function(weight)) {
        stop("'weight' must be a function of one numeric vector.",
          call. = FALSE
        )
      }
      wv <- weight(r$nodes)
      if (length(wv) != length(r$nodes) || anyNA(wv) || any(wv < 0)) {
        stop("'weight' must return one non-negative value per point.",
          call. = FALSE
        )
      }
      w <- w * wv
    }
    lb <- operator_eval(basis, r$nodes, op)
    g <- crossprod(sqrt(w) * lb)
  }
  g <- (g + t(g)) / 2
  dimnames(g) <- list(nm, nm)
  g
}


#' Is This the Package's Own Base Class?
#'
#' @description
#' Reports whether an S7 class is the abstract [basis] class. That is how
#' [basis_is_numerical()] tells a method registered on the base class, which
#' is a numerical fallback, from one a subclass supplied.
#'
#' @details
#' Identity is tried first, being the usual case and costing nothing, and the
#' name and package are compared when it fails. The second test is needed
#' because `identical()` can return `FALSE` for a class re-created from the
#' same definition, which happens when the code of a package is evaluated
#' instead of loaded, as coverage tools do.
#'
#' @param cls An S7 class object, as returned by `S7::S7_class()`.
#'
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [route_by_owner()], its only caller, and
#'   [basis_is_numerical()], which relies on it.
#'
#' @keywords internal
is_base_basis_class <- function(cls) {
  if (identical(cls, basis)) return(TRUE)
  identical(attr(cls, "name"), attr(basis, "name")) &&
    identical(attr(cls, "package"), attr(basis, "package"))
}


#' Which of a Basis's Methods Are Numerical
#'
#' @description
#' Reports, for each of [basis_deriv()], [basis_int()] and [basis_gram()],
#' whether the basis supplies its own method or falls back to the numerical
#' one on the abstract [basis] class. `TRUE` means the values come from finite
#' differences or quadrature, so they carry that method's error and cannot be
#' verified against a numerical reference.
#'
#' @details
#' # Use
#'
#' [check_basis()] reads it to decide which of its checks are informative.
#' Comparing a finite difference against a finite difference repeats the same
#' arithmetic, and the two agree however wrong the basis is, so the `deriv`
#' and `integral` checks are not run where this reports `TRUE`.
#' [print.basis()] shows the same information on its `Numerical:` line.
#'
#' # Computation
#'
#' For an ordinary class the result comes from [basis_numerical_route()],
#' whose default method is the owner test: the class on which each method is
#' registered, read through `attr(m, "signature")[[1]]` and tested with
#' [is_base_basis_class()]. A generic with no method at all counts as
#' numerical.
#'
#' Two classes are handled by delegation instead. A [TransformedBasis]
#' reports its parent's: all three of its methods are registered, but each one
#' calls the parent and multiplies, so what is exact about it is whatever was
#' exact about the parent. A [TensorBasis] reports `TRUE` for a quantity that
#' is numerical in **any** margin, every one of its methods being a product of
#' its margins'.
#'
#' # A method that branches
#'
#' The owner test reports where a method came from and not what it does, so
#' a method registered on a concrete class that itself calls the fallback
#' would be reported as exact. [basis_gram.FourierBasis()] does that when the
#' period is not the interval width, and this is why
#' [basis_numerical_route()] is a generic: a family whose route depends on its
#' own parameters registers a method, whose result replaces the owner test.
#' A basis written outside the package can do the same.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#'
#' @return A named logical vector of length 3, with elements `basis_deriv`,
#'   `basis_int` and `basis_gram`, `TRUE` where the numerical fallback is in
#'   force.
#'
#' @seealso [check_basis()], which reads it to decide what to test;
#'   [print.basis()], which prints it; `vignette("defining-a-basis")`, where a
#'   basis is taken from all three `TRUE` to all three `FALSE`.
#'
#' @examples
#' # The shipped families at their default arguments, and wrappers over
#' # them, are exact throughout.
#' basis_is_numerical(bspline_basis(dimension = 5))
#' basis_is_numerical(orthonorm_basis(fourier_basis(dimension = 5)))
#'
#' # A basis defined from its evaluation alone gives TRUE for all three.
#' Bumps <- S7::new_class("Bumps", parent = basis)
#' S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
#'   out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
#'                           "-")^2 / 0.12^2)
#'   colnames(out) <- basis_colnames(basis)
#'   out
#' }
#' bump <- Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)
#' basis_is_numerical(bump)
#'
#' # A product is exact only where every margin is.
#' basis_is_numerical(tensor_basis(bump, poly_basis(dimension = 3)))
#'
#' @export
basis_is_numerical <- function(basis) {
  # A transformed basis registers all three methods, but each of them delegates
  # to the parent and multiplies, so what is numerical about it is whatever was
  # numerical about the parent. Reporting its own methods would say a value was
  # exact when it was a finite difference arranged as a matrix.
  if (S7::S7_inherits(basis, TransformedBasis)) {
    return(basis_is_numerical(basis@parent_basis))
  }
  # A tensor product is exact exactly when its marginals are: every one of its
  # methods is a product of theirs.
  if (S7::S7_inherits(basis, TensorBasis)) {
    flags <- vapply(basis@marginals, basis_is_numerical, logical(3L))
    return(apply(matrix(flags, nrow = 3L), 1L, any) |>
      stats::setNames(c("basis_deriv", "basis_int", "basis_gram")))
  }
  basis_numerical_route(basis)
}


#' Which Route a Basis's Methods Take
#'
#' @description
#' Reports, for each of [basis_deriv()], [basis_int()] and [basis_gram()],
#' whether this basis computes the quantity numerically. [basis_is_numerical()]
#' calls it once it has handled the two wrapper classes, and a family
#' overrides it when its route depends on its own parameters instead of the
#' class on which its method is registered.
#'
#' @details
#' # Purpose of the generic
#'
#' The default method reads the class on which each method is registered,
#' which records where a method came from and not what it does. A method
#' registered on a concrete class may still call the fallback:
#' [basis_gram.FourierBasis()] does, whenever the period is not the interval
#' width, and the owner is `FourierBasis` in both branches. Such a family
#' reports its own route by registering a method here. The generic is
#' exported so that a basis written outside the package can do the same.
#'
#' The two wrapper classes, [TransformedBasis] and [TensorBasis], are handled
#' by [basis_is_numerical()] before this generic is reached; called directly
#' on a wrapper, the default method reports the wrapper's own methods.
#'
#' # Writing one
#'
#' A method takes the default through [S7::super()] and sets the entries
#' decided by its own branching:
#'
#' ```r
#' S7::method(basis_numerical_route, MyBasis) <- function(basis, ...) {
#'   out <- basis_numerical_route(S7::super(basis, basis7::basis))
#'   out[["basis_gram"]] <- !my_closed_form_applies(basis)
#'   out
#' }
#' ```
#'
#' The package's own override calls [route_by_owner()] instead, which is that
#' default under a name: inside a method the formal `basis` shadows the class
#' of the same name, so reaching the class through `S7::super()` would require
#' naming the package inside itself.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param ... Unused, and accepted so a method's signature can match.
#'
#' @return A named logical vector of length 3, with elements `basis_deriv`,
#'   `basis_int` and `basis_gram`, `TRUE` where the quantity is computed
#'   numerically.
#'
#' @seealso [basis_is_numerical()], the predicate consumers call;
#'   [route_by_owner()] for the default computation; and
#'   [basis_numerical_route.FourierBasis()] for the one family that overrides.
#'
#' @examples
#' # A B-spline basis uses the default method.
#' basis_numerical_route(bspline_basis(dimension = 5))
#'
#' # A Fourier basis whose period is not the interval width computes its Gram
#' # matrix by quadrature, which its own method reports and the owner test
#' # would not.
#' basis_numerical_route(fourier_basis(dimension = 5, omega = 0.7))
#'
#' @export
basis_numerical_route <- S7::new_generic(
  "basis_numerical_route", "basis",
  function(basis, ...) S7::S7_dispatch()
)


#' @title Default Numerical Route by Owner Test
#' @name basis_numerical_route.basis
#' @description
#' Returns [route_by_owner()]: the class on which each of [basis_deriv()],
#' [basis_int()] and [basis_gram()] is registered. That is correct for every
#' family whose methods take one route, which excludes the Fourier family,
#' and for a class that is not a wrapper; [basis_is_numerical()] handles the
#' wrappers itself.
#' @param basis A basis object.
#' @param ... Unused, and accepted so the signature matches the generic's.
#' @return The named logical vector [basis_numerical_route()] describes.
#' @seealso [basis_numerical_route.FourierBasis()], the one family that needs
#'   more than this.
#' @keywords internal
S7::method(basis_numerical_route, basis) <- function(basis, ...) {
  route_by_owner(basis)
}


#' The Owner Test Behind the Default Route
#'
#' @description
#' Reports which of [basis_deriv()], [basis_int()] and [basis_gram()] this
#' basis takes from the numerical fallback, by reading the class on which
#' each method is registered through `attr(m, "signature")[[1]]` and testing
#' it with
#' [is_base_basis_class()]. A generic with no method at all counts as
#' numerical, the base class being where the fallback lives.
#'
#' @details
#' It is the body of [basis_numerical_route.basis()], kept under a name so that
#' the package's own override can start from it. Inside a method the formal
#' `basis` shadows the class of the same name, so reaching that default
#' through `S7::super()` would require naming the package inside itself; a
#' basis written elsewhere uses `S7::super()`.
#'
#' The test does not see which branch a method takes once it has been
#' reached; [basis_numerical_route()] lets a family report that itself.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#'
#' @return A named logical vector of length 3, elements `basis_deriv`,
#'   `basis_int` and `basis_gram`.
#'
#' @seealso [basis_numerical_route()], the generic it is the default of, and
#'   [is_base_basis_class()], which it tests each owner with.
#'
#' @keywords internal
route_by_owner <- function(basis) {
  cls <- S7::S7_class(basis)
  gens <- list(
    basis_deriv = basis_deriv,
    basis_int = basis_int,
    basis_gram = basis_gram
  )
  vapply(gens, function(g) {
    m <- tryCatch(S7::method(g, cls), error = function(e) NULL)
    if (is.null(m)) return(TRUE)
    owner <- attr(m, "signature")[[1L]]
    is_base_basis_class(owner)
  }, logical(1))
}
