#' @include generics.R
NULL




#' Numerically Differentiate a Matrix-Valued Function
#'
#' @description
#' Estimates the `order`-th derivative of `f` at each point of `x` by a single
#' finite-difference stencil, symmetric where the interval leaves room for it
#' and one-sided where it does not. One stencil of the order wanted, never a
#' composition of lower-order differences, so the error is the truncation of
#' that one formula. The engine behind every numerical fallback in the
#' package.
#'
#' @details
#' # One stencil per point, chosen by the room available
#'
#' A basis is evaluated at its endpoints as readily as anywhere else, and a
#' symmetric stencil centered on an endpoint would ask for points outside the
#' interval, where the basis throws. Each point therefore gets the central
#' stencil when both sides have room, and otherwise the one-sided stencil that
#' points inward, with the same order and the same number of nodes. All three
#' weight vectors come from [numericals7::fd_weights()] on the offsets
#' [numericals7::fd_offsets()] supplies.
#'
#' Accuracy is not lost at the ends. Measured against exact Legendre
#' derivatives on \eqn{[0, 1]}, the first derivative agrees to 5.1e-10
#' relative at the two endpoints and to 3.9e-10 in the interior.
#'
#' # The step
#'
#' The step is [numericals7::fd_step()]'s,
#' \eqn{\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)}, which balances
#' truncation against rounding for order \eqn{d}. It is written nowhere in
#' this package, a second copy of a rule with one home being how two packages
#' come to disagree.
#'
#' It is then capped at `0.4 * (upper - lower) / (2 * reach)` so that the
#' whole stencil fits inside the interval: `0.2` of the width at orders 1 and
#' 2, where the stencil has three nodes, and `0.1` at orders 3 and 4, where it
#' has five.
#'
#' # What it costs in accuracy
#'
#' Measured against exact Legendre derivatives at 21 interior points of
#' \eqn{[0, 1]}, relative to the scale of the answer: 3.9e-10 at order 1,
#' 1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08 at order 4. That is what
#' a family which registers no derivative method gets, and the reason to write
#' a closed form where one exists.
#'
#' @param f A function of one numeric vector returning a numeric matrix with
#'   one row per element. Called `2 * reach + 1` times, so a costly `f` is
#'   evaluated three or five times over.
#' @param x A numeric vector of evaluation points. `NA` entries are given the
#'   central stencil and propagate to an `NA` row.
#' @param order The derivative order, a single positive whole number.
#' @param lower,upper The endpoints of the interval `f` is defined on, used to
#'   choose each point's stencil and to cap the step.
#' @param step_scale A multiplier on the step, default `1`. [fd_reference()]
#'   passes `0.5` to measure the reference's own uncertainty by how much the
#'   answer moves.
#'
#' @return A numeric matrix with `length(x)` rows and as many columns as `f`
#'   returns, with no dimnames; callers add them through [name_columns()].
#'
#' @seealso [numericals7::fd_weights()], [numericals7::fd_offsets()] and
#'   [numericals7::fd_step()], which supply the three pieces;
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

  w <- list(
    c = numericals7::fd_weights(off$central, order),
    f = numericals7::fd_weights(off$forward, order),
    b = numericals7::fd_weights(off$backward, order)
  )
  offs <- list(c = off$central, f = off$forward, b = off$backward)

  n <- length(x)
  nnode <- 2L * reach + 1L
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
#' They are computed at call time, which keeps every node count available.
#' That is what the exact spline rules need: one rule per knot interval, sized
#' from the degree and the derivative order, where a table would offer only
#' the counts someone thought to tabulate.
#'
#' @param n The number of nodes, a single positive whole number. `1` returns
#'   the midpoint rule directly. `n < 1` throws
#'   `'n' must be a positive integer.`
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
  n <- as.integer(n)
  if (n < 1L) stop("'n' must be a positive integer.", call. = FALSE)
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
#' Nothing is validated. `breaks` must be increasing and of length at least
#' two; the nodes come out in interval order, which is increasing when the
#' breaks are.
#'
#' @param breaks A numeric vector of at least two increasing breakpoints. The
#'   first and last are the ends of the span.
#' @param n The number of nodes per interval, a single positive whole number.
#'
#' @return A list of two numeric vectors of length
#'   `n * (length(breaks) - 1)`: `nodes` and `weights`, ordered interval by
#'   interval.
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
#' one finite-difference stencil of the order asked for, applied to
#' [basis_eval()]. A subclass supplying its evaluation alone therefore answers
#' [basis_deriv()] immediately, and registering a closed form later takes over
#' through dispatch with no change to calling code.
#'
#' @details
#' # Accuracy
#'
#' One stencil of the order wanted, never a chain of first differences.
#' Measured against exact Legendre derivatives at 21 interior points of
#' \eqn{[0, 1]}, relative to the scale of the answer: 3.9e-10 at order 1,
#' 1.5e-08 at order 2, 3.3e-09 at order 3 and 5.6e-08 at order 4. See
#' [numerical_deriv_matrix()] for the stencil, the step and the endpoint rule.
#'
#' # A basis of several variables
#'
#' The stencil differentiates along one coordinate at a time, replacing that
#' column of the points and holding the others. A mixed partial such as
#' `c(1, 1)` therefore throws: a stencil in the plane carries the product of
#' two errors, and the one family that needs mixed partials,
#' [tensor_basis()], computes them exactly from its margins.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points inside the basis interval: a numeric vector for
#'   a basis of one variable, or a matrix of [basis_nvar()] columns for a
#'   basis of several.
#' @param order The derivative order, a single non-negative whole number, or
#'   one per variable. More than one non-zero entry throws.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names [basis_colnames()].
#'
#' @seealso [numerical_deriv_matrix()], which does the work;
#'   [basis_is_numerical()] to ask an object whether its derivatives come from
#'   here; [basis_deriv()] for the generic.
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
#' the sorted evaluation points, so the whole set costs one pass where a
#' quadrature each would cost as many.
#'
#' @details
#' # How it is accumulated
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
#' `2 * nodes - 1` exactly on each segment, so on a polynomial family the
#' result is exact to rounding: measured against the closed-form Legendre
#' integrals at 21 points, 3.3e-16. On a family that is not polynomial the
#' error is the rule's own on each segment, which shrinks with the spacing of
#' the evaluation points; a single distant point is integrated by one rule
#' over the whole span.
#'
#' # One variable only
#'
#' A basis of several variables throws. The integral there is over a box, one
#' iterated integral per variable, and a family that wants it supplies its
#' own, as [basis_int.TensorBasis()] does from its margins.
#'
#' @param basis A basis object of one variable, of any class inheriting from
#'   [basis]. More than one variable throws.
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
#' [numerical_gram()], which is where the work is and which the other families
#' also call when their closed form does not apply.
#'
#' @details
#' Its accuracy is bounded by the derivative it integrates. On a polynomial
#' family the order-0 matrix agrees with the closed form to 5.2e-15, while at
#' order 2 the gap is 3.0e-08 relative, which is the finite-difference error
#' of [basis_deriv.basis()] carried through the integral, and no fault of the
#' quadrature.
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
#'   `ceiling(panels^(1/d))` panels, so 8 apiece at `panels = 50` on two
#'   variables.
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
#' whatever the integrand does. It is then symmetrized as `(G + t(G))/2`, the
#' two triangles of a crossproduct differing in their last bits.
#'
#' # Several variables
#'
#' The rule is a product over the box: the nodes are the lattice of the
#' marginal rules and the weights their products. Each coordinate gets
#' `max(2, ceiling(panels^(1/d)))` panels, so the total node count stays near
#' `panels * nodes^d` and does not grow as `panels^d`. At the defaults on two
#' variables that is 8 panels of 12 nodes per coordinate, 9216 points.
#'
#' # Accuracy
#'
#' Equal panels line up with nothing in particular, so a family whose
#' derivative has kinks is integrated less well than one whose does not. On a
#' polynomial family the order-0 matrix agrees with the closed form to
#' 5.2e-15; on a cubic B-spline at order 2, where the second derivative kinks
#' at knots the panels miss, the same comparison is 1.4e-03, or 2.1e-06
#' relative. Raising `panels` is the remedy where the breaks cannot be aligned.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param order The derivative order, default `0`. Passed through
#'   [check_order()], so a single non-zero order on a basis of several
#'   variables throws.
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


#' Is This the Package's Own Base Class?
#'
#' @description
#' Reports whether an S7 class is the abstract [basis] class. That is how
#' [basis_is_numerical()] tells a method registered on the base class, which
#' is a numerical fallback, from one a subclass supplied.
#'
#' @details
#' Identity is tried first, being the usual case and costing nothing, and the
#' name and package are compared when it fails. The second test is necessary:
#' `identical()` on an S7 class is object identity, so it is `FALSE` for a
#' class re-created from the same definition, and that is what happens
#' whenever a package's code is evaluated instead of loaded. Coverage tools do exactly
#' that, so an identity-only test passes every ordinary check and fails under
#' `covr` alone.
#'
#' The same defect in `linkfunctions7` made every fallback differentiate the
#' order below it, and the log link's fourth derivative came back wrong by a
#' factor of 900 while the whole five-platform check matrix stayed green.
#'
#' @param cls An S7 class object, as returned by `S7::S7_class()`.
#'
#' @return A single `TRUE` or `FALSE`.
#'
#' @seealso [basis_is_numerical()], its only caller.
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
#' # What it decides
#'
#' [check_basis()] reads it to decide which of its checks can mean anything.
#' Comparing a finite difference against a finite difference is the same
#' arithmetic twice, agreeing however wrong the basis is, so the `deriv` and
#' `integral` checks are not run at all where this reports `TRUE`.
#' [print.basis()] shows the same information on its `Numerical:` line.
#'
#' # How it is answered
#'
#' For an ordinary class, by asking which class each method is registered on
#' through `attr(m, "signature")[[1]]` and testing it with
#' [is_base_basis_class()]. A generic with no method at all counts as
#' numerical.
#'
#' Two classes are answered by delegation instead. A [TransformedBasis]
#' reports its parent's: all three of its methods are registered, but each one
#' calls the parent and multiplies, so what is exact about it is whatever was
#' exact about the parent. A [TensorBasis] reports `TRUE` for a quantity that
#' is numerical in **any** margin, every one of its methods being a product of
#' its margins'.
#'
#' # What it cannot see
#'
#' The question asked is which class the method is registered on. It says
#' nothing about which branch that method takes. A method registered on a concrete class
#' that itself calls the fallback is reported as exact:
#' [basis_gram.FourierBasis()] does that when the period is not the interval
#' width, and `basis_is_numerical()` reports `FALSE` for its Gram matrix
#' either way.
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
#' # Every shipped family, and every wrapper over one, is exact throughout.
#' basis_is_numerical(bspline_basis(dimension = 5))
#' basis_is_numerical(orthonorm_basis(fourier_basis(dimension = 5)))
#'
#' # A basis defined from its evaluation alone answers TRUE to all three.
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
