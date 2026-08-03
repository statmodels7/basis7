#' @include generics.R
NULL


#' Finite-Difference Weights for an Arbitrary Stencil
#'
#' @description
#' The weights that turn function values at \code{x + offsets * h} into the
#' \code{order}-th derivative at \code{x}, divided by \code{h^order}.
#'
#' @details
#' The weights solve the linear system that makes the combination exact on
#' polynomials up to the degree the stencil can carry. Writing \eqn{s_j} for
#' the offsets, the Taylor expansion of \eqn{\sum_j w_j f(x + s_j h)} has
#' \eqn{f^{(i)}(x)} multiplied by \eqn{h^i/i! \sum_j w_j s_j^i}, so requiring
#' \eqn{\sum_j w_j s_j^i = 0} for \eqn{i \neq d} and \eqn{= d!} for
#' \eqn{i = d} leaves exactly \eqn{h^d f^{(d)}(x)}. That is a Vandermonde
#' system in the offsets, solved once per stencil shape.
#'
#' Building the weights this way, rather than composing lower-order
#' differences, is what keeps a high order usable: each numerical
#' differentiation multiplies the error of the one before it, so a fourth
#' derivative reached by four nested first differences is noise. One stencil,
#' never nested.
#'
#' @param offsets A numeric vector of stencil offsets, in units of the step.
#' @param order The derivative order.
#'
#' @return A numeric vector of weights, the same length as \code{offsets}.
#'
#' @keywords internal
fd_weights <- function(offsets, order) {
  n <- length(offsets)
  # Row i is the offsets raised to the power i - 1: the exponent indexes the
  # rows, the offsets the columns.
  A <- outer(seq_len(n) - 1L, offsets, function(power, s) s^power)
  rhs <- numeric(n)
  rhs[order + 1L] <- factorial(order)
  solve(A, rhs)
}


#' Stencil Offsets for a Derivative Order
#'
#' @description
#' The symmetric offsets used away from the interval endpoints, and the
#' one-sided ones used where a symmetric stencil would not fit.
#'
#' @param order The derivative order.
#'
#' @return A list with the reach and the three offset vectors.
#'
#' @keywords internal
fd_offsets <- function(order) {
  r <- max(1L, as.integer(ceiling(order / 2)))
  list(
    reach = r,
    central = seq.int(-r, r),
    forward = seq.int(0, 2 * r),
    backward = seq.int(-2 * r, 0)
  )
}


#' Numerically Differentiate a Matrix-Valued Function
#'
#' @description
#' The \code{order}-th derivative of \code{f} at each point of \code{x}, by a
#' single finite-difference stencil, symmetric where the interval leaves room
#' for it and one-sided at the endpoints.
#'
#' @details
#' A basis is evaluated at its endpoints as readily as anywhere else, and a
#' symmetric stencil centred on an endpoint would ask for points outside the
#' interval, where the basis is not defined. Those points therefore get a
#' one-sided stencil of the same order and the same number of nodes, built by
#' \code{\link{fd_weights}} from shifted offsets.
#'
#' The step is \eqn{\varepsilon^{1/(d+2)}\max(1, \lvert x\rvert)}, which
#' balances truncation against rounding for order \eqn{d}, capped so that the
#' whole stencil fits inside the interval.
#'
#' @param f A function of a numeric vector returning a matrix with one row per
#'   element.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param lower,upper The endpoints of the interval \code{f} is defined on.
#' @param step_scale A factor applied to the step. Halving it is how
#'   \code{\link{fd_reference}} measures its own uncertainty.
#'
#' @return A numeric matrix with \code{length(x)} rows.
#'
#' @keywords internal
numerical_deriv_matrix <- function(f, x, order, lower, upper, step_scale = 1) {
  off <- fd_offsets(order)
  reach <- off$reach

  h <- .Machine$double.eps^(1 / (order + 2)) * pmax(1, abs(x)) * step_scale
  h <- pmin(h, 0.4 * (upper - lower) / (2 * reach))

  # Which stencil each point can afford: symmetric when both sides have room,
  # otherwise the one-sided stencil that points into the interval.
  room_left <- (x - lower) >= reach * h
  room_right <- (upper - x) >= reach * h
  kind <- ifelse(room_left & room_right, "c", ifelse(room_right, "f", "b"))
  kind[is.na(x)] <- "c"

  w <- list(
    c = fd_weights(off$central, order),
    f = fd_weights(off$forward, order),
    b = fd_weights(off$backward, order)
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
#' The \code{n}-point Gauss-Legendre rule on \eqn{[-1, 1]}, which integrates
#' polynomials of degree up to \eqn{2n - 1} exactly.
#'
#' @details
#' The nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix of
#' the Legendre recurrence and the weights come from the first component of
#' each eigenvector, which is the Golub-Welsch construction. Computing them
#' rather than tabulating them keeps any node count available, which the exact
#' spline rules need: a rule per knot interval, sized from the degree.
#'
#' @param n The number of nodes, a positive integer.
#'
#' @return A list with components \code{nodes} and \code{weights}.
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
#' Places an \code{n}-point Gauss-Legendre rule on each interval given by
#' consecutive breakpoints, and returns the pooled nodes and weights.
#'
#' @param breaks A numeric vector of at least two increasing breakpoints.
#' @param n The number of nodes per interval.
#'
#' @return A list with components \code{nodes} and \code{weights}.
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
#' @description
#' The default derivative method, applying one finite-difference stencil to
#' \code{\link{basis_eval}}. It is what a basis that registers no derivative
#' method of its own uses, so that implementing the evaluation is enough to
#' have a complete basis.
#' @details
#' See \code{\link{numerical_deriv_matrix}} for the stencil and the step, and
#' \code{\link{basis_is_numerical}} for asking an object whether its
#' derivatives come from here.
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
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
#' @description
#' The default integration method: composite Gauss-Legendre from the lower
#' endpoint, accumulated over the sorted evaluation points so that the whole
#' set costs one pass rather than one quadrature each.
#' @details
#' The rule is placed on the segments between consecutive evaluation points,
#' and the results are accumulated, which makes the value at the lower
#' endpoint exactly zero by construction rather than by cancellation.
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points.
#' @param nodes The number of quadrature nodes per segment.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
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
#' @description
#' The default inner-product method: composite Gauss-Legendre over the basis
#' interval, applied to the outer product of the requested derivatives.
#' @details
#' The rule is placed on \code{panels} equal subintervals. A basis with
#' closed-form inner products, or one that is piecewise polynomial and so can
#' be integrated exactly by a rule sized from its degree, registers its own
#' method and this one is not used.
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order.
#' @param panels The number of subintervals.
#' @param nodes The number of quadrature nodes per subinterval.
#' @param ... Unused.
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
#' @keywords internal
S7::method(basis_gram, basis) <- function(basis, order = 0L, at = NULL,
                                          weight = NULL, panels = 50L,
                                          nodes = 12L, ...) {
  numerical_gram(basis, order, panels, nodes)
}


#' Gram Matrix by Composite Quadrature
#'
#' @description
#' The inner products of the order-\code{order} derivatives, by Gauss-Legendre
#' on equal subintervals. Shared by the default method and by any basis whose
#' closed form does not apply to the arguments it was given.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param order The derivative order.
#' @param panels The number of subintervals.
#' @param nodes The number of quadrature nodes per subinterval.
#'
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
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
#' Asks whether an S7 class is the abstract \code{\link{basis}} class, which is
#' how a method registered on it is told apart from one a subclass supplied.
#'
#' @details
#' Identity is tried first because it is the usual case and costs nothing, then
#' the name and package, because identity is not preserved when a package's
#' code is re-evaluated rather than loaded. Coverage tools do exactly that, so
#' an identity-only test passes every ordinary check and fails there.
#'
#' @param cls An S7 class.
#'
#' @return \code{TRUE} or \code{FALSE}.
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
#' Reports, for each of the three derived generics, whether the basis supplies
#' its own method or falls back to the numerical one.
#'
#' @details
#' A value computed by finite differences or by quadrature cannot be checked
#' against a numerical reference: the comparison would be the same arithmetic
#' twice, agreeing however wrong the basis is. \code{\link{check_basis}} uses
#' this to report such an order as not checked rather than as passed, which is
#' the difference between a validator and a formality.
#'
#' @param basis An object inheriting from class \code{basis}.
#'
#' @return A named logical vector with elements \code{basis_deriv},
#'   \code{basis_int} and \code{basis_gram}, \code{TRUE} where the numerical
#'   fallback is in force.
#'
#' @examples
#' basis_is_numerical(bspline_basis(dimension = 5))
#'
#' @export
basis_is_numerical <- function(basis) {
  # A transformed basis registers all three methods, but each of them delegates
  # to the parent and multiplies, so what is numerical about it is whatever was
  # numerical about the parent. Reporting its own methods would say a value was
  # exact when it was a finite difference wearing a matrix.
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
