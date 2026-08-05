#' @include numerical_fallbacks.R
NULL


#' Linearly Transformed Basis
#'
#' @description
#' A basis obtained from another by a fixed linear map of its functions,
#' \eqn{\tilde{B}(x) = B(x)\,T}. Constructed by \code{\link{orthonorm_basis}},
#' \code{\link{constrain_basis}} or \code{\link{dr_basis}}.
#'
#' @details
#' Orthonormalizing a basis, restricting it to satisfy a linear constraint, and
#' rebuilding it so that it diagonalizes an inner product are the same
#' operation with different matrices, so they share one class. Derivatives and
#' integrals transform by the same \eqn{T}, because differentiation and
#' integration are linear and \eqn{T} does not depend on \eqn{x}; the Gram
#' matrix transforms by congruence, \eqn{T^\top G\, T}, so a parent with an
#' exact Gram matrix passes its exactness on.
#'
#' \eqn{T} may have fewer columns than rows, which is how a constraint reduces
#' the dimension.
#'
#' Transforms compose by multiplication rather than by nesting: transforming a
#' \code{TransformedBasis} again produces one object holding the product of the
#' two matrices, so a chain of transforms costs one matrix multiplication per
#' evaluation however long it is.
#'
#' @inheritParams basis
#' @param parent_basis The basis being transformed.
#' @param transform The matrix \eqn{T}, with one row per parent function.
#'
#' @return An object of class \code{TransformedBasis}.
#'
#' @seealso \code{\link{orthonorm_basis}}, \code{\link{constrain_basis}},
#'   \code{\link{dr_basis}}
#'
#' @examples
#' o <- orthonorm_basis(bspline_basis(dimension = 6))
#' S7::S7_inherits(o, TransformedBasis)
#' dim(o@transform)
#'
#' @export
TransformedBasis <- S7::new_class(
  "TransformedBasis",
  parent = basis,
  properties = list(
    parent_basis = basis,
    transform = S7::class_numeric
  ),
  validator = function(self) {
    tm <- self@transform
    if (!is.matrix(tm)) return("@transform must be a matrix")
    if (nrow(tm) != self@parent_basis@dimension) {
      return("@transform must have one row per parent basis function")
    }
    if (ncol(tm) != self@dimension) {
      return("@transform must have one column per function of the result")
    }
    NULL
  }
)


#' Build a Transformed Basis
#'
#' @description
#' Wraps a basis in a \code{\link{TransformedBasis}}, collapsing the transform
#' into the parent's when the parent is already one.
#'
#' @param basis The basis to transform.
#' @param transform The matrix \eqn{T}.
#' @param name The name of the resulting basis.
#' @param prefix The prefix for its column names.
#' @param params Extra entries for \code{basis_params}.
#'
#' @return An object of class \code{\link{TransformedBasis}}.
#'
#' @keywords internal
new_transformed <- function(basis, transform, name, prefix, params = list()) {
  if (S7::S7_inherits(basis, TransformedBasis)) {
    transform <- basis@transform %*% transform
    basis <- basis@parent_basis
  }
  TransformedBasis(
    basis_name = name,
    dimension = ncol(transform),
    lower = basis@lower,
    upper = basis@upper,
    basis_params = c(list(prefix = prefix), params),
    parent_basis = basis,
    transform = transform
  )
}


#' Column Names of a Transformed Basis
#'
#' @name basis_colnames.TransformedBasis
#' @description
#' The transform's prefix numbered from one. The parent's names cannot be kept:
#' each new function is a combination of all of them, and a transform may
#' produce fewer functions than it consumed.
#' @param basis A \code{\link{TransformedBasis}} object.
#' @param ... Unused.
#' @return A character vector of length \code{basis@dimension}.
#' @keywords internal
S7::method(basis_colnames, TransformedBasis) <- function(basis, ...) {
  paste0(basis@basis_params$prefix, seq_len(basis@dimension))
}


#' Evaluate a Transformed Basis
#'
#' @name basis_eval.TransformedBasis
#' @description The parent's evaluation, multiplied by the transform.
#' @param basis A \code{\link{TransformedBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows.
#' @keywords internal
S7::method(basis_eval, TransformedBasis) <- function(basis, x, ...) {
  name_columns(basis_eval(basis@parent_basis, x) %*% basis@transform, basis)
}


#' Derivatives of a Transformed Basis
#'
#' @name basis_deriv.TransformedBasis
#' @description
#' The parent's derivatives, multiplied by the transform: differentiation is
#' linear and the transform does not depend on \code{x}.
#' @param basis A \code{\link{TransformedBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows.
#' @keywords internal
S7::method(basis_deriv, TransformedBasis) <- function(basis, x, order = 1L, ...) {
  d <- basis_deriv(basis@parent_basis, x, order = order)
  name_columns(d %*% basis@transform, basis)
}


#' Integral of a Transformed Basis
#'
#' @name basis_int.TransformedBasis
#' @description
#' The parent's integral, multiplied by the transform. The anchor at the lower
#' endpoint survives, since a linear combination of columns that are all zero
#' there is zero there.
#' @param basis A \code{\link{TransformedBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows.
#' @keywords internal
S7::method(basis_int, TransformedBasis) <- function(basis, x, ...) {
  name_columns(basis_int(basis@parent_basis, x) %*% basis@transform, basis)
}


#' Gram Matrix of a Transformed Basis
#'
#' @name basis_gram.TransformedBasis
#' @description
#' The congruence \eqn{T^\top G\, T} of the parent's Gram matrix, so a parent
#' whose inner products are exact passes that on rather than falling back to
#' quadrature.
#' @param basis A \code{\link{TransformedBasis}} object.
#' @param order The derivative order.
#' @param at,weight Handled by the generic before dispatch; unused here.
#' @param ... Passed to the parent's method.
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
#' @keywords internal
S7::method(basis_gram, TransformedBasis) <- function(basis, order = 0L,
                                                     at = NULL, weight = NULL,
                                                     ...) {
  g <- basis_gram(basis@parent_basis, order = order, ...)
  out <- crossprod(basis@transform, g %*% basis@transform)
  out <- (out + t(out)) / 2
  nm <- basis_colnames(basis)
  dimnames(out) <- list(nm, nm)
  out
}


#' Cholesky Factorization, With the Rank Decided Before It
#'
#' @description
#' The Cholesky factor of a symmetric matrix, or \code{NULL} when the matrix is
#' not positive definite to the given relative tolerance.
#'
#' @details
#' The verdict comes from the eigenvalues rather than from whether
#' \code{\link[base]{chol}} raises. On a matrix with an exactly zero eigenvalue
#' the pivot that should be zero comes out positive or negative according to
#' rounding, so \code{chol()} succeeds on some platforms and fails on others,
#' and a construction that asks it whether a penalty is usable gets a different
#' answer on different machines. Comparing the smallest eigenvalue with the
#' largest is a statement about the matrix, and gives the same answer
#' everywhere.
#'
#' @param m A symmetric numeric matrix.
#' @param tol The relative tolerance below which the smallest eigenvalue counts
#'   as zero.
#'
#' @return The upper triangular Cholesky factor, or \code{NULL}.
#'
#' @keywords internal
chol_pd <- function(m, tol = 1e-12) {
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  if (!length(ev) || anyNA(ev)) return(NULL)
  if (max(ev) <= 0 || min(ev) <= tol * max(ev)) return(NULL)
  tryCatch(chol(m), error = function(e) NULL)
}


#' Orthonormalize a Basis
#'
#' @description
#' Returns the basis whose functions span the same space and are orthonormal in
#' \eqn{L^2} on the basis interval.
#'
#' @details
#' The transform is read off the Gram matrix rather than estimated on a grid.
#' Writing \eqn{G = R^\top R} for the Cholesky factorization of the Gram matrix,
#' the basis \eqn{B R^{-1}} has Gram matrix
#' \eqn{R^{-\top} R^\top R R^{-1} = I}. Because the parent's Gram matrix is
#' exact for the families that ship with the package, so is the
#' orthonormalization: there is no grid, no number of points to choose, and no
#' scale factor to correct.
#'
#' Orthonormalizing an already orthonormal basis returns it unchanged, up to
#' rounding, and the transforms collapse rather than nesting.
#'
#' @param basis The basis to orthonormalize.
#' @param order The derivative order whose inner products are made the
#'   identity. Zero, the default, orthonormalizes the functions themselves.
#'
#' @return An object of class \code{\link{TransformedBasis}}.
#'
#' @seealso \code{\link{basis_gram}}, \code{\link{constrain_basis}}
#'
#' @examples
#' o <- orthonorm_basis(bspline_basis(dimension = 6))
#' round(basis_gram(o), 12) # the identity, exactly
#'
#' @export
orthonorm_basis <- function(basis, order = 0L) {
  order <- check_order(order)
  g <- basis_gram(basis, order = order)
  r <- chol_pd(g)
  if (is.null(r)) {
    stop(
      "The Gram matrix is singular, so the basis functions are linearly ",
      "dependent and cannot be orthonormalized. Reduce 'dimension', or ",
      "orthonormalize at order 0.",
      call. = FALSE
    )
  }
  tm <- backsolve(r, diag(nrow(r)))
  new_transformed(basis, tm, paste0("orthonorm(", basis@basis_name, ")"), "on")
}


#' Restrict a Basis to a Linear Constraint
#'
#' @description
#' Returns the basis whose coefficient vectors are exactly those satisfying
#' \eqn{C\beta = 0}, with the dimension reduced by the rank of \eqn{C}.
#'
#' @details
#' The transform is an orthonormal basis of the null space of \eqn{C},
#' extracted from its singular value decomposition, so the constrained basis
#' spans precisely the admissible functions and no reparametrization of the
#' constraint changes the space it produces.
#'
#' What the package supplies is the mechanics. Which constraint a model term
#' should carry, whether a sum-to-zero condition for identifiability or
#' orthogonality to a linear part, is a decision that needs to know what the
#' term means, and belongs to the layer that does.
#'
#' @param basis The basis to restrict.
#' @param constraint A matrix with one column per basis function, or a vector
#'   for a single constraint.
#' @param tol The relative tolerance below which a singular value counts as
#'   zero when determining the rank.
#'
#' @return An object of class \code{\link{TransformedBasis}}.
#'
#' @seealso \code{\link{dr_basis}}, \code{\link{orthonorm_basis}}
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#'
#' # sum to zero over a grid: the usual identifiability constraint
#' x <- seq(0, 1, length.out = 200)
#' cs <- constrain_basis(b, colSums(basis_eval(b, x)))
#' cs@dimension
#' max(abs(colSums(basis_eval(cs, x))))
#'
#' @export
constrain_basis <- function(basis, constraint, tol = 1e-10) {
  cm <- if (is.matrix(constraint)) constraint else matrix(constraint, nrow = 1L)
  if (!is.numeric(cm) || ncol(cm) != basis@dimension) {
    stop(sprintf(
      "'constraint' must be numeric with %d columns, one per basis function.",
      basis@dimension
    ), call. = FALSE)
  }
  if (anyNA(cm)) stop("'constraint' must not contain missing values.", call. = FALSE)

  s <- svd(cm, nu = 0L, nv = basis@dimension)
  rank <- sum(s$d > tol * max(s$d, 0))
  if (rank >= basis@dimension) {
    stop(
      "The constraint leaves no functions: its rank equals the dimension of ",
      "the basis.",
      call. = FALSE
    )
  }
  tm <- s$v[, seq.int(rank + 1L, basis@dimension), drop = FALSE]
  new_transformed(
    basis, tm, paste0("constrained(", basis@basis_name, ")"), "cn",
    params = list(constraint_rank = rank)
  )
}


#' Demmler-Reinsch Basis
#'
#' @description
#' Returns the basis that simultaneously diagonalizes the empirical inner
#' product at the given points and the penalty, and is empirically orthogonal
#' to the constraint functions.
#'
#' @details
#' The construction has three steps. The constraint matrix
#' \eqn{C = (\mathbf{1}, x)^\top B} is formed and the basis restricted to its
#' null space \eqn{V_0}, which makes the remaining functions empirically
#' orthogonal to a constant and to \eqn{x}. The pencil
#' \eqn{(V_0^\top (B^\top B/n) V_0,\; V_0^\top P V_0)} is then diagonalized, and
#' the transform is \eqn{T = V_0 A}. The resulting design matrix
#' \eqn{Z = B T} has \eqn{Z^\top Z} diagonal, has \eqn{T^\top P T} equal to the
#' identity, and satisfies \eqn{(\mathbf{1}, x)^\top Z = 0} exactly.
#'
#' Three properties make this construction usable where others are not. It
#' factorizes only a \eqn{q \times K} and a \eqn{(K-q) \times (K-q)} matrix,
#' never anything of the size of the sample. It tolerates a rank-deficient
#' \eqn{B}, which equally spaced knots produce whenever the data leave a knot
#' span empty, because the matrix inverted is the penalty and not the design.
#' And the transform is kept, so prediction at new points is the parent's
#' evaluation multiplied by it, like any other transformed basis.
#'
#' The last property is what the separation of a linear from a nonlinear effect
#' needs: a reparametrization that does not satisfy
#' \eqn{(\mathbf{1}, x)^\top Z = 0} estimates the sum of the two correctly and
#' the split between them with bias.
#'
#' @param basis The basis to transform.
#' @param x The points the empirical inner product is taken at, normally the
#'   observed values of the covariate.
#' @param penalty A square penalty matrix with one row per basis function.
#'   Defaults to the second-derivative Gram matrix, the integrated squared
#'   second derivative. A discrete difference penalty is passed explicitly, for
#'   instance \code{crossprod(diff(diag(k), differences = 2))}.
#' @param constraints A matrix with one row per constraint and one column per
#'   evaluation point, whose row space the result is made empirically
#'   orthogonal to. Defaults to a constant and \code{x}, which is the
#'   separation of a linear from a nonlinear effect.
#' @param scale Whether to rescale so that
#'   \eqn{\mathrm{tr}(Z^\top Z/n) = 1}, which puts the bases of different terms
#'   on a common scale.
#'
#' @return An object of class \code{\link{TransformedBasis}}.
#'
#' @references
#' Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
#' smoothing. \emph{Numerische Mathematik} 24, 375-382.
#'
#' @seealso \code{\link{constrain_basis}}, \code{\link{basis_gram}}
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(200))
#' d <- dr_basis(bspline_basis(dimension = 12), x)
#'
#' z <- basis_eval(d, x)
#' round(max(abs(crossprod(z)[upper.tri(crossprod(z))])), 10) # diagonal
#' round(max(abs(crossprod(cbind(1, x), z))), 10) # orthogonal to 1 and x
#'
#' @export
dr_basis <- function(basis, x, penalty = NULL, constraints = NULL,
                     scale = TRUE) {
  x <- check_eval_points(basis, x)
  if (anyNA(x)) stop("'x' must not contain missing values.", call. = FALSE)
  k <- basis@dimension
  n <- length(x)

  if (is.null(penalty)) penalty <- basis_gram(basis, order = 2L)
  if (!is.matrix(penalty) || any(dim(penalty) != c(k, k))) {
    stop(sprintf("'penalty' must be a %d by %d matrix.", k, k), call. = FALSE)
  }

  b <- basis_eval(basis, x)
  cons <- if (is.null(constraints)) t(cbind(1, x)) else as.matrix(constraints)
  if (ncol(cons) != n) {
    stop("'constraints' must have one column per element of 'x'.", call. = FALSE)
  }

  # Step 1: restrict to the null space of the constraint, expressed in the
  # coefficients. This is where the empirical orthogonality comes from.
  cm <- cons %*% b
  s <- svd(cm, nu = 0L, nv = k)
  rank <- sum(s$d > 1e-10 * max(s$d, 0))
  if (rank >= k) {
    stop("The constraint leaves no functions in the basis.", call. = FALSE)
  }
  v0 <- s$v[, seq.int(rank + 1L, k), drop = FALSE]

  # Step 2: diagonalize the pencil. The PENALTY is the matrix factorized, not
  # the design, which is what lets a rank-deficient B through.
  g_tilde <- crossprod(b %*% v0) / n
  p_tilde <- crossprod(v0, penalty %*% v0)
  p_tilde <- (p_tilde + t(p_tilde)) / 2

  r <- chol_pd(p_tilde)
  if (is.null(r)) {
    stop(
      "The penalty is singular on the constrained space: some direction is ",
      "neither penalized nor identified there. Reduce 'dimension', or supply ",
      "a penalty with a smaller null space.",
      call. = FALSE
    )
  }
  r_inv <- backsolve(r, diag(nrow(r)))
  m <- crossprod(r_inv, g_tilde %*% r_inv)
  m <- (m + t(m)) / 2
  e <- eigen(m, symmetric = TRUE)
  a <- r_inv %*% e$vectors

  tm <- v0 %*% a
  if (scale) tm <- tm / sqrt(sum(e$values))

  new_transformed(
    basis, tm, paste0("dr(", basis@basis_name, ")"), "dr",
    params = list(
      empirical_variance = if (scale) e$values / sum(e$values) else e$values,
      constraint_rank = rank
    )
  )
}
