#' @include numerical_fallbacks.R
NULL


#' Linearly Transformed Basis
#'
#' @description
#' The S7 class of bases obtained from another by a fixed linear map of its
#' functions, \eqn{\tilde{B}(x) = B(x)\,T}. It carries the parent and the
#' matrix, and it is itself a basis, so every generic applies to a
#' constrained or rotated basis and code that receives it works in the same
#' way as for the parent. Constructed by [orthonorm_basis()],
#' [constrain_basis()] or [dr_basis()].
#'
#' @details
#' # Orthonormalization, constraints and rotations
#'
#' Orthonormalizing a basis, restricting it to satisfy a linear constraint,
#' and rebuilding it so that it diagonalizes an inner product are the same
#' operation with different matrices. Each generic follows from linearity:
#'
#' \deqn{\tilde{B}^{(d)}(x) = B^{(d)}(x)\,T, \qquad
#'   \int_a^x \tilde{B} = \left(\int_a^x B\right) T, \qquad
#'   \tilde{G} = T^\top G\, T.}
#'
#' Differentiation and integration are linear and \eqn{T} does not depend on
#' \eqn{x}, so both transform by the same matrix; the Gram matrix transforms
#' by congruence. A parent with exact derivatives and an exact Gram matrix
#' therefore passes its exactness on, and [basis_is_numerical()] reports the
#' parent's flags in place of those of this class's own methods.
#'
#' The anchoring of the integral is kept as well: a linear combination of
#' columns that are all zero at the lower endpoint is zero there.
#'
#' # Fewer columns than rows
#'
#' \eqn{T} may be \eqn{K \times m} with \eqn{m < K}, which is how
#' [constrain_basis()] reduces the dimension. `@dimension` is `ncol(transform)`
#' and the parent's is `nrow(transform)`.
#'
#' # Transforms compose by multiplication
#'
#' Transforming a `TransformedBasis` again produces one object holding the
#' product of the two matrices, with the original as its parent. A chain of
#' transformations therefore costs one matrix multiplication per evaluation
#' however long it is. The `@basis_name` still nests, so an orthonormalized
#' orthonormal basis reads `orthonorm(orthonorm(bspline))` while
#' `@parent_basis` is the B-spline.
#'
#' @inheritParams basis
#' @param parent_basis The basis being transformed, any object inheriting from
#'   [basis]. Kept whole, so that it can still be evaluated and inspected.
#' @param transform The matrix \eqn{T}, numeric, with one row per parent
#'   function and one column per function of the result.
#'
#' @return An object of class `TransformedBasis`, inheriting from [basis],
#'   with the five properties of a basis plus `parent_basis` and `transform`.
#'
#' @seealso [orthonorm_basis()], [constrain_basis()] and [dr_basis()], the
#'   three constructors; [new_transformed()], which they share.
#'
#' @examples
#' o <- orthonorm_basis(bspline_basis(dimension = 6))
#' S7::S7_inherits(o, TransformedBasis)
#' dim(o@transform)
#'
#' # A constraint takes a column away.
#' x <- seq(0, 1, length.out = 200)
#' b <- bspline_basis(dimension = 6)
#' dim(constrain_basis(b, colSums(basis_eval(b, x)))@transform)
#'
#' # Transforming twice keeps one matrix and the original parent.
#' oo <- orthonorm_basis(o)
#' dim(oo@transform)
#' class(oo@parent_basis)
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
#' Wraps a basis in a [TransformedBasis], collapsing the transform into the
#' parent's when the parent is already one, so that a chain of transformations
#' is stored as a single matrix. [orthonorm_basis()], [constrain_basis()]
#' and [dr_basis()] all build their result through it.
#'
#' @details
#' When `basis` is itself a `TransformedBasis`, the result holds
#' `basis@transform %*% transform` and the original's parent. The dimension of
#' the result is `ncol(transform)`, and its interval is the parent's, a linear
#' map of the functions leaving the domain alone.
#'
#' The function itself checks nothing. The validator of [TransformedBasis]
#' rejects a `transform` whose row count differs from the parent's dimension
#' when the object is built.
#'
#' @param basis The basis to transform, any object inheriting from [basis].
#' @param transform The matrix \eqn{T}, with one row per function of `basis`.
#' @param name The `basis_name` of the result. The callers pass
#'   `orthonorm(...)`, `constrained(...)` or `dr(...)` wrapped around the
#'   parent's own name, so the name nests even where the matrices collapse.
#' @param prefix The two-letter prefix its column names are numbered under:
#'   `"on"`, `"cn"` or `"dr"`.
#' @param params Extra entries for `basis_params`, merged with `prefix`.
#'   [constrain_basis()] adds `constraint_rank` and [dr_basis()] adds that and
#'   `empirical_variance`.
#'
#' @return An object of class [TransformedBasis].
#'
#' @seealso [TransformedBasis], the class it builds.
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
#'
#' @description
#' Numbers the columns under the two-letter prefix the transformation
#' recorded: `on1 ... onk` after [orthonorm_basis()], `cn1 ...` after
#' [constrain_basis()], `dr1 ...` after [dr_basis()].
#'
#' @details
#' The parent's names cannot be carried over. Each new function is a
#' combination of all the parent's, so none of the parent's names fits it,
#' and a constraint returns fewer functions than it consumed. The prefix
#' records which transformation produced the column. The default
#' [basis_colnames.basis()] takes the first two characters of `@basis_name`,
#' which here begins with the name of the transformation, and would give
#' `or1` for both an orthonormalization and its re-orthonormalization.
#'
#' @param basis A [TransformedBasis] object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A character vector of length `basis@dimension`.
#'
#' @seealso [basis_colnames()] for the generic.
#'
#' @keywords internal
S7::method(basis_colnames, TransformedBasis) <- function(basis, ...) {
  paste0(basis@basis_params$prefix, seq_len(basis@dimension))
}


#' Evaluate a Transformed Basis
#'
#' @name basis_eval.TransformedBasis
#'
#' @description
#' Evaluates the parent and multiplies by the transform,
#' \eqn{\tilde{B}(x) = B(x)\,T}, with no further computation. The result
#' spans a subspace of what the parent spans, and all of it when \eqn{T} is
#' square and invertible. Properties of the individual columns, such as the
#' partition of unity of a B-spline, are in general not kept.
#'
#' @details
#' Cost is one parent evaluation plus one `length(x)` by `K` by `m` matrix
#' product. Because [new_transformed()] collapses a chain into one matrix,
#' that is the cost however many transformations were composed.
#'
#' @param basis A [TransformedBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval,
#'   which is the parent's.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names `on1`, `cn1` or `dr1` and so on.
#'
#' @seealso [basis_eval()] for the generic; [TransformedBasis] for the class.
#'
#' @keywords internal
S7::method(basis_eval, TransformedBasis) <- function(basis, x, ...) {
  name_columns(basis_eval(basis@parent_basis, x) %*% basis@transform, basis)
}


#' Derivatives of a Transformed Basis
#'
#' @name basis_deriv.TransformedBasis
#'
#' @description
#' Returns the parent's `order`-th derivative multiplied by the transform,
#' \eqn{\tilde{B}^{(d)}(x) = B^{(d)}(x)\,T}. Differentiation is linear and
#' \eqn{T} does not depend on `x`, so the same matrix serves at every order
#' and no derivative of the transformation itself enters.
#'
#' @details
#' The accuracy is the parent's. Where the parent differentiates exactly so
#' does this; where the parent falls back to a stencil the result is that
#' stencil's result multiplied by \eqn{T}, which is why
#' [basis_is_numerical()] reports the parent's flags for a transformed basis
#' in place of those of its own registered methods.
#'
#' @param basis A [TransformedBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param order The derivative order, default `1`, passed to the parent
#'   unchanged: a single non-negative whole number, or one per variable for
#'   a parent of several variables.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns.
#'
#' @seealso [basis_deriv()] for the generic; [basis_is_numerical()] for whose
#'   accuracy this inherits.
#'
#' @keywords internal
S7::method(basis_deriv, TransformedBasis) <- function(basis, x, order = 1L, ...) {
  d <- basis_deriv(basis@parent_basis, x, order = order)
  name_columns(d %*% basis@transform, basis)
}


#' Integral of a Transformed Basis
#'
#' @name basis_int.TransformedBasis
#'
#' @description
#' Returns the parent's anchored integral multiplied by the transform. The
#' anchoring needs no correction: every column of the parent's integral is
#' zero at the lower endpoint, and a linear combination of zeros is zero.
#'
#' @details
#' Integration is linear and \eqn{T} does not depend on `x`, so the identity
#' is exactly the one [basis_deriv.TransformedBasis()] uses in the other
#' direction. The accuracy is again the parent's.
#'
#' @param basis A [TransformedBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, exactly zero in the row at `basis@lower`.
#'
#' @seealso [basis_int()] for the generic and the anchoring convention.
#'
#' @keywords internal
S7::method(basis_int, TransformedBasis) <- function(basis, x, ...) {
  name_columns(basis_int(basis@parent_basis, x) %*% basis@transform, basis)
}


#' Gram Matrix of a Transformed Basis
#'
#' @name basis_gram.TransformedBasis
#'
#' @description
#' Returns the congruence \eqn{T^\top G\,T} of the parent's Gram matrix. A
#' parent whose inner products are exact passes that exactness on, so the
#' Gram matrix of an orthonormalized B-spline is the identity up to rounding
#' and no quadrature is run.
#'
#' @details
#' The congruence preserves symmetry and positive semidefiniteness, and its
#' rank is at most the smaller of the rank of \eqn{G} and `ncol(T)`. The
#' result is symmetrized as `(G + t(G))/2` before it is returned.
#'
#' [orthonorm_basis()] relies on this: it chooses \eqn{T} so that
#' \eqn{T^\top G\,T} is the identity.
#'
#' @param basis A [TransformedBasis] object.
#' @param order The derivative order, a single non-negative whole number,
#'   default `0`, passed to the parent unchanged.
#' @param at,weight Handled in the body of [basis_gram()] before dispatch, so
#'   they never arrive here. Named only because S7 requires a method's formals
#'   to contain the generic's.
#' @param ... Passed to the parent's method, so `panels` and `nodes` reach a
#'   parent whose Gram matrix is a quadrature.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with [basis_colnames()] on both margins.
#'
#' @seealso [basis_gram()] for the generic; [orthonorm_basis()], which chooses
#'   \eqn{T} to make this the identity.
#'
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


#' Cholesky Factor with a Positive-Definiteness Test
#'
#' @description
#' Returns the upper triangular Cholesky factor of a symmetric matrix, or
#' `NULL` when the matrix is not positive definite to the given relative
#' tolerance. The test is made on the eigenvalues, so that its outcome is the
#' same on every platform.
#'
#' @details
#' On a matrix with an exactly zero eigenvalue the pivot that should be zero
#' comes out positive or negative according to rounding, so [base::chol()]
#' succeeds on some platforms and fails on others. A construction that relied
#' on `chol()` to decide whether a Gram matrix or a penalty is usable would
#' therefore give different results on different platforms.
#'
#' Comparing the smallest eigenvalue with the largest depends only on the
#' matrix and gives the same result everywhere. The eigendecomposition costs
#' little beside the other work of the callers.
#'
#' The `chol()` call is still wrapped, so a matrix that passes the eigenvalue
#' test and fails the factorization returns `NULL` instead of signalling an
#' error.
#'
#' @param m A symmetric numeric matrix. Symmetry is assumed, never checked:
#'   `eigen(symmetric = TRUE)` reads the lower triangle.
#' @param tol The relative tolerance below which the smallest eigenvalue
#'   counts as zero, default `1e-12`. The test is
#'   `min(ev) > tol * max(ev)`, so it is scale-free.
#'
#' @return The upper triangular Cholesky factor \eqn{R} with
#'   \eqn{m = R^\top R}, or `NULL` when `m` is empty, holds a value that is not
#'   finite, or is not positive definite to `tol`.
#'
#' @seealso [orthonorm_basis()], [dr_basis()] and [smoother_reparam()], its
#'   callers, which turn a `NULL` into an error stating what to change.
#'
#' @keywords internal
chol_pd <- function(m, tol = 1e-12) {
  # eigen() itself throws on an empty matrix and on a missing or infinite entry
  if (!length(m) || !all(is.finite(m))) return(NULL)
  ev <- eigen(m, symmetric = TRUE, only.values = TRUE)$values
  if (max(ev) <= 0 || min(ev) <= tol * max(ev)) return(NULL)
  tryCatch(chol(m), error = function(e) NULL)
}


#' Orthonormalize a Basis
#'
#' @description
#' Returns a basis spanning the same functions whose Gram matrix is the
#' identity, so that the functions are orthonormal in \eqn{L^2}. The
#' transform is read off the Gram matrix, so for a basis whose Gram matrix is
#' exact (the B-spline and Legendre families, and a Fourier basis over a whole
#' period) the orthonormalization is exact, with no grid and no number of
#' points to choose.
#'
#' @details
#' # The construction
#'
#' Write \eqn{G = R^\top R} for the Cholesky factorization of the Gram matrix.
#' The basis \eqn{B R^{-1}} then has Gram matrix
#' \eqn{R^{-\top} R^\top R\, R^{-1} = I}, so \eqn{T = R^{-1}}. The Gram
#' matrix of the result is the identity up to rounding.
#'
#' The span is unchanged, \eqn{R^{-1}} being invertible: a function the parent
#' can represent is fitted by the orthonormalized basis to rounding.
#'
#' # Orthonormal in which inner product
#'
#' `order` chooses it. At `0`, the default, the functions themselves are
#' orthonormal. Above that the `order`-th derivatives would be, and the Gram
#' matrix there is singular for every family, the constant differentiating
#' away, so the factorization fails and an error is signalled.
#' Orthonormalizing a derivative therefore needs a basis whose constant has
#' already been removed by [constrain_basis()].
#'
#' # Composing
#'
#' Orthonormalizing an already orthonormal basis returns it unchanged up to
#' rounding, and the two transforms collapse into one matrix; only
#' `@basis_name` records the second pass.
#'
#' @param basis The basis to orthonormalize, any object inheriting from
#'   [basis].
#' @param order The derivative order whose inner products are made the
#'   identity, a single non-negative whole number, or one per variable for a
#'   basis of several variables; default `0`. Any value above `0` signals an
#'   error for the shipped families, their higher-order Gram matrices being
#'   singular.
#'
#' @return An object of class [TransformedBasis] with `basis_name`
#'   `orthonorm(<parent>)` and column names `on1`, `on2`, and so on. Its
#'   dimension is the parent's.
#'
#' @seealso [basis_gram()], which supplies \eqn{G} and which the result makes
#'   the identity; [constrain_basis()] to remove a direction first;
#'   [dr_basis()] for a rotation that diagonalizes two matrices at once.
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' o <- orthonorm_basis(b)
#'
#' # The Gram matrix is the identity up to rounding.
#' round(basis_gram(o), 12)
#' max(abs(basis_gram(o) - diag(6)))
#'
#' # The span is unchanged: a function in the span of the parent is fitted
#' # exactly.
#' set.seed(1)
#' x <- seq(0, 1, length.out = 100)
#' f <- drop(basis_eval(b, x) %*% rnorm(6))
#' max(abs(lm.fit(basis_eval(o, x), f)$residuals))
#'
#' # Orthonormalizing again changes nothing, and keeps one matrix.
#' max(abs(basis_gram(orthonorm_basis(o)) - diag(6)))
#' class(orthonorm_basis(o)@parent_basis)
#'
#' # A higher order signals an error, that Gram matrix being singular.
#' try(orthonorm_basis(b, order = 2))
#'
#' @export
orthonorm_basis <- function(basis, order = 0L) {
  order <- check_order(order, basis_nvar(basis))
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


#' The Null Basis of a Constraint
#'
#' @description
#' An orthonormal basis of the null space of `cm`, from whose columns a
#' constrained basis is built. It is the transform **relative to the basis
#' being constrained**, which is not always the one stored by the resulting
#' object.
#'
#' @details
#' [new_transformed()] flattens a nested transform, so a basis that is
#' itself transformed comes back carrying the product against its own
#' parent: constraining a cyclic smoother's twelve periodic functions,
#' which are a transform of fifteen B-splines, gives an object whose
#' transform is 15 by 11 and not 12 by 11. The stored transform serves the
#' evaluation, and a matrix defined on the twelve functions needs the local
#' one, which is why it is returned here.
#'
#' @param cm The constraint, one row per direction removed and one column
#'   per basis function.
#' @param dimension The number of basis functions.
#' @param tol The relative tolerance at which a singular value counts as
#'   zero.
#'
#' @return A numeric matrix of `dimension` rows and one column per
#'   direction that survives the constraint.
#'
#' @seealso [constrain_basis()] and [smoother_reparam()], the two callers.
#'
#' @keywords internal
constraint_null <- function(cm, dimension, tol = 1e-10) {
  s <- svd(cm, nu = 0L, nv = dimension)
  rank <- sum(s$d > tol * max(s$d, 0))
  if (rank >= dimension) {
    stop(
      "The constraint leaves no functions: its rank equals the dimension of ",
      "the basis.",
      call. = FALSE
    )
  }
  s$v[, seq.int(rank + 1L, dimension), drop = FALSE]
}


#' Restrict a Basis to a Linear Constraint
#'
#' @description
#' Returns a basis whose functions are exactly those of the parent satisfying
#' \eqn{C\beta = 0}, with the dimension reduced by the rank of \eqn{C}. This
#' is how a basis that carries its own constant is made to sit beside an
#' intercept, and how a smooth term is made identifiable.
#'
#' @details
#' # The construction
#'
#' The transform is an orthonormal basis of the null space of \eqn{C}, taken
#' from the right singular vectors of its singular value decomposition
#' belonging to the zero singular values. The constrained basis therefore
#' spans precisely the admissible functions, and rescaling or reordering the
#' rows of \eqn{C} changes nothing about the space produced.
#'
#' The rank is counted as the number of singular values above
#' `tol * max(d)`, and the result has `basis@dimension - rank` columns, with
#' the rank recorded in `basis_params$constraint_rank`.
#'
#' # Writing the constraint
#'
#' `constraint` has one column per basis function, so a condition on the
#' fitted curve is expressed by first evaluating the basis. The usual
#' sum-to-zero identifiability constraint over a grid or over the observed
#' covariate is `colSums(basis_eval(b, x))`, which makes the fitted values sum
#' to zero there; several constraints are the rows of a matrix.
#'
#' The function applies a given constraint. The choice of constraint for a
#' model term, such as a sum-to-zero condition for identifiability or
#' orthogonality to a linear part, depends on the meaning of the term and is
#' made in the modeling layer.
#'
#' # Errors
#'
#' A `constraint` that is not numeric, or whose column count is not
#' `basis@dimension`, signals an error that states the required count. A
#' missing value signals an error, and so does a constraint of full rank,
#' which would leave a basis of zero columns.
#'
#' @param basis The basis to restrict, any object inheriting from [basis].
#' @param constraint A numeric matrix with one column per basis function, or a
#'   plain vector for a single constraint, which is taken as one row. Rows
#'   need not be independent; the rank is computed.
#' @param tol The relative tolerance below which a singular value counts as
#'   zero when the rank is determined, default `1e-10`, applied as
#'   `tol * max(d)`.
#'
#' @return An object of class [TransformedBasis] of dimension
#'   `basis@dimension - rank`, with `basis_name` `constrained(<parent>)`,
#'   column names `cn1`, `cn2`, and so on, and `basis_params$constraint_rank`.
#'
#' @seealso [orthonorm_basis()] for a transformation that keeps the dimension;
#'   [dr_basis()], which applies a constraint and a rotation together.
#'
#' @examples
#' b <- bspline_basis(dimension = 6)
#' x <- seq(0, 1, length.out = 200)
#'
#' # Sum to zero over a grid: the usual identifiability constraint, which
#' # removes the constant a B-spline basis carries.
#' cs <- constrain_basis(b, colSums(basis_eval(b, x)))
#' cs@dimension
#' cs@basis_params$constraint_rank
#' max(abs(colSums(basis_eval(cs, x))))
#'
#' # Two constraints take two columns: sum to zero and orthogonal to x.
#' C <- rbind(colSums(basis_eval(b, x)), colSums(basis_eval(b, x) * x))
#' constrain_basis(b, C)@dimension
#'
#' # A constraint of full rank leaves no functions and signals an error.
#' try(constrain_basis(b, diag(6)))
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

  tm <- constraint_null(cm, basis@dimension, tol)
  rank <- basis@dimension - ncol(tm)
  new_transformed(
    basis, tm, paste0("constrained(", basis@basis_name, ")"), "cn",
    params = list(constraint_rank = rank)
  )
}


#' Demmler-Reinsch Basis
#'
#' @description
#' Returns the basis that simultaneously diagonalizes the empirical inner
#' product at the given points and a penalty, and is empirically orthogonal to
#' a constant and to `x`. Its columns are ordered from the smoothest to the
#' wiggliest, each carrying a known share of the empirical variance, so a
#' penalized fit against it is a shrinkage of independent coordinates and the
#' linear part of the effect is separated from the nonlinear part exactly.
#'
#' @details
#' # The construction
#'
#' Three steps. The constraint matrix \eqn{C = (\mathbf{1}, x)^\top B} is
#' formed and the basis restricted to its null space \eqn{V_0}, which makes
#' the remaining functions empirically orthogonal to a constant and to
#' \eqn{x}. The pencil
#' \eqn{(V_0^\top (B^\top B/n) V_0,\; V_0^\top P V_0)} is then diagonalized,
#' and the transform is \eqn{T = V_0 A}.
#'
#' The resulting design matrix \eqn{Z = B T} then has \eqn{Z^\top Z} diagonal
#' and satisfies \eqn{(\mathbf{1}, x)^\top Z = 0}, both up to rounding.
#'
#' # What `scale` does to the penalty
#'
#' \eqn{T^\top P T} is **proportional to** the identity, and equal to it only
#' when `scale = FALSE`. With `scale = TRUE`, the default, \eqn{T} is divided
#' by \eqn{\sqrt{\sum_j \lambda_j}}, so \eqn{T^\top P T} is
#' \eqn{(\sum_j \lambda_j)^{-1} I} and `tr(Z'Z/n)` is 1. At `scale = FALSE`
#' \eqn{T^\top P T} is the identity and `tr(Z'Z/n)` is
#' \eqn{\sum_j \lambda_j}.
#'
#' The scaling puts the bases of different terms on a common footing, so one
#' smoothing parameter means the same thing across them. A consumer that
#' needs the penalty to be exactly \eqn{I} passes `scale = FALSE`.
#'
#' `basis_params$empirical_variance` holds the eigenvalues, normalized to sum
#' to one when `scale = TRUE`. They equal `diag(Z'Z/n)` up to rounding, and
#' are the share of empirical variance each column carries, falling from the
#' smoothest column to the wiggliest.
#'
#' # Cost and properties
#'
#' It factorizes only a \eqn{q \times K} matrix and a
#' \eqn{(K-q) \times (K-q)} one, never anything of the size of the sample. It
#' tolerates a rank-deficient \eqn{B}, which equally spaced knots produce when
#' the data leave `degree + 1` or more consecutive knot spans of a B-spline
#' empty, the matrix inverted being the penalty. The transform is kept, so prediction at new points is the
#' parent's evaluation multiplied by it, as for any other transformed basis.
#'
#' The orthogonality \eqn{(\mathbf{1}, x)^\top Z = 0} is what separates a
#' linear from a nonlinear effect: a reparametrization without it estimates
#' the sum of the two correctly and the split between them with bias.
#'
#' # Errors
#'
#' A missing value in `x`, a `penalty` that is not `K` by `K`, or a
#' `constraints` with the wrong number of columns each signal an error. A
#' penalty singular on the constrained space signals an error that names two
#' remedies: there is then a direction that is neither penalized nor
#' identified, and a rotation of the basis cannot remove it.
#'
#' @param basis The basis to transform, any object inheriting from [basis].
#' @param x The points the empirical inner product is taken at, normally the
#'   observed covariate. A numeric vector inside the basis interval, with no
#'   missing values. The construction depends on where they lie, so a basis
#'   built for one sample is not the basis for another.
#' @param penalty A square numeric penalty matrix with one row and column per
#'   basis function. `NULL`, the default, uses `basis_gram(basis, order = 2)`,
#'   the integrated squared second derivative. A discrete difference penalty
#'   is passed explicitly, for instance
#'   `crossprod(diff(diag(k), differences = 2))`.
#' @param constraints A matrix with one row per constraint and one column per
#'   element of `x`, whose row space the result is made empirically orthogonal
#'   to. `NULL`, the default, uses `t(cbind(1, x))`, which is the separation
#'   of a linear from a nonlinear effect.
#' @param scale Whether to rescale so that
#'   \eqn{\mathrm{tr}(Z^\top Z/n) = 1}, default `TRUE`. See above for what it
#'   does to \eqn{T^\top P T}.
#'
#' @return An object of class [TransformedBasis] of dimension
#'   `basis@dimension - rank(C)`, with `basis_name` `dr(<parent>)`, column
#'   names `dr1`, `dr2`, and so on ordered from smoothest to wiggliest, and
#'   `basis_params` holding `empirical_variance` and `constraint_rank`.
#'
#' @references
#' Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
#' smoothing. *Numerische Mathematik* **24**, 375-382.
#'
#' @seealso [constrain_basis()], which does its first step alone;
#'   [basis_gram()], which supplies the default penalty;
#'   [orthonorm_basis()] for a rotation that diagonalizes one matrix.
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(200))
#' d <- dr_basis(bspline_basis(dimension = 12), x)
#'
#' # Two columns went to the constraint, leaving ten.
#' c(d@dimension, d@basis_params$constraint_rank)
#'
#' z <- basis_eval(d, x)
#'
#' # Z'Z is diagonal, and Z is orthogonal to a constant and to x.
#' max(abs(crossprod(z)[upper.tri(crossprod(z))]))
#' max(abs(crossprod(cbind(1, x), z)))
#'
#' # The columns run from smoothest to wiggliest, and the recorded shares of
#' # empirical variance are the diagonal of Z'Z/n.
#' round(d@basis_params$empirical_variance, 6)
#' max(abs(d@basis_params$empirical_variance - diag(crossprod(z)) / length(x)))
#'
#' # Scaled, the penalty is a multiple of the identity and the trace is one.
#' P <- basis_gram(bspline_basis(dimension = 12), order = 2)
#' round(unique(round(diag(crossprod(d@transform, P %*% d@transform)), 6)), 4)
#' sum(diag(crossprod(z))) / length(x)
#'
#' # Unscaled, the penalty is the identity instead.
#' du <- dr_basis(bspline_basis(dimension = 12), x, scale = FALSE)
#' max(abs(crossprod(du@transform, P %*% du@transform) - diag(du@dimension)))
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
