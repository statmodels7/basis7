#' A Smoother: the Four Decisions of a Penalized Smooth
#' @name smoother
#'
#' @description
#' A smoother is a recipe for a penalized smooth. It carries the four
#' decisions a smooth is made of, which are independent of one another and
#' which a single basis does not determine:
#'
#' 1. which functions span the space, the **basis**;
#' 2. what counts as roughness, the **penalty**;
#' 3. which directions the penalty leaves alone and what becomes of them,
#'    the **null space**;
#' 4. which coordinates the coefficients live in, the **reparametrization**.
#'
#' It is a recipe rather than a built object because the third and fourth
#' decisions need the data: the Demmler-Reinsch rotation diagonalizes the
#' pencil of the empirical Gram matrix against the penalty, and the interval
#' of the default basis is read from the covariate. [smoother_build()]
#' resolves a smoother at a vector of covariate values and returns the block
#' and its penalty matrix; [smoother_apply()] reapplies the transform
#' recorded there at new values.
#'
#' # Why one object rather than a basis and a penalty
#'
#' The two are not independently choosable. The null space is a property of
#' the **pair**: the second-derivative Gram matrix of a cubic B-spline has a
#' two-dimensional null space, that of a Fourier basis is one-dimensional at
#' every order because the basis contains no linear function, and the Gram
#' matrix of a B-spline of degree \eqn{d} at an order above \eqn{d} is
#' identically zero, which penalizes nothing. Passing a basis and a penalty
#' separately would leave the caller to satisfy that compatibility at every
#' call site. A smoother is one object, so each constructor validates its own
#' arguments where the caller wrote them.
#'
#' # A penalty of your own
#'
#' `penalty` replaces the roughness matrix with a penalty built by a factory
#' of the coefficient count: `bspline_smooth(penalty =
#' penalties7::lasso_penalty)`. It is a factory and not a built penalty
#' because how many coefficients a smooth has is settled by the data, the
#' constraint and the null space moving with them.
#'
#' A smoother **stores the function and never calls it**. This package sits
#' at the bottom of the dependency graph and imports \pkg{numericals7}
#' alone, so it cannot name \pkg{penalties7} and cannot ask whether what the
#' function returns is a penalty; [check_penalty()] asks only that it be a
#' function of one argument. Whichever layer builds the term calls it, at
#' the count only the data settle, and checks the result there.
#'
#' The construction is unaffected. [smoother_build()] returns the roughness
#' matrix in `S` whether or not a factory is given, because the
#' reparametrization reads that matrix: it is what orders the coordinates
#' from the smoothest to the most wiggly and makes the penalty on them the
#' identity, and that ordering is the reason a penalty of another shape is
#' worth reaching for. `unpenalized` counts the columns the roughness leaves
#' free, which a caller building a separable penalty needs, since such a
#' penalty has no zero row with which to leave a column alone.
#'
#' `smoother` is abstract: construct one through a family, of which
#' [bspline_smooth()] is the first.
#'
#' @param smoother_name A single string naming the family, printed by
#'   [print.smoother()]. Not read by any computation.
#' @param dimension The number of basis functions before any constraint or
#'   reparametrization, a single positive integer. Must be of storage mode
#'   integer.
#' @param order The order of derivative the penalty integrates, a single
#'   positive integer. It says what a strongly penalized fit contracts
#'   toward: a constant at 1, a straight line at 2, a parabola at 3.
#' @param measure The measure the roughness is integrated against.
#'   `"lebesgue"` is the length measure on the interval.
#' @param constrain The directions the smooth is made orthogonal to, or
#'   `NULL` for the null space of the basis and the penalty together. The
#'   form it takes belongs to the family, a periodic basis having no reading
#'   for "polynomials up to degree c".
#' @param null_space What becomes of the directions the penalty does not
#'   see: `"keep"` leaves them as free columns, `"drop"` removes them,
#'   `"shrink"` penalizes them under the same smoothing parameter.
#' @param reparam The coordinates the coefficients live in, a single string.
#'   `"dr"` is the Demmler-Reinsch rotation.
#' @param penalty `NULL` for the quadratic roughness matrix, or a factory
#'   building a \pkg{penalties7} penalty from a coefficient count.
#' @param lower,upper The endpoints of the interval, each a single finite
#'   number or `NULL` to read it from the data at build. Given, both must be
#'   given, with `lower < upper`.
#' @param smoother_params A named list of whatever else the family needs.
#'
#' @seealso [bspline_smooth()] for the B-spline family, [smoother_build()]
#'   for what a smoother produces at data, [basis] for the object a smoother
#'   builds on.
#'
#' @return An S7 class object. `smoother` itself is abstract and cannot be
#'   constructed; a family constructor such as [bspline_smooth()] returns an
#'   object inheriting from it.
#'
#' @examples
#' sm <- bspline_smooth(k = 10)
#' S7::S7_inherits(sm, smoother)
#' sm@order
#'
#' # Abstract: there is no direct constructor.
#' try(smoother())
#' @export
smoother <- S7::new_class(
  "smoother",
  abstract = TRUE,
  properties = list(
    smoother_name = S7::class_character,
    dimension = S7::class_integer,
    order = S7::class_integer,
    measure = S7::class_any,
    constrain = S7::class_any,
    null_space = S7::class_character,
    reparam = S7::class_character,
    penalty = S7::class_any,
    lower = S7::class_any,
    upper = S7::class_any,
    smoother_params = S7::class_list
  ),
  validator = function(self) {
    if (length(self@smoother_name) != 1L || is.na(self@smoother_name)) {
      return("@smoother_name must be a single string")
    }
    if (length(self@dimension) != 1L || is.na(self@dimension) ||
      self@dimension < 1L) {
      return("@dimension must be a single positive integer")
    }
    if (length(self@order) != 1L || is.na(self@order) || self@order < 1L) {
      return("@order must be a single positive integer")
    }
    if (length(self@null_space) != 1L ||
      !self@null_space %in% c("keep", "drop", "shrink")) {
      return("@null_space must be one of \"keep\", \"drop\", \"shrink\"")
    }
    if (length(self@reparam) != 1L || is.na(self@reparam)) {
      return("@reparam must be a single string")
    }
    bad <- function(v) {
      !is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v))
    }
    if (bad(self@lower) || bad(self@upper)) {
      return("@lower and @upper must each be NULL or one finite number")
    }
    if (!is.null(self@lower) && !is.null(self@upper) &&
      self@lower >= self@upper) {
      return("@lower must be strictly less than @upper")
    }
    NULL
  }
)


#' The B-Spline Smoother Class
#' @name BsplineSmoother
#'
#' @description
#' The class [bspline_smooth()] returns: a B-spline basis with a roughness
#' penalty, the Demmler-Reinsch reparametrization and the linear column kept
#' free. It adds `degree` to the properties of [smoother].
#'
#' @inheritParams smoother
#' @param degree The degree of the B-spline pieces, a single positive
#'   integer, at most `dimension - 1`.
#'
#' @return An S7 object of class `BsplineSmoother`, inheriting from
#'   [smoother]. Construct one with [bspline_smooth()], which validates its
#'   arguments; the class constructor does not.
#'
#' @seealso [bspline_smooth()], which is the way to build one.
#'
#' @examples
#' sm <- bspline_smooth(k = 12, degree = 3)
#' c(class = class(sm)[1], degree = sm@degree, dimension = sm@dimension)
#' @export
BsplineSmoother <- S7::new_class(
  "BsplineSmoother",
  parent = smoother,
  properties = list(degree = S7::class_integer),
  validator = function(self) {
    if (length(self@degree) != 1L || is.na(self@degree) || self@degree < 1L) {
      return("@degree must be a single positive integer")
    }
    if (self@dimension < self@degree + 1L) {
      return("@dimension must be at least @degree + 1")
    }
    NULL
  }
)


#' A B-Spline Smoother
#'
#' @description
#' The B-spline smoother: `k` B-spline functions of degree `degree` over an
#' interval, penalized by the integrated squared derivative of order `order`,
#' rotated to the Demmler-Reinsch coordinates and carrying the unpenalized
#' direction as a free column. It is the construction [modelterms7::s()] has
#' always used, named as an object.
#'
#' @details
#' # The construction
#'
#' At a vector of covariate values [smoother_build()] performs, in order:
#'
#' 1. a [bspline_basis()] of `k` functions over the interval, which is
#'    `lower` and `upper` when both are given and otherwise the range of the
#'    data padded by a thousandth of its width;
#' 2. [dr_basis()], which restricts the basis to the orthogonal complement of
#'    the constant and the linear function over the observed values, then
#'    diagonalizes the pencil of the empirical Gram matrix against the
#'    roughness matrix. The result is a basis whose columns are orthogonal
#'    over the data and ordered from the smoothest to the most oscillatory,
#'    and whose penalty is the identity;
#' 3. with `null_space = "keep"`, the standardized covariate prepended as a
#'    free column, so the penalty is `diag(0, 1, ..., 1)` and a strongly
#'    penalized fit contracts to a straight line rather than to a constant.
#'
#' A basis of `k` functions therefore gives `k - 1` columns when the null
#' space is kept and `k - 2` when it is dropped: the constraint against the
#' constant and the linear function removes two directions, and the free
#' column adds one back.
#'
#' # What `order` means
#'
#' `order` is the order of derivative the penalty integrates, and it fixes
#' what a fit contracts toward as the smoothing parameter grows: a constant
#' at `order = 1`, a straight line at `order = 2`, a parabola at `order = 3`.
#' The null space of the penalty is the polynomials of degree below `order`,
#' so `order` may not exceed `degree`; above it the roughness matrix is
#' identically zero and penalizes nothing.
#'
#' @param k The number of basis functions, a whole number of at least 3. Two
#'   directions are removed by the constraint, so a smaller `k` leaves
#'   nothing to smooth.
#' @param degree The degree of the B-spline pieces, a whole number of at
#'   least 1. `3`, the default, is the cubic spline.
#' @param order The order of derivative the penalty integrates. See the
#'   section above.
#' @param measure The measure the roughness is integrated against.
#' @param constrain The directions the smooth is made orthogonal to. `NULL`,
#'   the default, is the null space of the basis and the penalty together.
#' @param null_space What becomes of the directions the penalty does not
#'   see: `"keep"` leaves them as free columns, `"drop"` removes them,
#'   `"shrink"` penalizes them under the same smoothing parameter.
#' @param reparam The coordinates the coefficients live in. `"dr"`, the
#'   default, is the Demmler-Reinsch rotation.
#' @param penalty `NULL` for the quadratic roughness penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param lower,upper The interval. `NULL`, the default for each, reads it
#'   from the data at build; give both to fix it, which is what a prediction
#'   outside the observed range needs.
#'
#' @return An S7 object of class [BsplineSmoother], inheriting from
#'   [smoother]. It is a recipe: pass it to [smoother_build()] with the
#'   covariate to obtain the block and its penalty matrix.
#'
#' @references
#' Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
#' smoothing. *Numerische Mathematik*, 24, 375--382.
#'
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11, 89--121.
#'
#' @seealso [smoother_build()] for what it produces at data, [smoother] for
#'   the four decisions it carries, [bspline_basis()] for the basis alone.
#'
#' @examples
#' sm <- bspline_smooth(k = 10)
#' sm
#'
#' # k = 10 gives nine columns: the free linear column and eight deviations.
#' set.seed(1)
#' x <- sort(runif(200, -2, 3))
#' out <- smoother_build(sm, x)
#' dim(out$X)
#' out$names
#'
#' # The penalty is diag(0, 1, ..., 1), and one column is unpenalized.
#' round(diag(out$S), 6)
#' out$unpenalized
#'
#' # The free column is the standardized covariate, orthogonal to the rest
#' # over the observed values.
#' round(cor(out$X[, 1], x), 12)
#' max(abs(crossprod(out$X[, 1], out$X[, -1])))
#'
#' # Dropping the null space removes it.
#' dim(smoother_build(bspline_smooth(k = 10, null_space = "drop"), x)$X)
#'
#' # A penalty factory is STORED AND NEVER CALLED here: one that raises
#' # still builds, because it is the model layer that calls it.
#' sm2 <- bspline_smooth(k = 10, penalty = function(n_coef) stop("not here"))
#' out2 <- smoother_build(sm2, x)
#' identical(out2$S, out$S)
#'
#' # It must be a function of the count, and it cannot be combined with a
#' # shrunk null space, which is a weight inside the matrix it replaces.
#' try(bspline_smooth(k = 10, penalty = 3))
#' try(bspline_smooth(k = 10, penalty = function(n) n, null_space = "shrink"))
#'
#' # 'k' must leave something after the constraint.
#' try(bspline_smooth(k = 2))
#' @export
bspline_smooth <- function(k = 10, degree = 3, order = 2,
                           measure = "lebesgue", constrain = NULL,
                           null_space = "keep", reparam = "dr",
                           penalty = NULL, lower = NULL, upper = NULL) {
  # the floor is 2 because how many directions the constraint removes depends
  # on 'order' and 'constrain'; the check that k leaves something to smooth
  # is made below, once both are settled
  k <- check_whole(k, "k", 2L)
  degree <- check_whole(degree, "degree", 1L)
  order <- check_whole(order, "order", 1L)

  if (k < degree + 1L) {
    stop(sprintf(paste0(
      "'k' (%d) is too small for 'degree' (%d): a B-spline basis of degree",
      " m\n  needs at least m + 1 functions."
    ), k, degree), call. = FALSE)
  }
  # ABOVE THE DEGREE THE GRAM MATRIX IS IDENTICALLY ZERO, so the penalty
  # would penalize nothing and every direction would be free. Reported here,
  # where the caller wrote the two numbers, rather than at build.
  if (order > degree) {
    stop(sprintf(paste0(
      "'order' (%d) exceeds 'degree' (%d): the derivative of that order of",
      " a\n  spline of degree m is zero, so the penalty would be the zero",
      " matrix."
    ), order, degree), call. = FALSE)
  }

  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("dr", "none", "orthonorm"))
  check_interval(lower, upper)
  measure <- check_measure(measure)
  penalty <- check_penalty(penalty)
  constrain <- check_constrain(constrain, order)
  # the constraint removes one direction per degree it spans, and a basis
  # with nothing left after it is an error several frames down
  ncon <- if (is.null(constrain)) order else constrain + 1L
  if (k <= ncon) {
    stop(sprintf(paste0(
      "'k' (%d) leaves nothing to smooth: the constraint removes %d",
      " directions.\n  Raise 'k' above %d."
    ), k, ncon, ncon), call. = FALSE)
  }

  sm <- BsplineSmoother(
    smoother_name = "bspline",
    dimension = k,
    degree = degree,
    order = order,
    measure = measure,
    constrain = constrain,
    null_space = null_space,
    reparam = reparam,
    penalty = penalty,
    lower = lower,
    upper = upper,
    smoother_params = list()
  )
  check_available(sm)
  sm
}


#' The Basis a Smoother Builds On
#' @name smoother_basis
#'
#' @description
#' Returns the [basis] a smoother expands the covariate in, before any
#' constraint or reparametrization. It is the one step of the construction
#' that differs between families, so a family is added by registering a
#' method here rather than by rewriting [smoother_build()].
#'
#' @param sm A [smoother].
#' @param x The covariate, a numeric vector. A family whose interval is not
#'   fixed on the object reads it from these values.
#' @param ... Passed to methods.
#'
#' @return A [basis] of `sm@dimension` functions.
#'
#' @seealso [smoother_build()], which calls it.
#'
#' @examples
#' set.seed(1)
#' x <- runif(50)
#' b <- smoother_basis(bspline_smooth(k = 8), x)
#' c(dimension = b@dimension, lower = b@lower, upper = b@upper)
#'
#' # The interval is the range of the data padded by a thousandth of its
#' # width, unless the smoother fixes it.
#' bf <- smoother_basis(bspline_smooth(k = 8, lower = 0, upper = 1), x)
#' c(lower = bf@lower, upper = bf@upper)
#' @export
smoother_basis <- S7::new_generic(
  "smoother_basis", "sm",
  function(sm, x, ...) S7::S7_dispatch()
)

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, BsplineSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  bspline_basis(
    lower = int[[1L]], upper = int[[2L]],
    dimension = sm@dimension, degree = sm@degree
  )
}


#' What a Smoother Removes and What It Gives Back
#' @name smoother_span
#'
#' @description
#' Returns the directions a smoother constrains its block against, and the
#' ones `null_space = "keep"` restores as free columns. Both are evaluated at
#' the covariate, one column per direction.
#'
#' @details
#' # Why the pair and not the basis
#'
#' The null space is a property of the basis and the penalty **together**,
#' not of the basis alone. The second-derivative Gram matrix of a cubic
#' B-spline has the polynomials of degree below 2 in its null space and the
#' order-3 matrix the polynomials of degree below 3, while a Fourier basis
#' has the constant alone at every order, the basis containing no linear
#' function. Measured, the null function of a Fourier Gram matrix has a
#' standard deviation of exactly zero at orders 1, 2 and 3. So the answer
#' belongs to the smoother, which carries both, and a family declares it by
#' registering a method here.
#'
#' # The default, which serves every polynomial family
#'
#' The base method removes the polynomials of degree below `order`, or up to
#' `constrain` when that is given and larger, and restores all but the
#' constant. The constant is not restored because a model carrying an
#' intercept already spans it, which is the convention [dr_basis()] has
#' always followed.
#'
#' The restored columns are the raw powers made orthogonal to one another and
#' to the constant over the observed covariate, then standardized. The first
#' is `(x - mean(x)) / sd(x)`, so at `order = 2` there is exactly one and it
#' is the standardized covariate.
#'
#' Directions of the constraint **beyond** the null space are removed and not
#' restored. Returning them as free columns would put a parametric term
#' inside the smooth instead of in the formula, where a reader sees it.
#'
#' @param sm A [smoother].
#' @param x The covariate, a numeric vector.
#' @param params For `smoother_span_apply()`, the `params` element of the
#'   `smoother_span()` result.
#' @param newx For `smoother_span_apply()`, the new covariate values.
#' @param ... Passed to methods.
#'
#' @return `smoother_span()` returns a list of three elements: `constraint`,
#'   a numeric matrix of one column per direction removed, or `NULL` for
#'   none; `free`, a numeric matrix of one column per direction restored,
#'   with zero columns where none is; and `params`, what
#'   `smoother_span_apply()` needs to rebuild `free` at new values.
#'   `smoother_span_apply()` returns that matrix at `newx`.
#'
#' @seealso [smoother_build()], which calls both.
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(60))
#' sp <- smoother_span(bspline_smooth(k = 8), x)
#' c(removed = ncol(sp$constraint), restored = ncol(sp$free))
#'
#' # the restored column is the standardized covariate
#' round(c(mean = mean(sp$free), sd = sd(sp$free), cor = cor(sp$free, x)), 12)
#'
#' # reapplied at new values rather than rebuilt
#' i <- c(2L, 20L, 55L)
#' max(abs(smoother_span_apply(bspline_smooth(k = 8), sp$params, x[i]) -
#'         sp$free[i, , drop = FALSE]))
#' @export
smoother_span <- S7::new_generic(
  "smoother_span", "sm",
  function(sm, x, ...) S7::S7_dispatch()
)

#' @name smoother_span
#' @keywords internal
S7::method(smoother_span, smoother) <- function(sm, x, ...) {
  ncon <- if (is.null(sm@constrain)) sm@order else sm@constrain + 1L
  # the raw powers, which at order 2 are cbind(1, x): the constraint
  # dr_basis() builds for itself when it is given none
  cons <- outer(x, seq.int(0L, ncon - 1L), "^")
  free <- poly_free(x, sm@order - 1L)
  list(constraint = cons, free = free$free, params = free$params)
}

#' @name smoother_span
#' @export
smoother_span_apply <- S7::new_generic(
  "smoother_span_apply", "sm",
  function(sm, params, newx, ...) S7::S7_dispatch()
)

#' @name smoother_span
#' @keywords internal
S7::method(smoother_span_apply, smoother) <- function(sm, params, newx, ...) {
  poly_free_apply(params, newx)
}


#' The Free Columns of a Polynomial Null Space
#'
#' @description
#' Builds the columns `null_space = "keep"` restores for a family whose null
#' space is the polynomials of degree below `order`: the powers
#' \eqn{x, \ldots, x^{m}} made orthogonal to the constant and to one another
#' over the observed covariate, then standardized.
#'
#' @details
#' The first column is written as `(x - mean(x)) / sd(x)` rather than as the
#' general regression it is a case of. The two agree to the last bit for the
#' quantities themselves, and the literal form is what the construction has
#' always computed at `order = 2`, which is every smooth the toolkit has
#' fitted; a change of arithmetic there would move fits that are not being
#' asked to move.
#'
#' @param x The covariate, a numeric vector.
#' @param m How many columns, `order - 1`. Zero gives a matrix of no columns.
#'
#' @return A list of `free`, the matrix of `m` columns, and `params`, what
#'   [poly_free_apply()] needs to rebuild it at new values.
#'
#' @keywords internal
poly_free <- function(x, m) {
  n <- length(x)
  if (m < 1L) {
    return(list(free = matrix(numeric(0), n, 0L), params = list(steps = list())))
  }
  free <- matrix(0, n, m)
  steps <- vector("list", m)
  for (j in seq_len(m)) {
    if (j == 1L) {
      ctr <- mean(x)
      scl <- stats::sd(x)
      if (!is.finite(scl) || scl == 0) scl <- 1
      free[, 1L] <- (x - ctr) / scl
      steps[[1L]] <- list(center = ctr, scale = scl)
    } else {
      q <- cbind(1, free[, seq_len(j - 1L), drop = FALSE])
      a <- qr.solve(q, x^j)
      res <- x^j - as.vector(q %*% a)
      scl <- stats::sd(res)
      if (!is.finite(scl) || scl == 0) scl <- 1
      free[, j] <- res / scl
      steps[[j]] <- list(coef = a, scale = scl)
    }
  }
  list(free = free, params = list(steps = steps))
}


#' The Free Columns of a Polynomial Null Space at New Values
#'
#' @description
#' Rebuilds what [poly_free()] produced, at new covariate values, from the
#' centering, the orthogonalization coefficients and the scales recorded
#' there. Nothing is recomputed from `newx`.
#'
#' @param params The `params` element of a [poly_free()] result.
#' @param newx The new covariate values, a numeric vector.
#'
#' @return A numeric matrix of `length(newx)` rows and one column per step.
#'
#' @keywords internal
poly_free_apply <- function(params, newx) {
  st <- params$steps
  m <- length(st)
  out <- matrix(0, length(newx), m)
  for (j in seq_len(m)) {
    if (j == 1L) {
      out[, 1L] <- (newx - st[[1L]]$center) / st[[1L]]$scale
    } else {
      q <- cbind(1, out[, seq_len(j - 1L), drop = FALSE])
      out[, j] <- (newx^j - as.vector(q %*% st[[j]]$coef)) / st[[j]]$scale
    }
  }
  out
}


#' The Roughness Matrix a Smoother Penalizes With
#'
#' @description
#' The Gram matrix of the derivative of order `sm@order`, integrated against
#' the measure the smoother carries: the length measure on the interval for
#' `"lebesgue"`, the empirical measure of the covariate for `"empirical"`,
#' and the density a function gives otherwise.
#'
#' @details
#' The measure is not decorative. Measured on a cubic B-spline of twelve
#' functions over \eqn{[-2, 2]} at order 2, the correlation between the
#' Lebesgue Gram matrix and the one weighted by a Gaussian of standard
#' deviation 0.25 is 0.11: they are different penalties, and a fit under one
#' is not a fit under the other.
#'
#' @param sm A [smoother].
#' @param b The basis [smoother_basis()] returned.
#' @param x The covariate, for the empirical measure.
#' @param ... Passed to methods.
#'
#' @return A symmetric numeric matrix of `b@dimension` rows and columns.
#'
#' @seealso [smoother_build()], which calls it, and [basis_gram()], the
#'   generic it reads.
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(100, -2, 2))
#' sm <- bspline_smooth(k = 8, lower = -2, upper = 2)
#' g <- smoother_gram(sm, smoother_basis(sm, x), x)
#' dim(g)
#'
#' # the measure is a different penalty, not a detail
#' sme <- bspline_smooth(k = 8, lower = -2, upper = 2,
#'                       measure = "empirical")
#' round(cor(as.vector(g),
#'           as.vector(smoother_gram(sme, smoother_basis(sme, x), x))), 3)
#' @export
smoother_gram <- S7::new_generic(
  "smoother_gram", "sm",
  function(sm, b, x, ...) S7::S7_dispatch()
)

#' @name smoother_gram
#' @keywords internal
S7::method(smoother_gram, smoother) <- function(sm, b, x, ...) {
  ms <- sm@measure
  if (is.function(ms)) return(basis_gram(b, order = sm@order, weight = ms))
  if (identical(ms, "empirical")) {
    return(basis_gram(b, order = sm@order, at = x))
  }
  basis_gram(b, order = sm@order)
}


#' Build a Smoother at Data
#'
#' @description
#' Resolves a [smoother] at a vector of covariate values and returns the
#' design block, the penalty matrix and what a caller needs to reapply the
#' construction at new values. It is the whole of a penalized smooth's
#' arithmetic: what comes back is the pair \eqn{(X, S)} of
#' \eqn{(X'X + \lambda S)^{-1} X'y}.
#'
#' @details
#' The block carries no column names, no `by` expansion and no label: those
#' belong to whichever model layer places the smooth in a formula. The
#' column names the family gives its own coordinates are returned separately
#' in `names`.
#'
#' The transform is computed on `x` and recorded in `blueprint`, so a
#' prediction at new values reapplies it through [smoother_apply()] rather
#' than rebuilding it. Rebuilding would give a basis over a different
#' interval and a rotation of a different Gram matrix, which is a different
#' function of the covariate.
#'
#' @param sm A [smoother], such as one from [bspline_smooth()].
#' @param x The covariate, a numeric vector with no missing values.
#' @param ... Passed to methods.
#'
#' @return A list of five elements:
#'   \describe{
#'     \item{`X`}{The design block, a numeric matrix of `length(x)` rows and
#'       one column per coordinate, with no dimnames.}
#'     \item{`S`}{The penalty matrix, square of `ncol(X)`.}
#'     \item{`unpenalized`}{The number of leading columns the penalty does
#'       not cover, as an integer.}
#'     \item{`blueprint`}{What [smoother_apply()] needs to reapply the
#'       construction at new values.}
#'     \item{`names`}{One name per column of `X`.}
#'   }
#'
#' @seealso [smoother_apply()] for the block at new values,
#'   [bspline_smooth()] for the smoother this builds.
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(150, 0, 1))
#' out <- smoother_build(bspline_smooth(k = 10), x)
#' vapply(out, function(e) class(e)[1], character(1))
#' dim(out$X)
#'
#' # (X, S) is everything a penalized least-squares fit needs.
#' y <- sin(2 * pi * x) + rnorm(150, sd = 0.2)
#' b <- solve(crossprod(out$X) + 0.1 * out$S, crossprod(out$X, y))
#' round(sqrt(mean((out$X %*% b - sin(2 * pi * x))^2)), 3)
#'
#' # The penalty leaves the first column free.
#' out$unpenalized
#' round(diag(out$S), 6)
#' @export
smoother_build <- S7::new_generic(
  "smoother_build", "sm",
  function(sm, x, ...) S7::S7_dispatch()
)

#' @name smoother_build
#' @keywords internal
S7::method(smoother_build, smoother) <- function(sm, x, ...) {
  check_available(sm)
  x <- check_smoother_x(x)
  b <- smoother_basis(sm, x)
  sp <- smoother_span(sm, x)

  # The Demmler-Reinsch rotation diagonalizes the pencil of the empirical
  # Gram matrix against the roughness matrix, so the penalty on what
  # survives the constraint is the identity. Both arguments are what
  # dr_basis() builds for itself at order 2 with no constraint, which is
  # what keeps the default construction the one the toolkit has always run.
  g <- smoother_gram(sm, b, x)
  cons <- sp$constraint
  rp <- smoother_reparam(sm, b, x, g, cons)
  # evaluated from the basis object smoother_apply() will evaluate, so the
  # two are the same arithmetic rather than the same formula
  z <- basis_eval(rp$basis, x)
  s_pen <- rp$S
  nm <- paste0("z", seq_len(ncol(z)))
  free <- sp$free
  nfree <- ncol(free)

  if (sm@null_space %in% c("keep", "shrink") && nfree > 0L) {
    z <- cbind(free, z)
    nm <- c(free_names(nfree), nm)
    # SHRINK gives the free directions a weight of their own rather than
    # zero, so a strongly penalized fit contracts to nothing and the term
    # can leave the model. See shrink_weight() for the number and where it
    # comes from.
    w <- if (identical(sm@null_space, "shrink")) shrink_weight() else 0
    s_pen <- rbind(
      cbind(diag(w, nfree), matrix(0, nfree, ncol(s_pen))),
      cbind(matrix(0, ncol(s_pen), nfree), s_pen)
    )
  } else {
    nfree <- 0L
  }

  dimnames(z) <- NULL
  dimnames(s_pen) <- NULL
  list(
    X = z,
    S = s_pen,
    # the LEADING columns no entry of the penalty touches. Counted from the
    # matrix rather than from the branch that built it, so the two cannot
    # disagree; cumprod stops at the first column that is penalized.
    unpenalized = as.integer(sum(cumprod(colSums(abs(s_pen)) == 0))),
    blueprint = list(kind = sm@reparam, basis = rp$basis, nfree = nfree,
                     span = sp$params),
    names = nm
  )
}


#' The Weight a Shrunk Null Space Carries
#'
#' @description
#' The penalty `null_space = "shrink"` puts on the directions the roughness
#' matrix does not see: one tenth of what a penalized direction carries.
#'
#' @details
#' It is \pkg{mgcv}'s rule translated rather than a number chosen here.
#' Measured on `mgcv::s(bs = "ts")`, the shrinkage construction leaves the
#' positive eigenvalues of the penalty exactly as they were -- 7708.76 down
#' to 20.82, identical to the unshrunk `"tp"` -- and replaces each zero with
#' 2.082, which is a tenth of the smallest positive one. In Demmler-Reinsch
#' coordinates every penalized direction has eigenvalue exactly 1, so the
#' rule is the single number below.
#'
#' What the tenth buys is measured. Both a weight of 0.1 and a weight of 1
#' let the term leave the model, the fitted values reaching a standard
#' deviation under 1e-6 at a large smoothing parameter. They differ in rate:
#' on a genuinely linear truth at a smoothing parameter of 100, the root
#' mean square error against that truth is 0.0723 at 0.1 and 0.4049 at 1,
#' so the heavier weight destroys a real linear trend at a smoothing
#' parameter chosen to smooth the wiggles.
#'
#' @return A single number.
#'
#' @keywords internal
shrink_weight <- function() 0.1


#' The Coordinates a Smoother's Coefficients Live In
#' @name smoother_reparam
#'
#' @description
#' Applies the constraint and the reparametrization, returning the block, the
#' penalty matrix on it, and the basis object [smoother_apply()] reapplies.
#'
#' @details
#' The three routes differ in what the coefficients mean, and the fit they
#' express is the same span in each.
#'
#' \describe{
#'   \item{`"dr"`}{[dr_basis()] diagonalizes the pencil of the empirical Gram
#'     matrix against the roughness matrix, so the columns are orthogonal
#'     over the data, ordered from the smoothest to the most oscillatory, and
#'     the penalty is the identity. A separable penalty is available under a
#'     diagonal map and not under a general one, so a sparse or heavy-tailed
#'     prior on a smooth is computationally reachable precisely because the
#'     basis is rotated this way.}
#'   \item{`"none"`}{The constrained basis as it stands, with the penalty the
#'     congruence of the roughness matrix. The coefficients are the basis
#'     coefficients, which is what a difference penalty is written on.}
#'   \item{`"orthonorm"`}{The constrained basis rotated so that it satisfies
#'     \eqn{X'X = I} over the observed covariate. The orthonormality is
#'     against the **empirical** measure, which is what makes the design
#'     orthonormal; [orthonorm_basis()] orthonormalizes against the
#'     \eqn{L^2} inner product instead and remains available as an operation
#'     on a basis. It is the **reparametrized part** that is orthonormal:
#'     with `null_space = "keep"` a free column is prepended afterwards and
#'     the whole block is then not orthonormal, while `null_space = "drop"`
#'     gives \eqn{X'X = I} for the block itself, measured at 3.1e-15.}
#' }
#'
#' The three describe the same space, so an unpenalized fit cannot tell them
#' apart: measured at `k = 12` over 300 observations, the fitted values of
#' the three agree to 1.8e-15. What differs is what a coefficient means and
#' therefore what the penalty is.
#'
#' @param sm A [smoother].
#' @param b The basis [smoother_basis()] returned.
#' @param x The covariate.
#' @param g The roughness matrix [smoother_gram()] returned.
#' @param cons The constraint, one column per direction removed, or `NULL`.
#'
#' @return A list of `S`, the penalty on the reparametrized block, and
#'   `basis`, the object whose evaluation **is** that block and which
#'   [smoother_apply()] evaluates at new values.
#'
#' @keywords internal
smoother_reparam <- function(sm, b, x, g, cons) {
  if (identical(sm@reparam, "dr")) {
    d <- dr_basis(b, x, penalty = g,
                  constraints = if (is.null(cons)) NULL else t(cons))
    # the penalty IS the identity there, by the construction; forming the
    # congruence instead would give it only to rounding
    return(list(S = diag(1, d@dimension), basis = d))
  }

  # the constraint in COEFFICIENT space: C'B beta = 0 says the fitted
  # function is orthogonal to C over the observed covariate, which is the
  # same condition dr_basis() imposes in its first step
  cb <- constrain_basis(b, crossprod(cons, basis_eval(b, x)))

  out <- if (identical(sm@reparam, "none")) {
    cb
  } else {
    r <- chol_pd(crossprod(basis_eval(cb, x)))
    if (is.null(r)) {
      stop(paste0(
        "the block is rank deficient over these covariate values, so it",
        " cannot be\n  orthonormalized. Reduce 'k'."
      ), call. = FALSE)
    }
    new_transformed(cb, backsolve(r, diag(nrow(r))),
                    paste0("orthonorm(", cb@basis_name, ")"), "on")
  }

  # THE PENALTY IS READ OFF THE BASIS THE BLOCK COMES FROM, not off the
  # local product. new_transformed() flattens a nested transform, so
  # B %*% (T1 T2) is what a later evaluation computes while (B T1) T2 is
  # what building it here in two steps would; the two agree in exact
  # arithmetic and not in the last bit, and taking both from one object is
  # what makes smoother_apply() reproduce the block exactly.
  tm <- out@transform
  s <- crossprod(tm, g %*% tm)
  list(S = (s + t(s)) / 2, basis = out)
}


#' A Smoother's Block at New Values
#'
#' @description
#' Evaluates the construction recorded by [smoother_build()] at new covariate
#' values. The basis, the constraint and the reparametrization are the ones
#' computed on the original data: nothing is recomputed from `newx`.
#'
#' @param sm A [smoother], the one passed to [smoother_build()].
#' @param blueprint The `blueprint` element of that call's result.
#' @param newx The new covariate values, a numeric vector.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix of `length(newx)` rows and as many columns as the
#'   block the blueprint came from, with no dimnames.
#'
#' @seealso [smoother_build()], which records the blueprint.
#'
#' @examples
#' set.seed(1)
#' x <- sort(runif(150, 0, 1))
#' sm <- bspline_smooth(k = 10)
#' out <- smoother_build(sm, x)
#'
#' # Reapplied at a subset, the rows are the rows of the original block.
#' i <- c(3, 40, 91)
#' max(abs(smoother_apply(sm, out$blueprint, x[i]) - out$X[i, ]))
#' @export
smoother_apply <- S7::new_generic(
  "smoother_apply", "sm",
  function(sm, blueprint, newx, ...) S7::S7_dispatch()
)

#' @name smoother_apply
#' @keywords internal
S7::method(smoother_apply, smoother) <- function(sm, blueprint, newx, ...) {
  if (!is.list(blueprint) || is.null(blueprint$kind)) {
    stop("'blueprint' must be the blueprint element of smoother_build().",
      call. = FALSE
    )
  }
  newx <- check_smoother_x(newx)
  z <- basis_eval(blueprint$basis, newx)
  if (blueprint$nfree > 0L) {
    free <- smoother_span_apply(sm, blueprint$span, newx)
    z <- cbind(free[, seq_len(blueprint$nfree), drop = FALSE], z)
  }
  dimnames(z) <- NULL
  z
}


#' The Names of the Free Columns
#'
#' @description
#' Names the columns `null_space = "keep"` restores: `lin` for the linear
#' one, which is the only one at `order = 2`, and `poly2`, `poly3` and so on
#' for the higher powers an order above 2 leaves unpenalized.
#'
#' @param n How many.
#'
#' @return A character vector of length `n`.
#'
#' @keywords internal
free_names <- function(n) {
  if (n < 1L) return(character(0))
  c("lin", if (n > 1L) paste0("poly", seq.int(2L, n)))
}


#' The Interval a Smoother Expands Over
#'
#' @description
#' Returns the interval a smoother's basis is built on: the endpoints stored
#' on the object when both are given, and otherwise the range of `x` padded
#' by a thousandth of its width. The padding keeps the observed values
#' strictly inside the interval, which is what a basis whose validator
#' requires an open interval needs at its endpoints.
#'
#' @param sm A [smoother].
#' @param x The covariate, a numeric vector.
#'
#' @return A list of two numbers, the lower and upper endpoints.
#'
#' @keywords internal
smoother_interval <- function(sm, x) {
  if (!is.null(sm@lower) && !is.null(sm@upper)) {
    return(list(sm@lower, sm@upper))
  }
  r <- range(x)
  pad <- diff(r) * 0.001 + .Machine$double.eps
  list(
    if (is.null(sm@lower)) r[[1L]] - pad else sm@lower,
    if (is.null(sm@upper)) r[[2L]] + pad else sm@upper
  )
}


#' Check a Smoother's Covariate
#'
#' @description
#' Validates the covariate a smoother is built or reapplied at: one finite
#' numeric value per observation, with no missing values.
#'
#' @param x The covariate.
#'
#' @return `x` as a numeric vector.
#'
#' @keywords internal
check_smoother_x <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 1L) {
    stop("'x' must have at least one value.", call. = FALSE)
  }
  if (anyNA(x) || !all(is.finite(x))) {
    stop("'x' must be finite and must not contain missing values.",
      call. = FALSE
    )
  }
  x
}


#' Check a Whole-Number Argument
#'
#' @description
#' Validates one of the whole-number arguments of a smoother constructor and
#' returns it as an integer.
#'
#' @param v The value.
#' @param nm The argument's name, for the message.
#' @param lo The smallest value admitted.
#' @param why An explanation appended to the message, or `""`.
#'
#' @return `v` as a length-one integer.
#'
#' @keywords internal
check_whole <- function(v, nm, lo, why = "") {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) ||
    v != round(v) || v < lo) {
    stop(sprintf(
      "'%s' must be a whole number of at least %d%s.", nm, lo,
      if (nzchar(why)) paste0(":\n  ", why) else ""
    ), call. = FALSE)
  }
  as.integer(v)
}


#' Check a Smoother's Interval Arguments
#'
#' @description
#' Validates the `lower` and `upper` arguments of a smoother constructor.
#' Each may be `NULL`, meaning the endpoint is read from the data at build.
#'
#' @param lower,upper The endpoints.
#'
#' @return `NULL`, invisibly. Called for the error it signals.
#'
#' @keywords internal
check_interval <- function(lower, upper) {
  for (p in list(list(lower, "lower"), list(upper, "upper"))) {
    v <- p[[1L]]
    if (!is.null(v) && (!is.numeric(v) || length(v) != 1L || !is.finite(v))) {
      stop(sprintf("'%s' must be NULL or a single finite number.", p[[2L]]),
        call. = FALSE
      )
    }
  }
  if (!is.null(lower) && !is.null(upper) && lower >= upper) {
    stop("'lower' must be strictly less than 'upper'.", call. = FALSE)
  }
  invisible(NULL)
}


#' Check What a Smoother Asks For Against What Is Built
#'
#' @description
#' Signals an error for a smoother whose settings contradict each other or
#' reach arithmetic this version does not write. The arguments are on the
#' constructor because they are part of its interface, and an argument
#' accepted and ignored would report a fit of a model the caller did not ask
#' for.
#'
#' @param sm A [smoother].
#'
#' @return `NULL`, invisibly. Called for the error it signals.
#'
#' @keywords internal
check_available <- function(sm) {
  # A FACTORY AND A SHRUNK NULL SPACE CONTRADICT EACH OTHER, and the
  # contradiction is in the arithmetic rather than in the vocabulary.
  # "shrink" is a weight written INSIDE the roughness matrix -- one tenth of
  # what a penalized direction carries, which is a ratio against that
  # matrix's own eigenvalues -- and a factory replaces the matrix with a
  # penalty that has no such eigenvalue to be a tenth of. The free columns
  # would keep a weight in a matrix the fit no longer reads.
  if (!is.null(sm@penalty) && identical(sm@null_space, "shrink")) {
    stop(paste0(
      "'penalty' and null_space = \"shrink\" cannot both be given: the",
      " shrinkage is\n  a weight inside the roughness matrix, and a penalty",
      " factory replaces that\n  matrix. Use null_space = \"keep\" to leave",
      " the unpenalized directions free,\n  or \"drop\" to remove them."
    ), call. = FALSE)
  }
  invisible(NULL)
}


#' Check a Smoother's Penalty Factory
#'
#' @description
#' Validates the `penalty` argument of a smoother constructor: `NULL`, or a
#' function of one argument.
#'
#' @details
#' The check is deliberately weak, and the reason is the dependency graph.
#' \pkg{basis7} sits at the bottom of it and imports \pkg{numericals7} alone,
#' so it cannot name \pkg{penalties7} and cannot ask whether what the
#' function returns is a penalty. It stores the function and never calls it.
#' Whichever layer builds the term calls it, at the coefficient count only
#' the data settle, and checks the result there.
#'
#' @param penalty The value given.
#'
#' @return `penalty`, unchanged.
#'
#' @keywords internal
check_penalty <- function(penalty) {
  if (is.null(penalty)) return(penalty)
  if (is.function(penalty) && length(formals(args(penalty)))) return(penalty)
  stop(paste0(
    "'penalty' must be NULL, or a function of the number of coefficients",
    " giving a\n  penalty. A penalties7 constructor passes bare --",
    " penalty = penalties7::lasso_penalty\n  -- and anything else is",
    " written out, as function(n_coef) my_penalty(n_coef)."
  ), call. = FALSE)
}


#' Check a Smoother's Measure
#'
#' @description
#' Validates the `measure` argument of a smoother constructor: one of the
#' names the package integrates against, or a function of one numeric vector
#' returning one non-negative weight per point.
#'
#' @param measure The value given.
#'
#' @return `measure`, unchanged.
#'
#' @keywords internal
check_measure <- function(measure) {
  if (is.function(measure)) {
    if (length(formals(measure)) < 1L) {
      stop("a 'measure' function must take one argument, the points.",
        call. = FALSE
      )
    }
    return(measure)
  }
  if (!is.character(measure) || length(measure) != 1L ||
    !measure %in% c("lebesgue", "empirical")) {
    stop(paste0(
      "'measure' must be \"lebesgue\", \"empirical\", or a function of one",
      "\n  numeric vector giving a non-negative weight per point."
    ), call. = FALSE)
  }
  measure
}


#' Check a Polynomial Family's Constraint Argument
#'
#' @description
#' Validates `constrain` for a family whose null space is the polynomials of
#' degree below `order`. The constraint must **contain** that null space, so
#' a degree below `order - 1` is rejected.
#'
#' @details
#' The requirement is not a convention. A direction the penalty does not see
#' and the constraint does not remove is neither penalized nor identified,
#' and [dr_basis()] signals an error there: with `order = 3` and a constraint
#' spanning only the constant and the linear function, the quadratic
#' direction is left free and unpenalized. Reported here, where the two
#' numbers were written, rather than several frames down.
#'
#' Constraining **beyond** the null space is legitimate. It buys
#' orthogonality to a parametric term written in the formula. Measured on
#' `y ~ 1 + x + x^2 + s(x)` at `k = 20` and 400 observations, the largest
#' correlation between a column of the block and `x^2` falls from 0.995 at
#' the default to 1.9e-15 at `constrain = 2`, which is exact by
#' construction, and the standard error of the quadratic coefficient falls
#' with it by more than an order of magnitude. The cost is one dimension.
#'
#' @param constrain The value given, `NULL` or a whole number.
#' @param order The penalty's derivative order.
#'
#' @return `constrain` as a length-one integer, or `NULL`.
#'
#' @keywords internal
check_constrain <- function(constrain, order) {
  if (is.null(constrain)) return(NULL)
  constrain <- check_whole(constrain, "constrain", 0L)
  if (constrain < order - 1L) {
    stop(sprintf(paste0(
      "'constrain' (%d) does not contain the null space of an order-%d",
      " penalty.\n  The polynomials of degree below %d are neither",
      " penalized nor identified\n  unless the constraint removes them, so",
      " 'constrain' must be at least %d."
    ), constrain, order, order, order - 1L), call. = FALSE)
  }
  constrain
}


#' @name print.smoother
#' @title Print a Smoother
#'
#' @description
#' Prints the family, the number of basis functions and the four decisions
#' the smoother carries.
#'
#' @param x A [smoother].
#' @param ... Ignored.
#'
#' @return `x`, invisibly. Called for the printing.
#'
#' @examples
#' bspline_smooth(k = 12)
#' bspline_smooth(k = 12, null_space = "drop")
#' @keywords internal
S7::method(print, smoother) <- function(x, ...) {
  cat(sprintf(
    "<%s> %d functions", class(x)[1L], x@dimension
  ))
  if (S7::S7_inherits(x, BsplineSmoother)) {
    cat(sprintf(", degree %d", x@degree))
  }
  cat("\n")
  int <- if (is.null(x@lower) || is.null(x@upper)) {
    "from the data"
  } else {
    sprintf("[%g, %g]", x@lower, x@upper)
  }
  cat(sprintf(
    "  penalty: derivative of order %d, %s measure\n", x@order,
    if (is.character(x@measure)) x@measure else "supplied"
  ))
  cat(sprintf(
    "  null space: %s     coordinates: %s\n", x@null_space, x@reparam
  ))
  cat(sprintf("  interval: %s\n", int))
  invisible(x)
}
