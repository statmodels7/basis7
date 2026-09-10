#' The Fourier Smoother Class
#' @name FourierSmoother
#'
#' @description
#' The class [fourier_smooth()] returns: a Fourier basis with a roughness
#' penalty and the Demmler-Reinsch reparametrization. It has no `degree` and
#' no `constrain`, neither of which a periodic family has a reading for.
#'
#' @inheritParams smoother
#' @param omega The period. `NULL` takes the width of the interval.
#'
#' @return An S7 object of class `FourierSmoother`, inheriting from
#'   [smoother]. Construct one with [fourier_smooth()].
#'
#' @seealso [fourier_smooth()], which is the way to build one.
#'
#' @examples
#' sm <- fourier_smooth(k = 9)
#' c(class = class(sm)[1], dimension = sm@dimension)
#' @export
FourierSmoother <- S7::new_class(
  "FourierSmoother",
  parent = smoother,
  properties = list(omega = S7::class_any)
)


#' A Fourier Smoother
#'
#' @description
#' The periodic smoother: `k` Fourier functions over an interval, penalized
#' by the integrated squared derivative of order `order`, rotated to the
#' Demmler-Reinsch coordinates. Every column of the block is periodic, so a
#' fit built on it takes the same value at the two ends of the interval.
#'
#' @details
#' # What it does not have, and why
#'
#' `degree` is a B-spline's, and a Fourier basis has none.
#'
#' `constrain` in the form "the polynomials up to degree c" has no reading
#' here: a periodic basis contains no linear function, so its penalty's null
#' space is the **constant at every order**, not a space growing with the
#' order. Measured on the Gram matrix of a nine-function Fourier basis at
#' orders 1, 2 and 3, the null function has a standard deviation of exactly
#' zero, which is to say it is constant, and the null space is
#' one-dimensional in all three.
#'
#' The constant is removed, so a model carrying an intercept spans the level
#' and the smooth carries the shape. Nothing is restored as a free column,
#' which is why `null_space` has no effect here and `k` functions give
#' `k - 1` columns.
#'
#' # The interval is the period
#'
#' A periodic basis needs its interval fixed by the modeller, not read from
#' the data: the period is a fact about the covariate, and the range of one
#' sample is not it. Give `lower` and `upper`, as in
#' `fourier_smooth(k = 9, lower = 0, upper = 365)` for a day of the year.
#' Without them the interval is the padded range of the data, which makes
#' the basis periodic on an interval nothing else knows about.
#'
#' @param k The number of basis functions, an **odd** whole number of at
#'   least 3: a Fourier basis holds a constant plus complete sine-cosine
#'   pairs.
#' @param order The order of derivative the penalty integrates.
#' @param measure The measure the roughness is integrated against.
#' @param null_space What becomes of the null space. It is the constant
#'   alone here, and the constant is removed whatever this says, so the
#'   argument is accepted for symmetry with the other families and changes
#'   nothing.
#' @param reparam The coordinates the coefficients live in.
#' @param penalty `NULL` for the quadratic roughness penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param omega The period, `NULL` for the width of the interval.
#' @param lower,upper The interval, which for a periodic basis is the period.
#'   See the section above: give both.
#'
#' @return An S7 object of class [FourierSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [smoother_build()] for what it produces at data,
#'   [bspline_smooth()] for the non-periodic family, [fourier_basis()] for
#'   the basis alone.
#'
#' @examples
#' # A periodic covariate, and a truth that is periodic on [0, 1].
#' set.seed(3)
#' x <- sort(runif(300))
#' f <- function(u) sin(2 * pi * u) + 0.4 * cos(4 * pi * u)
#' y <- f(x) + rnorm(300, sd = 0.2)
#'
#' sm <- fourier_smooth(k = 9, lower = 0, upper = 1)
#' out <- smoother_build(sm, x)
#' dim(out$X)
#'
#' # THE FIT IS PERIODIC, which is the property the basis was chosen for.
#' b <- solve(crossprod(out$X) + 0.01 * out$S, crossprod(out$X, y))
#' ends <- smoother_apply(sm, out$blueprint, c(0, 1)) %*% b
#' format(diff(as.vector(ends)), digits = 3)
#' @export
fourier_smooth <- function(k = 9, order = 2, measure = "lebesgue",
                           null_space = "keep", reparam = "dr",
                           penalty = NULL, omega = NULL,
                           lower = NULL, upper = NULL) {
  k <- check_whole(k, "k", 3L)
  if (k %% 2L == 0L) {
    stop(sprintf(paste0(
      "'k' must be odd: a Fourier basis holds a constant plus complete",
      " sine-cosine\n  pairs, so %d would leave half a pair. Use %d or %d."
    ), k, k - 1L, k + 1L), call. = FALSE)
  }
  order <- check_whole(order, "order", 1L)
  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("dr", "none", "orthonorm"))
  check_interval(lower, upper)
  measure <- check_measure(measure)
  penalty <- check_penalty(penalty)
  if (!is.null(omega) && (!is.numeric(omega) || length(omega) != 1L ||
    !is.finite(omega) || omega <= 0)) {
    stop("'omega' must be NULL or a single positive number.", call. = FALSE)
  }

  sm <- FourierSmoother(
    smoother_name = "fourier",
    dimension = k,
    order = order,
    measure = measure,
    constrain = NULL,
    null_space = null_space,
    reparam = reparam,
    penalty = penalty,
    omega = omega,
    lower = lower,
    upper = upper,
    smoother_params = list()
  )
  check_available(sm)
  sm
}

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, FourierSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  fourier_basis(
    lower = int[[1L]], upper = int[[2L]],
    dimension = sm@dimension, omega = sm@omega
  )
}

#' @name smoother_span
#' @keywords internal
S7::method(smoother_span, FourierSmoother) <- function(sm, x, ...) {
  # THE CONSTANT AT EVERY ORDER, and nothing restored. A periodic basis
  # contains no linear function, so the null space does not grow with the
  # order the way a polynomial family's does; and the constant is removed
  # because a model carrying an intercept already spans it. Restoring a
  # non-periodic column here is what makes a fit on a periodic basis lose
  # the property the basis was chosen for.
  list(
    constraint = matrix(1, length(x), 1L),
    free = matrix(numeric(0), length(x), 0L),
    params = list(steps = list())
  )
}


#' The Legendre Smoother Class
#' @name LegendreSmoother
#'
#' @description
#' The class [legendre_smooth()] returns: an orthogonal polynomial basis with
#' a roughness penalty and the Demmler-Reinsch reparametrization. Its null
#' space is the polynomials of degree below `order`, as a B-spline's is, so
#' it inherits the default [smoother_span()].
#'
#' @inheritParams smoother
#'
#' @return An S7 object of class `LegendreSmoother`, inheriting from
#'   [smoother]. Construct one with [legendre_smooth()].
#'
#' @seealso [legendre_smooth()], which is the way to build one.
#'
#' @examples
#' sm <- legendre_smooth(k = 8)
#' c(class = class(sm)[1], dimension = sm@dimension)
#' @export
LegendreSmoother <- S7::new_class("LegendreSmoother", parent = smoother)


#' A Legendre Smoother
#'
#' @description
#' The orthogonal polynomial smoother: `k` shifted Legendre polynomials over
#' an interval, penalized by the integrated squared derivative of order
#' `order`, rotated to the Demmler-Reinsch coordinates. It is a global basis
#' where a B-spline is local: every function is supported on the whole
#' interval, so a feature at one end moves the fit at the other.
#'
#' @details
#' The null space of the roughness matrix is the polynomials of degree below
#' `order`, exactly as for a B-spline, so `order` and `constrain` mean what
#' they mean there and `null_space = "keep"` restores `order - 1` free
#' columns. Measured, the null space of the order-2 Gram matrix of an
#' eight-function Legendre basis is spanned by the constant and the linear
#' function to an R-squared of 1.0000000000.
#'
#' @param k The number of polynomials, a whole number of at least 2. The
#'   highest degree is `k - 1`.
#' @param order The order of derivative the penalty integrates, at most
#'   `k - 1`.
#' @param measure The measure the roughness is integrated against.
#' @param constrain The directions the smooth is made orthogonal to, `NULL`
#'   for the null space of the penalty.
#' @param null_space What becomes of the directions the penalty does not see.
#' @param reparam The coordinates the coefficients live in.
#' @param penalty `NULL` for the quadratic roughness penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param lower,upper The interval, `NULL` to read it from the data.
#'
#' @return An S7 object of class [LegendreSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [bspline_smooth()] for the local family, [poly_basis()] for the
#'   basis alone.
#'
#' @examples
#' set.seed(3)
#' x <- sort(runif(200, -1, 1))
#' out <- smoother_build(legendre_smooth(k = 8), x)
#' dim(out$X)
#' out$names
#'
#' # 'order' may not exceed what the basis can differentiate away.
#' try(legendre_smooth(k = 3, order = 4))
#' @export
legendre_smooth <- function(k = 8, order = 2, measure = "lebesgue",
                            constrain = NULL, null_space = "keep",
                            reparam = "dr", penalty = NULL,
                            lower = NULL, upper = NULL) {
  k <- check_whole(k, "k", 2L)
  order <- check_whole(order, "order", 1L)
  # ABOVE THE HIGHEST DEGREE the Gram matrix is identically zero, so the
  # penalty would penalize nothing
  if (order > k - 1L) {
    stop(sprintf(paste0(
      "'order' (%d) exceeds the highest degree the basis carries (%d): the",
      "\n  derivative of that order is zero, so the penalty would be the",
      " zero matrix."
    ), order, k - 1L), call. = FALSE)
  }
  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("dr", "none", "orthonorm"))
  check_interval(lower, upper)
  measure <- check_measure(measure)
  penalty <- check_penalty(penalty)
  constrain <- check_constrain(constrain, order)
  ncon <- if (is.null(constrain)) order else constrain + 1L
  if (k <= ncon) {
    stop(sprintf(paste0(
      "'k' (%d) leaves nothing to smooth: the constraint removes %d",
      " directions.\n  Raise 'k' above %d."
    ), k, ncon, ncon), call. = FALSE)
  }

  sm <- LegendreSmoother(
    smoother_name = "legendre",
    dimension = k,
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

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, LegendreSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  poly_basis(lower = int[[1L]], upper = int[[2L]], dimension = sm@dimension)
}
