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


#' The Cyclic Smoother Class
#' @name CyclicSmoother
#'
#' @description
#' The class [cyclic_smooth()] returns: a periodic B-spline basis with a
#' roughness penalty and the Demmler-Reinsch reparametrization. It carries
#' `degree` beside the properties every [smoother] has.
#'
#' @inheritParams smoother
#' @param degree The degree of the underlying B-spline, `3` for a cubic.
#'
#' @return An S7 object of class `CyclicSmoother`, inheriting from
#'   [smoother]. Construct one with [cyclic_smooth()].
#'
#' @examples
#' sm <- cyclic_smooth(k = 10)
#' c(class = class(sm)[1], dimension = sm@dimension)
#' @export
CyclicSmoother <- S7::new_class(
  "CyclicSmoother",
  parent = smoother,
  properties = list(degree = S7::class_integer)
)


#' A Cyclic Smoother
#'
#' @description
#' The local periodic smoother: `k` B-spline functions over one period,
#' constrained so that the fit and its first `degree - 1` derivatives take
#' the same value at the two ends of the interval, penalized by the
#' integrated squared derivative of order `order` and rotated to the
#' Demmler-Reinsch coordinates. It is to [fourier_smooth()] what
#' [bspline_smooth()] is to [legendre_smooth()]: a local basis where the
#' other is global, so a feature at one point of the cycle leaves the rest of
#' it alone.
#'
#' @details
#' # The construction
#'
#' A spline of degree \eqn{d} on \eqn{[a, b]} is periodic exactly when
#' \deqn{f^{(j)}(a) = f^{(j)}(b), \qquad j = 0, 1, \ldots, d - 1,}
#' which is \eqn{d} linear conditions on its coefficients. The periodic
#' splines are therefore a subspace of an ordinary spline space, and the
#' whole construction is [constrain_basis()] applied to a B-spline of
#' dimension `k + degree`, whose null space has dimension `k`.
#'
#' Building it this way rather than by folding a widened knot sequence is
#' what makes every quantity exact. The parent basis lives on \eqn{[a, b]},
#' the period itself, so its Gram matrix integrates over one period and the
#' penalty needs no numerical fallback: the roughness matrix is
#' \eqn{T^\top G T} with \eqn{G} the parent's own exact Gram. A folded
#' construction puts its parent on a widened interval, and its Gram then
#' integrates over more than one period.
#'
#' # What it buys over a Fourier basis, and where it buys nothing
#'
#' Both families are periodic, so the choice between them is locality and
#' nothing else, and what it is worth depends entirely on the truth. The
#' measurement sweeps the concentration of a seasonal peak,
#' \eqn{\exp(\kappa \cos t)}, at 400 observations with both bases ten columns
#' wide and both smoothing parameters chosen by \eqn{\mathrm{REML}}, three
#' samples per setting, reporting the root mean square error against the
#' truth:
#'
#' \tabular{lrrr}{
#'   \eqn{\kappa} \tab cyclic \tab fourier \tab ratio \cr
#'   1  \tab 0.0600 \tab 0.0612 \tab 1.02 \cr
#'   3  \tab 0.0652 \tab 0.0655 \tab 1.00 \cr
#'   8  \tab 0.0717 \tab 0.1783 \tab \strong{2.49} \cr
#'   20 \tab 0.3962 \tab 0.6156 \tab 1.55 \cr
#'   50 \tab 0.9282 \tab 1.1168 \tab 1.20
#' }
#'
#' At a low concentration the two are the same fit: \eqn{\exp(\kappa \cos t)}
#' is a von Mises density, whose Fourier coefficients are modified Bessel
#' functions and decay geometrically, so it is the function a global basis is
#' best at. The gain appears where the peak becomes narrow against what ten
#' columns can carry, and falls back again at \eqn{\kappa = 50}, where
#' neither basis of that width represents the spike at all and the comparison
#' stops being about locality. So the family is worth choosing for a
#' localized seasonal feature and is worth nothing for a smooth one, which is
#' the same reading [bspline_smooth()] and [legendre_smooth()] have against
#' each other away from the circle.
#'
#' Measured on a cubic with `k = 9` over \eqn{[0, 1]}: the basis functions
#' agree at the two ends to 5.6e-17 in value, 7.1e-15 in the first
#' derivative and 8.5e-14 in the second, where an ordinary B-spline of the
#' same `k` disagrees by 4.12. The roughness matrix agrees with a fine
#' trapezoid of the second derivatives over one period, the gap falling by
#' exactly 4.00 at each halving of the step, which is the reference's own
#' order of convergence rather than an error of the matrix. Integrating over
#' 99 per cent of the period instead moves that matrix by 167.6 against a
#' size of 2820.4.
#'
#' # What it does not have, and why
#'
#' `constrain` in the form "the polynomials up to degree c" has no reading
#' here, for the reason it has none on [fourier_smooth()]: a non-constant
#' periodic function is never a polynomial, so the null space of the
#' roughness matrix is the constant at every order rather than a space
#' growing with the order. Measured on the pair at orders 1, 2 and 3, the
#' null space is one-dimensional in all three and its function is constant to
#' 2.7e-15.
#'
#' The constant is removed, so a model carrying an intercept spans the level
#' and the smooth carries the shape, and nothing is restored as a free
#' column: a linear column is not periodic, and restoring one is what makes a
#' fit lose the property the basis was chosen for. A block therefore carries
#' `k - 1` columns.
#'
#' # The interval is the period
#'
#' `lower` and `upper` are the ends of one cycle, and they are a property of
#' the problem rather than of the sample: day-of-year data run over
#' \eqn{[0, 365]} whether or not an observation falls on the first day of the
#' year. Left `NULL` they are read from the data, which makes the fitted
#' period the observed range and is almost never what a periodic model means.
#'
#' @param k The number of periodic functions, a whole number of at least 3.
#'   The block carries `k - 1` columns, the constant being removed.
#' @param degree The degree of the underlying B-spline, `3` for a cubic. The
#'   fit matches at the two ends in its value and its first `degree - 1`
#'   derivatives.
#' @param order The order of derivative the penalty integrates, at most
#'   `degree`.
#' @param measure The measure the roughness is integrated against.
#' @param null_space What becomes of the directions the penalty does not see.
#'   A periodic basis restores none, so the settings differ only in what they
#'   refuse.
#' @param reparam The coordinates the coefficients live in.
#' @param penalty `NULL` for the quadratic roughness penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param lower,upper The ends of one period, `NULL` to read them from the
#'   data. See the section above.
#'
#' @return An S7 object of class [CyclicSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [fourier_smooth()] for the global periodic family,
#'   [bspline_smooth()] for the local non-periodic one.
#'
#' @examples
#' set.seed(3)
#' doy <- sort(runif(200, 0, 365))
#' sm <- cyclic_smooth(k = 10, lower = 0, upper = 365)
#' out <- smoother_build(sm, doy)
#' dim(out$X)
#'
#' # the block takes the same value at the two ends of the period
#' ends <- smoother_apply(sm, out$blueprint, c(0, 365))
#' max(abs(ends[1, ] - ends[2, ]))
#'
#' # 'order' may not exceed the degree.
#' try(cyclic_smooth(k = 10, degree = 3, order = 4))
#' @export
cyclic_smooth <- function(k = 10, degree = 3, order = 2,
                          measure = "lebesgue", null_space = "keep",
                          reparam = "dr", penalty = NULL,
                          lower = NULL, upper = NULL) {
  # THREE is the floor the sibling periodic family uses, and for the same
  # reason: the constant is removed, so k = 2 leaves one column and k = 1
  # leaves none. The constraints themselves never lose rank -- measured, the
  # rank is exactly 'degree' at every k from 1 to 12 and every degree from 1
  # to 5 -- so nothing here is a guard against that.
  k <- check_whole(k, "k", 3L)
  degree <- check_whole(degree, "degree", 1L)
  order <- check_whole(order, "order", 1L)
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

  sm <- CyclicSmoother(
    smoother_name = "cyclic",
    dimension = k,
    order = order,
    measure = measure,
    constrain = NULL,
    null_space = null_space,
    reparam = reparam,
    penalty = penalty,
    degree = degree,
    lower = lower,
    upper = upper,
    smoother_params = list()
  )
  check_available(sm)
  sm
}

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, CyclicSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  parent <- bspline_basis(
    lower = int[[1L]], upper = int[[2L]],
    dimension = sm@dimension + sm@degree, degree = sm@degree
  )
  out <- constrain_basis(parent, periodic_constraint(parent, sm@degree))
  # THE DIMENSION IS ASSERTED rather than assumed. constrain_basis() removes
  # one direction per unit of rank, so a constraint that lost rank would
  # return a basis wider than the caller asked for, silently, and the whole
  # block with it.
  if (out@dimension != sm@dimension) {
    stop(sprintf(paste0(
      "the periodicity constraint has rank %d where %d was expected, so the",
      "\n  basis would carry %d functions rather than the %d asked for."
    ), parent@dimension - out@dimension, sm@degree,
    out@dimension, sm@dimension), call. = FALSE)
  }
  out
}

#' @name smoother_span
#' @keywords internal
S7::method(smoother_span, CyclicSmoother) <- function(sm, x, ...) {
  # THE CONSTANT AT EVERY ORDER, and nothing restored -- the same reading as
  # fourier_smooth()'s, and for the same reason: a non-constant periodic
  # function is never a polynomial, so the null space does not grow with the
  # order, and a linear column restored here would not be periodic.
  list(
    constraint = matrix(1, length(x), 1L),
    free = matrix(numeric(0), length(x), 0L),
    params = list(steps = list())
  )
}

#' The Periodicity Constraint of a Spline Basis
#'
#' @description
#' Returns the matrix whose rows say that a function of `basis` and its first
#' `degree - 1` derivatives take the same value at the two ends of the
#' interval: row \eqn{j} is \eqn{B^{(j)}(a) - B^{(j)}(b)} for
#' \eqn{j = 0, \ldots, d - 1}. Its null space is the periodic splines, which
#' is what [cyclic_smooth()] builds on.
#'
#' @details
#' The rows are read from [basis_deriv()] at the two endpoints, order zero
#' included, so the constraint is evaluated by the same arithmetic that
#' evaluates the basis rather than from a formula about knots.
#'
#' @param basis The parent basis, of dimension `k + degree`.
#' @param degree The degree of the spline, one row per derivative order below
#'   it.
#'
#' @return A numeric matrix with `degree` rows and `basis@dimension` columns.
#'
#' @keywords internal
periodic_constraint <- function(basis, degree) {
  ends <- c(basis@lower, basis@upper)
  rows <- lapply(seq_len(degree) - 1L, function(j) {
    m <- basis_deriv(basis, ends, order = j)
    m[1L, ] - m[2L, ]
  })
  do.call(rbind, rows)
}
