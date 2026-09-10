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


#' The P-spline Smoother Class
#' @name PsplineSmoother
#'
#' @description
#' The class [pspline_smooth()] returns: a B-spline basis with a difference
#' penalty on its coefficients. It carries `degree` and `diff` beside the
#' properties every [smoother] has.
#'
#' @inheritParams smoother
#' @param degree The degree of the B-spline.
#' @param diff The order of difference the penalty takes.
#'
#' @return An S7 object of class `PsplineSmoother`, inheriting from
#'   [smoother]. Construct one with [pspline_smooth()].
#'
#' @examples
#' sm <- pspline_smooth(k = 20)
#' c(class = class(sm)[1], dimension = sm@dimension)
#' @export
PsplineSmoother <- S7::new_class(
  "PsplineSmoother",
  parent = smoother,
  properties = list(degree = S7::class_integer, diff = S7::class_integer)
)


#' A P-spline Smoother
#'
#' @description
#' The Eilers-Marx smoother: a B-spline basis of `k` functions over equally
#' spaced knots, penalized by the sum of squared `diff`-th differences of its
#' coefficients rather than by an integrated squared derivative. A rich basis
#' and a cheap penalty, which is the construction's own argument: `k` is
#' chosen large enough not to matter and the smoothing parameter does the
#' rest.
#'
#' @details
#' # The penalty is a functional of the coefficients
#'
#' Where [bspline_smooth()] integrates \eqn{f^{(m)}} against a measure, this
#' penalizes \eqn{\sum_j (\Delta^d c_j)^2}, so the roughness matrix is
#' \eqn{D_d^\top D_d} with \eqn{D_d} the difference operator on the
#' coefficient vector. Nothing is integrated, which is why the family has no
#' `measure` and no `order`: there is no measure to integrate against and no
#' derivative whose order to name.
#'
#' That is also the reason the argument belongs to this family and not to
#' every one. A difference penalty reads the coefficients as an ordered
#' sequence in which neighbours are comparable, which a B-spline's are and a
#' Fourier basis's are not -- there "adjacent" is a sine, a cosine and the
#' next sine, and their difference means nothing.
#'
#' # What it contracts to
#'
#' \eqn{D_d c = 0} exactly when the coefficients are a polynomial of degree
#' below \eqn{d} in their index, and the null space of the roughness matrix
#' has dimension exactly `diff`: measured at `k = 20`, `degree = 3`, it is 1,
#' 2 and 3 at `diff` of 1, 2 and 3.
#'
#' That null space is only APPROXIMATELY the polynomials, and the reason is
#' the boundary knots. Marsden's identity gives
#' \eqn{\sum_j \xi_j B_j(x) = x} with \eqn{\xi_j} the Greville abscissae, so
#' coefficients affine in the index give a straight line exactly where
#' \eqn{\xi_j} is itself affine in \eqn{j} -- which fails at the ends of a
#' clamped sequence, whose boundary knots are repeated. Measured, the
#' \eqn{R^2} of \eqn{\xi_j} against \eqn{j} is 0.9893, 0.9979 and 0.9997 at
#' `k` of 10, 20 and 40, the departure being a fixed number of knots out of
#' `k`; and the functions spanning the null space are the polynomials of
#' degree below `diff` to an \eqn{R^2} of 1.0000000000, 0.9994 and 0.9957.
#'
#' ⚠️ It costs the construction nothing, which is the measurement that
#' matters rather than the one above. [smoother_span()] constrains the block
#' against the exact polynomials, so the Demmler-Reinsch rotation runs on
#' their complement, where the difference penalty is positive definite: the
#' built penalty is the identity to 1.0000000000 on every one of its 23
#' penalized directions, exactly as [bspline_smooth()]'s is, and a strongly
#' penalized fit contracts to a straight line with an \eqn{R^2} against
#' \eqn{(1, x)} of 1.0000000000 at \eqn{\lambda = 10^{10}} for both.
#'
#' # Against the integrated penalty on the same basis
#'
#' The two are different penalties and neither contains the other. Measured
#' at `k = 25`, `degree = 3` over 300 observations, the raw roughness
#' matrices correlate at 0.5075 and their scales differ by four orders --
#' 6 against 2.556e+05 -- the difference operator carrying no factor of the
#' knot spacing.
#'
#' ⚠️ The smoothing parameters nevertheless mean the same thing, and that is
#' the reparametrization doing what it is for. At matched effective degrees
#' of freedom of 5, 8 and 12 the two smoothing parameters stand in a ratio
#' of 1.0, 0.9 and 0.8, not the four orders the raw matrices differ by,
#' because after the Demmler-Reinsch rotation both penalties are the
#' identity. The fitted functions differ by a root mean square of 0.0135,
#' 0.0243 and 0.0180 against a signal whose own standard deviation is
#' 0.7061.
#'
#' @param k The number of basis functions, a whole number greater than
#'   `degree` and greater than `diff`.
#' @param degree The degree of the B-spline, `3` for a cubic.
#' @param diff The order of difference the penalty takes, `2` for the usual
#'   construction. The fit contracts to a polynomial of degree `diff - 1`.
#' @param constrain The directions the smooth is made orthogonal to, `NULL`
#'   for the null space of the penalty.
#' @param null_space What becomes of the directions the penalty does not see.
#' @param reparam The coordinates the coefficients live in.
#' @param penalty `NULL` for the quadratic difference penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param lower,upper The interval, `NULL` to read it from the data.
#'
#' @return An S7 object of class [PsplineSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [bspline_smooth()] for the integrated-derivative penalty on the
#'   same basis.
#'
#' @references
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines
#' and penalties. \emph{Statistical Science} 11(2), 89-121.
#'
#' @examples
#' set.seed(3)
#' x <- sort(runif(200))
#' out <- smoother_build(pspline_smooth(k = 20), x)
#' dim(out$X)
#'
#' # the roughness matrix is a difference operator, so it has no measure
#' try(pspline_smooth(k = 20, measure = "empirical"))
#'
#' # and 'diff' may not reach the number of functions
#' try(pspline_smooth(k = 4, degree = 3, diff = 4))
#' @export
pspline_smooth <- function(k = 20, degree = 3, diff = 2, constrain = NULL,
                           null_space = "keep", reparam = "dr",
                           penalty = NULL, lower = NULL, upper = NULL) {
  k <- check_whole(k, "k", 2L)
  degree <- check_whole(degree, "degree", 1L)
  # base::diff is shadowed by the argument from here on, so the difference
  # matrix below is built with the qualified name
  diff <- check_whole(diff, "diff", 1L)
  if (k < degree + 1L) {
    stop(sprintf(paste0(
      "'k' (%d) is too small for 'degree' (%d): a B-spline basis of degree",
      " m\n  needs at least m + 1 functions."
    ), k, degree), call. = FALSE)
  }
  # AT diff = k the difference matrix has no rows at all and the penalty is
  # the zero matrix; at diff = k - 1 it has one, and the null space is
  # everything but one direction
  if (diff >= k) {
    stop(sprintf(paste0(
      "'diff' (%d) must be smaller than 'k' (%d): the %d-th difference of",
      " %d\n  coefficients has no rows, so the penalty would be the zero",
      " matrix."
    ), diff, k, diff, k), call. = FALSE)
  }
  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("dr", "none", "orthonorm"))
  check_interval(lower, upper)
  penalty <- check_penalty(penalty)
  constrain <- check_constrain(constrain, diff)
  ncon <- if (is.null(constrain)) diff else constrain + 1L
  if (k <= ncon) {
    stop(sprintf(paste0(
      "'k' (%d) leaves nothing to smooth: the constraint removes %d",
      " directions.\n  Raise 'k' above %d."
    ), k, ncon, ncon), call. = FALSE)
  }

  sm <- PsplineSmoother(
    smoother_name = "pspline",
    dimension = k,
    # ORDER RECORDS THE DIFFERENCE ORDER, because it is what the null space
    # is governed by and what smoother_span() reads: the null space of the
    # d-th difference is the polynomials of degree below d, the same one an
    # integrated d-th derivative has. Nothing here integrates anything.
    order = diff,
    measure = "lebesgue",
    constrain = constrain,
    null_space = null_space,
    reparam = reparam,
    penalty = penalty,
    degree = degree,
    diff = diff,
    lower = lower,
    upper = upper,
    smoother_params = list()
  )
  check_available(sm)
  sm
}

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, PsplineSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  bspline_basis(
    lower = int[[1L]], upper = int[[2L]],
    dimension = sm@dimension, degree = sm@degree
  )
}

#' @name smoother_gram
#' @keywords internal
S7::method(smoother_gram, PsplineSmoother) <- function(sm, b, x, ...) {
  # THE ROUGHNESS MATRIX IS NOT A GRAM MATRIX HERE. It is a functional of
  # the coefficients rather than of the functions, so nothing is integrated
  # and the basis enters only through its dimension.
  d <- base::diff(diag(1, b@dimension), differences = sm@diff)
  out <- crossprod(d)
  nm <- basis_colnames(b)
  dimnames(out) <- list(nm, nm)
  out
}


#' The Adaptive Smoother Class
#' @name AdaptiveSmoother
#'
#' @description
#' The class [adaptive_smooth()] returns: a B-spline basis whose difference
#' penalty carries a weight that varies along the covariate, so that one
#' stretch may be smoothed harder than another. It adds `degree`, `diff` and
#' `m` to the properties of [smoother].
#'
#' @inheritParams smoother
#' @param degree The degree of the B-spline pieces.
#' @param diff The order of difference the penalty takes.
#' @param m The number of components the weight profile is built from, which
#'   is also the number of smoothing parameters.
#'
#' @return An S7 object of class `AdaptiveSmoother`, inheriting from
#'   [smoother]. Construct one with [adaptive_smooth()], which validates its
#'   arguments; the class constructor does not.
#'
#' @seealso [adaptive_smooth()], which is the way to build one.
#'
#' @examples
#' sm <- adaptive_smooth(k = 40, m = 5)
#' c(class = class(sm)[1], dimension = sm@dimension, m = sm@m)
#' @export
AdaptiveSmoother <- S7::new_class(
  "AdaptiveSmoother",
  parent = smoother,
  properties = list(
    degree = S7::class_integer,
    diff = S7::class_integer,
    m = S7::class_integer
  )
)


#' An Adaptive Smoother
#'
#' @description
#' A difference penalty whose weight varies along the covariate, so that a
#' function may be smoothed hard where it is quiet and left free where it is
#' not. Where [pspline_smooth()] penalizes every difference alike under one
#' smoothing parameter, this one carries `m` of them and lets the data say
#' how the roughness is distributed.
#'
#' @details
#' # The construction
#'
#' Write \eqn{D} for the matrix of `diff`-th differences of the coefficients.
#' A P-spline penalizes \eqn{\lambda\, c^\top D^\top D\, c}; this one gives
#' each difference a weight of its own,
#' \deqn{c^\top D^\top \mathrm{diag}(w)\, D\, c, \qquad
#'       w = \sum_{i=1}^{m} \lambda_i v_i,}
#' where \eqn{v_1, \ldots, v_m} is a B-spline basis evaluated over the
#' coefficient INDEX. Because \eqn{\mathrm{diag}} is linear the whole penalty
#' is the sum \eqn{\sum_i \lambda_i S_i} with
#' \eqn{S_i = D^\top \mathrm{diag}(v_i) D}, so [smoother_build()] answers with
#' those `m` components and whichever layer places the smooth turns them into
#' one penalty with `m` smoothing parameters. The weight profile is itself a
#' spline, and its coefficients are those smoothing parameters.
#'
#' The index basis is built over the index's own range, so an affine
#' relabelling of the index -- \eqn{1, \ldots, n}, or \eqn{i/n}, or
#' \pkg{mgcv}'s \eqn{i/k} -- gives the identical profile, measured to 1e-16.
#' The construction carries no arbitrary constant.
#'
#' # At equal smoothing parameters it IS a P-spline
#'
#' The weight functions are a B-spline basis, hence a partition of unity, so
#' \eqn{\sum_i v_i = 1} and therefore \eqn{\sum_i S_i = D^\top D} exactly --
#' measured at 8.9e-16 to 2.7e-15 for `m` from 2 to 12. Holding every
#' \eqn{\lambda_i} at one value gives the P-spline penalty at that value, so
#' [pspline_smooth()] is the interior point of this family rather than a
#' different construction, and the extra freedom is spent only where the data
#' pay for it.
#'
#' # What it buys, and what it costs
#'
#' Measured against a single-lambda P-spline through \pkg{mgcv}'s REML, which
#' shares no code with this package, on 400 observations at `k = 40` and
#' eight seeds: on a truth of variable roughness -- a sine with a narrow bump
#' -- the adaptive wins on 8 seeds of 8, at a median root mean square error
#' of 0.0357 against 0.0443, and does so at FEWER effective degrees of
#' freedom, 17.4 against 23.6. On a truth of constant roughness it wins on 0
#' seeds of 8, 0.0315 against 0.0310: about 1.6 per cent worse, which is what
#' the extra smoothing parameters cost where there is nothing to adapt to.
#' That second measurement is the control, without which the first would show
#' only that more parameters fit better.
#'
#' Where the gain comes from is not where it is first looked for. Split by
#' region on one sample, the adaptive is better on the quiet left (0.0242
#' against 0.0305) and on the smooth right (0.0437 against 0.0476) and
#' slightly WORSE at the feature itself (0.0553 against 0.0503). What it buys
#' is not a sharper peak but less noise chasing where the function is quiet.
#'
#' # The coordinates
#'
#' `reparam = "dr"` is rejected, and by construction rather than by choice:
#' Demmler-Reinsch diagonalizes the pencil of the Gram matrix against a
#' single penalty, and here there are `m` of them, so a rotation making one
#' component the identity leaves the others arbitrary. The default is
#' `"none"`, which the measurement prefers on both axes it can be judged on.
#' At a spread of smoothing parameters an adaptive really reaches -- 1.3e8,
#' measured on \pkg{mgcv} -- the condition number of the system solved is
#' 3.4e3 in the raw coordinates against 4.7e4 orthonormalized, and the
#' components' scales, which decide whether the `m` smoothing parameters are
#' comparable with one another, spread by a factor of 1.7 against 2.8.
#'
#' # The rank is the family's
#'
#' Each \eqn{S_i} is nearly all null space, its weight vanishing off the
#' support of its own weight function: measured at `k = 40`, `m = 5`, the
#' five components have null dimensions 20, 1, 1, 1 and 19 out of 38. The
#' null space of the SUM is the intersection of theirs and does not move with
#' the smoothing parameters, which is why the rank must be read from the
#' components rather than from the assembled \eqn{S(\lambda)}. Measured, a
#' rank counted off the assembled matrix reads 38, 38, 26 and 19 as one
#' parameter is raised through 1, 1e6, 1e12 and 1e24, where the family's own
#' is 38 throughout.
#'
#' @param k The number of basis functions. A rich basis is the point of the
#'   construction, so the default is larger than [pspline_smooth()]'s.
#' @param degree The degree of the B-spline, `3` for a cubic.
#' @param diff The order of difference the penalty takes. The fit contracts
#'   to a polynomial of degree `diff - 1`.
#' @param m The number of weight components, and so the number of smoothing
#'   parameters. \pkg{mgcv} calls it `m` as well, and takes the same default.
#'   Two gives one weight rising and one falling across the index; more give
#'   a finer profile at the price of a smoothing parameter each.
#' @param constrain The directions the smooth is made orthogonal to, `NULL`
#'   for the null space of the penalty.
#' @param null_space What becomes of the directions the penalty does not see.
#'   `"shrink"` is rejected here; see the note below.
#' @param reparam The coordinates the coefficients live in, `"none"` or
#'   `"orthonorm"`. `"dr"` is rejected.
#' @param penalty Rejected here, and accepted only to say so: the penalty of
#'   this family is the sum of its components, and a factory replaces it with
#'   one penalty over the coefficients, which is [pspline_smooth()] with that
#'   factory.
#' @param lower,upper The interval, `NULL` to read it from the data.
#'
#' @return An S7 object of class [AdaptiveSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [pspline_smooth()], which is this family at equal smoothing
#'   parameters, and [smoother_build()], which returns the components.
#'
#' @references
#' Ruppert, D. and Carroll, R. J. (2000). Spatially-adaptive penalties for
#' spline fitting. \emph{Australian and New Zealand Journal of Statistics}
#' 42(2), 205-223.
#'
#' Krivobokova, T., Crainiceanu, C. M. and Kauermann, G. (2008). Fast
#' adaptive penalized splines. \emph{Journal of Computational and Graphical
#' Statistics} 17(1), 1-20.
#'
#' @examples
#' set.seed(4)
#' x <- sort(runif(300))
#' out <- smoother_build(adaptive_smooth(k = 30, m = 4), x)
#'
#' # the penalty comes back as one component per smoothing parameter
#' length(out$S)
#' dim(out$S[[1]])
#'
#' # and at equal smoothing parameters their sum is the P-spline penalty
#' ps <- smoother_build(pspline_smooth(k = 30, reparam = "none"), x)
#' max(abs(Reduce(`+`, out$S) - ps$S))
#'
#' # the coordinates cannot be Demmler-Reinsch: there are several pencils
#' try(adaptive_smooth(k = 30, m = 4, reparam = "dr"))
#' @export
adaptive_smooth <- function(k = 40, degree = 3, diff = 2, m = 5,
                            constrain = NULL, null_space = "keep",
                            reparam = "none", penalty = NULL,
                            lower = NULL, upper = NULL) {
  k <- check_whole(k, "k", 2L)
  degree <- check_whole(degree, "degree", 1L)
  # base::diff is shadowed by the argument from here on, as in
  # pspline_smooth(), so the difference matrix is built with the qualified
  # name where it is built
  diff <- check_whole(diff, "diff", 1L)
  # ONE component is a P-spline and is that family's business, not a
  # degenerate case of this one
  m <- check_whole(m, "m", 2L,
                   "a single component is pspline_smooth(), on the same basis.")
  if (k < degree + 1L) {
    stop(sprintf(paste0(
      "'k' (%d) is too small for 'degree' (%d): a B-spline basis of degree",
      " m\n  needs at least m + 1 functions."
    ), k, degree), call. = FALSE)
  }
  if (diff >= k) {
    stop(sprintf(paste0(
      "'diff' (%d) must be smaller than 'k' (%d): the %d-th difference of",
      " %d\n  coefficients has no rows, so the penalty would be the zero",
      " matrix."
    ), diff, k, diff, k), call. = FALSE)
  }
  # THE WEIGHT PROFILE LIVES ON THE DIFFERENCES, of which there are
  # k - diff, so it cannot hold more functions than there are of them. It is
  # mgcv's own condition, which stops at the same place.
  if (m >= k - diff) {
    stop(sprintf(paste0(
      "'m' (%d) is too large for 'k' (%d) at diff = %d: the weight profile",
      " is\n  carried on the %d differences of the coefficients, so it must",
      " hold fewer\n  than %d functions. Raise 'k' or lower 'm'."
    ), m, k, diff, k - diff, k - diff), call. = FALSE)
  }
  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("none", "orthonorm", "dr"))
  if (identical(reparam, "dr")) {
    stop(paste0(
      "reparam = \"dr\" is not available for an adaptive smoother:",
      " Demmler-Reinsch\n  diagonalizes the pencil of the Gram matrix",
      " against ONE penalty, and this family\n  carries several, so a",
      " rotation making one of them the identity leaves the\n  others",
      " arbitrary. Use \"none\" (the default) or \"orthonorm\"."
    ), call. = FALSE)
  }
  # THE SHRINKAGE IS A RATIO AGAINST ONE MATRIX'S EIGENVALUES -- a tenth of
  # what a penalized direction carries -- and with several components the
  # free direction would be shrunk once per component, so each smoothing
  # parameter would silently also govern the polynomial part. Refused rather
  # than given a meaning that nothing calibrates.
  if (identical(null_space, "shrink")) {
    stop(paste0(
      "null_space = \"shrink\" is not available for an adaptive smoother:",
      " the shrinkage\n  is one tenth of what a penalized direction",
      " carries, which is a ratio against a\n  single roughness matrix, and",
      " here the penalty is a sum of components. Use\n  \"keep\" to leave",
      " the unpenalized directions free, or \"drop\" to remove them."
    ), call. = FALSE)
  }
  if (!is.null(penalty)) {
    stop(paste0(
      "'penalty' is not available for an adaptive smoother: its penalty IS",
      " the sum\n  of its components, and a factory replaces that sum with",
      " one penalty over the\n  coefficients. That model is",
      " pspline_smooth(penalty = ...), on the same basis."
    ), call. = FALSE)
  }
  check_interval(lower, upper)
  constrain <- check_constrain(constrain, diff)
  ncon <- if (is.null(constrain)) diff else constrain + 1L
  if (k <= ncon) {
    stop(sprintf(paste0(
      "'k' (%d) leaves nothing to smooth: the constraint removes %d",
      " directions.\n  Raise 'k' above %d."
    ), k, ncon, ncon), call. = FALSE)
  }

  sm <- AdaptiveSmoother(
    smoother_name = "adaptive",
    dimension = k,
    # ORDER RECORDS THE DIFFERENCE ORDER, as it does for pspline_smooth():
    # it is what smoother_span() reads to build the constraint, and the null
    # space of the SUM of the components is the null space of the difference
    # operator, the weights being a partition of unity.
    order = diff,
    measure = "lebesgue",
    constrain = constrain,
    null_space = null_space,
    reparam = reparam,
    penalty = NULL,
    degree = degree,
    diff = diff,
    m = m,
    lower = lower,
    upper = upper,
    smoother_params = list()
  )
  check_available(sm)
  sm
}

#' @name smoother_basis
#' @keywords internal
S7::method(smoother_basis, AdaptiveSmoother) <- function(sm, x, ...) {
  int <- smoother_interval(sm, x)
  bspline_basis(
    lower = int[[1L]], upper = int[[2L]],
    dimension = sm@dimension, degree = sm@degree
  )
}

#' @name smoother_gram
#' @keywords internal
S7::method(smoother_gram, AdaptiveSmoother) <- function(sm, b, x, ...) {
  # ONE COMPONENT PER SMOOTHING PARAMETER. Nothing is integrated here either:
  # the penalty is a functional of the coefficients, and the basis enters
  # only through its dimension.
  d <- base::diff(diag(1, b@dimension), differences = sm@diff)
  # THE WEIGHT PROFILE IS A SPLINE OVER THE COEFFICIENT INDEX, built over the
  # index's own range, so any affine relabelling of the index gives the same
  # profile to the last bit. Its degree falls with m only where m is too
  # small to carry a cubic: at m = 2 the profile is one weight rising and one
  # falling, which is a partition of unity where mgcv's own two-component
  # case, cbind(1, index), is not.
  idx <- seq_len(nrow(d))
  v <- basis_eval(
    bspline_basis(min(idx), max(idx), dimension = sm@m,
                  degree = min(3L, sm@m - 1L)),
    idx
  )
  nm <- basis_colnames(b)
  lapply(seq_len(sm@m), function(i) {
    out <- crossprod(d, as.numeric(v[, i]) * d)
    out <- (out + t(out)) / 2
    dimnames(out) <- list(nm, nm)
    out
  })
}
