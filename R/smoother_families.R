#' The Fourier Smoother Class
#' @name FourierSmoother
#'
#' @description
#' The class [fourier_smooth()] returns: a Fourier basis with a roughness
#' penalty and the Demmler-Reinsch reparametrization. It has no `degree` and
#' no `constrain`, which have no meaning for a periodic family.
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
#' by a linear differential operator and rotated to the Demmler-Reinsch
#' coordinates. Every column of the block is periodic, so a fit built on it
#' takes the same value at the two ends of the interval.
#'
#' @details
#' # The penalty is the harmonic acceleration operator
#'
#' The default `order` is [harmonic_operator()] and not a derivative, and
#' the period that it uses is the width of the interval. The two differ in
#' what a strongly penalized fit contracts **to**. On a periodic basis a
#' derivative penalty has the constant as its only null direction, so the
#' fit contracts to a constant; the harmonic operator leaves the level and
#' the fundamental cycle unpenalized and penalizes every higher harmonic, so
#' the fit contracts to a sinusoid.
#'
#' The default null space is `"shrink"` for this family and `"keep"` for the
#' others. With `"keep"` the sine and cosine of the fundamental are free
#' columns, so they are spent whether or not the covariate carries a cycle;
#' with `"shrink"` they are penalized lightly and can leave the model. The
#' first suits a cycle known to be present and estimated, the second a cycle
#' whose presence is a hypothesis.
#'
#' `order = 2` gives the integrated squared second derivative.
#'
#' # Arguments the family does not take
#'
#' `degree` is a B-spline's, and a Fourier basis has none.
#'
#' `constrain` in the form "the polynomials up to degree c" has no reading
#' here: a periodic basis contains no linear function, so a derivative
#' penalty's null space is the **constant at every order**, not a space
#' growing with the order: the null space of the Gram matrix of a Fourier
#' basis is one-dimensional at every order, spanned by the constant.
#'
#' # Accepted operators
#'
#' The constant is always removed, including when the operator penalizes it,
#' as [oscillator_operator()] does, so a model carrying an intercept spans
#' the level and the smooth carries the shape. Every **other** function of
#' the operator's null space is restored as a free column, and it must be
#' periodic on the period of the basis, which means a pure sine or cosine at
#' a whole multiple of the fundamental frequency. An operator whose null
#' space holds \eqn{t}, or \eqn{e^{at}}, or a frequency that is not a
#' multiple, is rejected, because such a column would make the fit
#' non-periodic.
#'
#' So `harmonic_operator()` and `oscillator_operator()` are accepted at any
#' number of harmonics that leaves at least one function to penalize, which
#' [fourier_smooth()] checks at construction; `deriv_operator(m)` is
#' accepted and restores nothing; and `deriv_operator(2) *
#' oscillator_operator(p)`, with a period `p`, is rejected because its null
#' space holds \eqn{t}, which suits a B-spline and not a periodic basis.
#'
#' # The interval is the period
#'
#' The interval of a periodic basis is fixed by the modeler and not read
#' from the data, because the period is a property of the covariate and the
#' range of one sample does not determine it. `lower` and `upper` give it, as
#' in `fourier_smooth(k = 9, lower = 0, upper = 365)` for a day of the year.
#' Without them the interval is the range of the data padded by a thousandth
#' of its width, and the basis is periodic on that interval.
#'
#' @param k The number of basis functions, an **odd** whole number of at
#'   least 3: a Fourier basis holds a constant plus complete sine-cosine
#'   pairs.
#' @param order What the penalty measures: a [LinearOperator], or a whole
#'   number `m` as the shorthand for `deriv_operator(m)`. The default is the
#'   harmonic acceleration operator at the period of the interval.
#' @param measure The measure the roughness is integrated against.
#' @param null_space What happens to the directions that the penalty does
#'   not see, other than the constant, which is always removed. `NULL` takes
#'   `"shrink"`, or `"keep"` where a `penalty` factory is given, the two
#'   being incompatible. See the section above.
#' @param reparam The coordinates in which the coefficients are expressed.
#' @param penalty `NULL` for the quadratic roughness penalty, or a factory
#'   building a penalty from a coefficient count. See the section on the
#'   smoother's own page.
#' @param omega The period **of the basis**, `NULL` for the width of the
#'   interval. The period of the operator is carried by
#'   [harmonic_operator()].
#' @param lower,upper The interval, which for a periodic basis is the period.
#'   See the section above: give both.
#'
#' @return An S7 object of class [FourierSmoother], inheriting from
#'   [smoother].
#'
#' @seealso [harmonic_operator()] for the default penalty,
#'   [smoother_build()] for what it produces at data, [bspline_smooth()] for
#'   the non-periodic family, [fourier_basis()] for the basis alone.
#'
#' @examples
#' # A periodic covariate, and a truth that is periodic on [0, 1].
#' set.seed(3)
#' x <- sort(runif(300))
#' f <- function(u) sin(2 * pi * u) + 0.4 * cos(4 * pi * u)
#' y <- f(x) + rnorm(300, sd = 0.2)
#'
#' sm <- fourier_smooth(k = 9, lower = 0, upper = 1)
#' sm
#' out <- smoother_build(sm, x)
#' dim(out$X)
#'
#' # THE FIT IS PERIODIC, which is the property the basis was chosen for.
#' b <- solve(crossprod(out$X) + 0.01 * out$S, crossprod(out$X, y))
#' ends <- smoother_apply(sm, out$blueprint, c(0, 1)) %*% b
#' format(diff(as.vector(ends)), digits = 3)
#'
#' # THE FUNDAMENTAL IS WHAT A STRONG PENALTY LEAVES. With the null space
#' # kept it is free, and the heavily penalized fit is a pure sinusoid.
#' sk <- fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep")
#' ok <- smoother_build(sk, x)
#' ok$unpenalized
#' bk <- solve(crossprod(ok$X) + 1e8 * ok$S, crossprod(ok$X, y - mean(y)))
#' fv <- as.vector(ok$X %*% bk)
#' round(sqrt(mean(resid(lm(fv ~ sin(2 * pi * x) + cos(2 * pi * x)))^2)), 8)
#'
#' # The old construction, in one argument.
#' fourier_smooth(k = 9, order = 2, lower = 0, upper = 1)
#'
#' # Two harmonics left alone instead of one.
#' fourier_smooth(k = 13, lower = 0, upper = 365,
#'                order = harmonic_operator(harmonics = 2))
#'
#' # An operator whose null space a periodic basis cannot carry.
#' try(smoother_build(
#'   fourier_smooth(k = 9, lower = 0, upper = 1,
#'                  order = deriv_operator(2) * oscillator_operator(1)), x))
#' @export
fourier_smooth <- function(k = 9, order = harmonic_operator(),
                           measure = "lebesgue",
                           null_space = NULL, reparam = "dr",
                           penalty = NULL, omega = NULL,
                           lower = NULL, upper = NULL) {
  # A PENALTY FACTORY AND "shrink" ARE REFUSED TOGETHER by check_available(),
  # so the family's own default cannot be "shrink" unconditionally: it would
  # turn fourier_smooth(penalty = f), which has always worked, into an
  # error. The two arguments genuinely interact and the resolution says how.
  if (is.null(null_space)) {
    null_space <- if (is.null(penalty)) "shrink" else "keep"
  }
  k <- check_whole(k, "k", 3L)
  if (k %% 2L == 0L) {
    stop(sprintf(paste0(
      "'k' must be odd: a Fourier basis holds a constant plus complete",
      " sine-cosine\n  pairs, so %d would leave half a pair. Use %d or %d."
    ), k, k - 1L, k + 1L), call. = FALSE)
  }
  order <- as_operator(order)
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
  check_periodic_size(order, k)
  sm
}


#' Check That a Periodic Basis Is Wider Than What the Operator Removes
#'
#' @description
#' Signals an error at construction when the constant and the functions of
#' the operator's null space take all `k` functions of the basis, leaving
#' none to penalize. The count is the order of the operator, plus one where
#' the constant is not in its null space, as for [oscillator_operator()].
#'
#' @param op A [LinearOperator], resolved or not.
#' @param k The number of basis functions.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @keywords internal
check_periodic_size <- function(op, k) {
  if (is_deriv_operator(op)) return(invisible(NULL))
  kind <- op@operator_params$kind
  const_in_null <- if (identical(kind, "harmonic")) {
    TRUE
  } else if (identical(kind, "oscillator")) {
    FALSE
  } else {
    isTRUE(op@weights[[1L]] == 0)
  }
  removed <- operator_order(op) + as.integer(!const_in_null)
  if (removed >= k) {
    kmin <- removed + 1L
    if (kmin %% 2L == 0L) kmin <- kmin + 1L
    stop(sprintf(paste0(
      "the constant and the operator's null space take %d functions, and",
      " the basis has %d,\n  which leaves none to penalize. Use 'k' of at",
      " least %d."
    ), removed, k, kmin), call. = FALSE)
  }
  invisible(NULL)
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
  op <- smoother_operator(sm, x)
  if (is_deriv_operator(op)) {
    # THE CONSTANT AT EVERY ORDER, and nothing restored. A periodic basis
    # contains no linear function, so a derivative penalty's null space does
    # not grow with the order the way a polynomial family's does; and the
    # constant is removed because a model carrying an intercept already
    # spans it.
    return(list(
      constraint = matrix(1, length(x), 1L),
      free = matrix(numeric(0), length(x), 0L),
      params = list(steps = list())
    ))
  }
  check_periodic_null(sm, op, x)
  sp <- operator_span(op, x)
  # the constant is removed even where the operator penalizes it, as
  # oscillator_operator() does: the basis contains it, and a model carrying
  # an intercept would otherwise hold it twice
  if (all(sp$params$keep)) {
    sp$constraint <- cbind(1, sp$constraint)
  }
  sp
}


#' Check That the Null Space of an Operator Is Periodic
#'
#' @description
#' Every function of the operator's null space other than the constant is
#' restored as a free column of a periodic block, so each must itself be
#' periodic on the basis's period: a pure sine or cosine at a whole multiple
#' of the fundamental frequency, with no power of \eqn{t} and no exponential
#' in front of it.
#'
#' @details
#' Restoring anything else gives a block whose columns are not all periodic,
#' and a fit on it no longer takes the same value at the two ends of the
#' interval, the property for which a Fourier basis is chosen. The check is
#' made on the null space of the operator, read analytically, and not on the
#' rank of the assembled penalty, which depends on the tolerance.
#'
#' @param sm A [FourierSmoother].
#' @param op A [LinearOperator], with its period resolved.
#' @param x The covariate, from which the interval is taken where the
#'   smoother does not fix it.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @keywords internal
check_periodic_null <- function(sm, op, x) {
  int <- smoother_interval(sm, x)
  period <- if (is.null(sm@omega)) int[[2L]] - int[[1L]] else sm@omega
  nu <- 2 * pi / period
  nl <- operator_null(op)
  bad <- nl$degree > 0L | abs(nl$rate) > 1e-8 |
    abs(nl$freq / nu - round(nl$freq / nu)) > 1e-6
  if (any(bad)) {
    stop(sprintf(paste0(
      "a periodic basis cannot carry every function of this operator's null",
      " space.\n  It would restore %s as free columns, and a column that is",
      " not periodic on\n  the basis's period costs the fit the property",
      " the basis was chosen for.\n  Use harmonic_operator() or",
      " oscillator_operator() at this period, a whole\n  'order', or a",
      " non-periodic family such as bspline_smooth()."
    ), paste(nl$label[bad], collapse = ", ")), call. = FALSE)
  }
  invisible(NULL)
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
#' columns.
#'
#' @param k The number of polynomials, a whole number of at least 2. The
#'   highest degree is `k - 1`.
#' @param order What the penalty measures: a [LinearOperator] from
#'   [deriv_operator()], [harmonic_operator()], [oscillator_operator()] or
#'   [linear_operator()], or a whole number `m` as the shorthand for
#'   `deriv_operator(m)`. It determines the functions toward which a
#'   strongly penalized fit contracts: a constant for `m = 1`, a straight
#'   line for 2, a parabola for 3, and the functions of [operator_null()]
#'   for any operator. Its order is at most `k - 1`, the highest degree the
#'   basis carries.
#' @param measure The measure against which the roughness is integrated.
#' @param constrain The directions to which the smooth is made orthogonal,
#'   `NULL` for the null space of the penalty.
#' @param null_space What happens to the directions that the penalty does
#'   not see.
#' @param reparam The coordinates in which the coefficients are expressed.
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
  order <- as_operator(order)
  # ABOVE THE HIGHEST DEGREE the Gram matrix is identically zero, so the
  # penalty would penalize nothing
  m_ord <- operator_order(order)
  if (m_ord > k - 1L) {
    stop(sprintf(paste0(
      "'order' (%d) exceeds the highest degree the basis carries (%d): the",
      "\n  derivative of that order is zero, so the penalty would be the",
      " zero matrix."
    ), m_ord, k - 1L), call. = FALSE)
  }
  null_space <- match.arg(null_space, c("keep", "drop", "shrink"))
  reparam <- match.arg(reparam, c("dr", "none", "orthonorm"))
  check_interval(lower, upper)
  measure <- check_measure(measure)
  penalty <- check_penalty(penalty)
  constrain <- check_constrain(constrain, order)
  ncon <- if (is.null(constrain)) m_ord else constrain + 1L
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
#' Building it this way instead of by folding a widened knot sequence makes
#' every quantity exact. The parent basis lives on \eqn{[a, b]},
#' the period itself, so its Gram matrix integrates over one period and the
#' penalty needs no numerical fallback: the roughness matrix is
#' \eqn{T^\top G T} with \eqn{G} the parent's own exact Gram. A folded
#' construction puts its parent on a widened interval, and its Gram then
#' integrates over more than one period.
#'
#' # Comparison with the Fourier basis
#'
#' Both families are periodic and differ in locality: a Fourier basis is
#' global and a periodic B-spline basis is local. For a smooth periodic
#' function, such as \eqn{\exp(\kappa \cos t)} at a low concentration
#' \eqn{\kappa}, whose Fourier coefficients are modified Bessel functions
#' and decay rapidly, the two give the same fit. A narrow seasonal feature,
#' at a higher concentration, is represented better by the local basis at
#' the same number of columns.
#'
#' At the two ends of the period the basis functions agree in value and in
#' their first `degree - 1` derivatives, up to rounding.
#'
#' # Arguments the family does not take
#'
#' `constrain` in the form "the polynomials up to degree c" has no meaning
#' here, for the same reason as on [fourier_smooth()]: a non-constant
#' periodic function is never a polynomial, so the null space of the
#' roughness matrix of a derivative penalty is the constant at every order,
#' instead of a space growing with the order.
#'
#' The constant is removed, so a model carrying an intercept spans the level
#' and the smooth carries the shape, and no free column is restored: a
#' linear column is not periodic. A block therefore carries `k - 1`
#' columns.
#'
#' # The interval is the period
#'
#' `lower` and `upper` are the ends of one cycle, and they are a property of
#' the problem rather than of the sample: day-of-year data run over
#' \eqn{[0, 365]} whether or not an observation falls on the first day of the
#' year. Left `NULL` they are read from the data, which makes the period the
#' observed range.
#'
#' @param k The number of periodic functions, a whole number of at least 3.
#'   The block carries `k - 1` columns, the constant being removed.
#' @param degree The degree of the underlying B-spline, `3` for a cubic. The
#'   fit matches at the two ends in its value and its first `degree - 1`
#'   derivatives.
#' @param order What the penalty measures: a [LinearOperator] from
#'   [deriv_operator()], [harmonic_operator()], [oscillator_operator()] or
#'   [linear_operator()], or a whole number `m` as the shorthand for
#'   `deriv_operator(m)`. The block is constrained against the constant
#'   only, so with a derivative penalty a strongly penalized fit contracts to
#'   a constant at every order. Its order is at most `degree`.
#' @param measure The measure against which the roughness is integrated.
#' @param null_space What happens to the directions that the penalty does
#'   not see. A periodic basis restores no free column, so the three settings
#'   give the same block, and `"shrink"` differs only in that it cannot be
#'   combined with a `penalty` factory.
#' @param reparam The coordinates in which the coefficients are expressed.
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
  order <- as_operator(order)
  m_ord <- operator_order(order)
  if (m_ord > degree) {
    stop(sprintf(paste0(
      "'order' (%d) exceeds 'degree' (%d): the derivative of that order of",
      " a\n  spline of that degree is zero, so the penalty would be the zero",
      " matrix."
    ), m_ord, degree), call. = FALSE)
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
      "\n  basis would carry %d functions instead of the %d requested."
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
#' Returns the matrix whose rows state that a function of `basis` and its
#' first `degree - 1` derivatives take the same value at the two ends of the
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
#' coefficients instead of an integrated squared derivative. The basis is
#' rich and the penalty cheap to compute: `k` is chosen large enough not to
#' limit the fit, and the smoothing parameter controls the roughness.
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
#' The argument `diff` belongs to this family alone. A difference penalty
#' reads the coefficients as an ordered sequence in which neighbors are
#' comparable. This holds for the coefficients of a B-spline and not for
#' those of a Fourier basis, where adjacent coefficients are a sine, a cosine
#' and the next sine, and their difference has no meaning.
#'
#' # The null space
#'
#' \eqn{D_d c = 0} exactly when the coefficients are a polynomial of degree
#' below \eqn{d} in their index, so the null space of the roughness matrix
#' has dimension `diff`.
#'
#' The differences are taken on the coefficients of the basis of Eilers and
#' Marx, whose knots are equally spaced with the same step beyond the
#' interval as inside it, and the penalty is carried onto the clamped basis
#' the package evaluates (the two span the same splines on the interval, so
#' only the coordinates change). On those knots the Greville abscissae are
#' equally spaced, and Marsden's identity makes the null space exactly the
#' polynomials of degree below `diff`. It is also the penalty of
#' \pkg{mgcv}'s `bs = "ps"` smooths, up to the normalization that
#' \pkg{mgcv} applies to its penalty matrices.
#'
#' # Comparison with the integrated penalty on the same basis
#'
#' The two are different penalties, and neither contains the other. The
#' difference operator carries no factor of the knot spacing, so the raw
#' P-spline roughness matrix is on a scale several orders of magnitude
#' smaller than the integrated-derivative matrix on the same basis. After the
#' Demmler-Reinsch rotation both penalties are the identity on the penalized
#' directions, so equal smoothing parameters give similar amounts of
#' smoothing.
#'
#' @param k The number of basis functions, a whole number greater than
#'   `degree` and greater than `diff`.
#' @param degree The degree of the B-spline, `3` for a cubic.
#' @param diff The order of difference the penalty takes, `2` for the usual
#'   construction. The fit contracts to a polynomial of degree `diff - 1`.
#' @param constrain The directions to which the smooth is made orthogonal,
#'   `NULL` for the null space of the penalty.
#' @param null_space What happens to the directions that the penalty does
#'   not see.
#' @param reparam The coordinates in which the coefficients are expressed.
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
      " %d\n  needs at least %d functions."
    ), k, degree, degree, degree + 1L), call. = FALSE)
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
  constrain <- check_constrain(constrain, deriv_operator(diff))
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
    order = deriv_operator(diff),
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
  # the coefficients rather than of the functions, so nothing is integrated.
  # The differences are taken on the coefficients of the Eilers-Marx basis,
  # whose knots are equally spaced beyond the interval as well as inside it,
  # and carried onto the clamped basis the package evaluates. The two bases
  # span the same splines, so only the coordinates of the penalty change.
  d <- base::diff(diag(1, b@dimension), differences = sm@diff) %*%
    pspline_map(b)
  out <- crossprod(d)
  nm <- basis_colnames(b)
  dimnames(out) <- list(nm, nm)
  out
}

#' The Eilers-Marx Coordinates of a Clamped B-Spline Basis
#'
#' @description
#' Returns the matrix \eqn{M} with \eqn{a = M c}, where \eqn{c} are the
#' coefficients of a function on the clamped B-spline basis `b` and \eqn{a}
#' are the coefficients of the same function on the basis of Eilers and Marx
#' (1996), whose knots are equally spaced with the same step beyond the
#' interval as inside it.
#'
#' @details
#' The two bases have the same degree and the same breakpoints inside the
#' interval, so on the interval they span the same splines and \eqn{M} is
#' square and invertible. It is computed by evaluating both bases on a grid
#' of ten points per basis function and solving the least-squares system,
#' which is exact to rounding. A difference penalty \eqn{\lVert D a\rVert^2}
#' then reads \eqn{c^\top M^\top D^\top D M c} on the clamped coefficients.
#' On the Eilers-Marx knots the Greville abscissae are equally spaced, so
#' coefficients that are a polynomial of degree below `diff` in their index
#' give exactly a polynomial of that degree.
#'
#' @param b A [BsplineBasis] built by [bspline_basis()], with equally spaced
#'   interior knots.
#'
#' @return A square numeric matrix of `b@dimension` rows and columns.
#'
#' @references
#' Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with B-splines
#' and penalties. *Statistical Science*, 11(2), 89-121.
#'
#' @examples
#' b <- bspline_basis(0, 1, dimension = 8)
#' M <- basis7:::pspline_map(b)
#' dim(M)
#' @keywords internal
pspline_map <- function(b) {
  k <- b@dimension
  m <- b@basis_params$degree
  h <- (b@upper - b@lower) / (k - m)
  knots <- b@lower + h * seq(-m, k)
  x <- seq(b@lower, b@upper, length.out = 10L * k + 1L)
  be <- splines::splineDesign(knots, x, ord = m + 1L, outer.ok = TRUE)
  bc <- as.matrix(basis_eval(b, x))
  qr.solve(be, bc)
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
#' smoothing parameter, this family carries `m` of them, and the data
#' determine how the roughness is distributed.
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
#' \eqn{S_i = D^\top \mathrm{diag}(v_i) D}, so [smoother_build()] returns
#' those `m` components and the layer that places the smooth combines them
#' into one penalty with `m` smoothing parameters. The weight profile is itself a
#' spline, and its coefficients are those smoothing parameters.
#'
#' The index basis is built over the range of the index, so an affine
#' relabeling of the index (\eqn{1, \ldots, n}, \eqn{i/n}, or the \eqn{i/k}
#' of \pkg{mgcv}) gives the same profile up to rounding.
#'
#' # Equal smoothing parameters
#'
#' The weight functions are a B-spline basis, hence a partition of unity, so
#' \eqn{\sum_i v_i = 1} and therefore \eqn{\sum_i S_i = D^\top D}, up to
#' rounding. Holding every \eqn{\lambda_i} at one value gives the P-spline
#' penalty at that value, so [pspline_smooth()] is a special case of this
#' family, and the additional freedom is used only where the data require
#' it.
#'
#' # Comparison with a single smoothing parameter
#'
#' On a function whose roughness varies along the covariate the adaptive
#' penalty smooths the quiet stretches more strongly than a P-spline with
#' one smoothing parameter, and follows the noise less there. On a function
#' of constant roughness its additional smoothing parameters have nothing to
#' adapt to, and the fit is slightly worse than the P-spline fit.
#'
#' The differences are taken in the Eilers-Marx coordinates of
#' [pspline_smooth()], as in the adaptive P-spline of \pkg{mgcv}'s
#' `bs = "ad"`. At equal smoothing parameters the two sums are the same
#' P-spline penalty; the weight functions differ, \pkg{mgcv} building them
#' from a different basis over the index, so the two fits are close and not
#' identical.
#'
#' # The coordinates
#'
#' `reparam = "dr"` is rejected: the Demmler-Reinsch rotation diagonalizes
#' the pencil of the Gram matrix against a single penalty, this family has
#' `m` penalties, and a rotation that makes one component the identity
#' leaves the others unspecified. The default is `"none"`, and `"orthonorm"`
#' is also available.
#'
#' # The rank of the penalty
#'
#' Each \eqn{S_i} has a large null space, made of the directions on which its
#' weight function vanishes. The null space of the sum is the intersection
#' of the null spaces of the components and does not depend on the
#' smoothing parameters, so the rank is read from the components and not
#' from the assembled \eqn{S(\lambda)}, whose numerical rank falls as one
#' smoothing parameter grows and depends on the tolerance used.
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
#' @param constrain The directions to which the smooth is made orthogonal,
#'   `NULL` for the null space of the penalty.
#' @param null_space What happens to the directions that the penalty does
#'   not see. `"shrink"` is rejected, because the shrinkage is a fraction of
#'   the eigenvalues of one roughness matrix and this penalty is a sum of
#'   components.
#' @param reparam The coordinates in which the coefficients are expressed,
#'   `"none"` or `"orthonorm"`. `"dr"` is rejected.
#' @param penalty Any value other than `NULL` signals an error. The penalty of
#'   this family is the sum of its components, and a factory would replace it
#'   with one penalty over the coefficients, which is [pspline_smooth()] with
#'   that factory.
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
                   "a single component is pspline_smooth(), on the same basis")
  if (k < degree + 1L) {
    stop(sprintf(paste0(
      "'k' (%d) is too small for 'degree' (%d): a B-spline basis of degree",
      " %d\n  needs at least %d functions."
    ), k, degree, degree, degree + 1L), call. = FALSE)
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
  constrain <- check_constrain(constrain, deriv_operator(diff))
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
    order = deriv_operator(diff),
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
  # the differences are taken in the Eilers-Marx coordinates, as for
  # pspline_smooth(), so that equal smoothing parameters give its penalty
  d <- base::diff(diag(1, b@dimension), differences = sm@diff) %*%
    pspline_map(b)
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
