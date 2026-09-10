#' The Linear Differential Operator Class
#' @name LinearOperator
#'
#' @description
#' A linear differential operator with constant coefficients,
#' \deqn{L x = w_0 x + w_1 Dx + \cdots + w_{m-1} D^{m-1} x + D^m x,}
#' which a smoother uses to say what its penalty measures. The leading
#' coefficient is 1 by construction, so an operator is determined by the
#' weights \eqn{w_0, \ldots, w_{m-1}} and by nothing else.
#'
#' @details
#' Build one with [deriv_operator()], [harmonic_operator()] or
#' [linear_operator()], and compose two with `*`. The properties are not
#' meant to be read directly: [operator_order()], [operator_weights()] and
#' [operator_null()] are the accessors, and they answer for an operator whose
#' period is still to be resolved where reading the property would not.
#'
#' @param weights The coefficients \eqn{w_0, \ldots, w_{m-1}}, a numeric
#'   vector, or `NA` of the right length for an operator whose period has not
#'   been resolved yet.
#' @param operator_name A short name used when the operator prints.
#' @param operator_params A list of whatever the constructor recorded, such
#'   as the period and the number of harmonics.
#'
#' @return An S7 object of class `LinearOperator`.
#'
#' @seealso [deriv_operator()], [harmonic_operator()] and
#'   [linear_operator()], which build one; [operator_null()] for the space a
#'   maximal penalty leaves untouched.
#'
#' @examples
#' op <- harmonic_operator(365)
#' c(class = class(op)[1], order = operator_order(op))
#' @export
LinearOperator <- S7::new_class(
  "LinearOperator",
  properties = list(
    weights = S7::class_numeric,
    operator_name = S7::class_character,
    operator_params = S7::class_list
  ),
  validator = function(self) {
    if (length(self@operator_name) != 1L || is.na(self@operator_name)) {
      return("@operator_name must be a single string")
    }
    if (length(self@weights) < 1L) {
      return("@weights must have at least one entry")
    }
    if (!all(is.na(self@weights)) && !all(is.finite(self@weights))) {
      return("@weights must be finite, or NA throughout")
    }
    NULL
  }
)


#' The Derivative Operator
#'
#' @description
#' The operator \eqn{L x = D^m x}, whose penalty
#' \eqn{\int (D^m x)^2} is the integrated squared derivative every smoother
#' in this package penalized with before operators existed. Writing
#' `order = m` on a smoother is the shorthand for `order = deriv_operator(m)`
#' and builds the identical construction.
#'
#' @details
#' Its null space is the polynomials of degree below \eqn{m}: the
#' characteristic polynomial is \eqn{r^m}, a root at zero of multiplicity
#' \eqn{m}, which contributes \eqn{1, t, \ldots, t^{m-1}}. A strongly
#' penalized fit therefore contracts to a constant at `m = 1`, to a straight
#' line at `m = 2` and to a parabola at `m = 3`.
#'
#' @param m The order of the derivative, a whole number of at least 1.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @seealso [harmonic_operator()] for the periodic one, [linear_operator()]
#'   for the general form, and [operator_null()] for what a maximal penalty
#'   leaves.
#'
#' @examples
#' deriv_operator(2)
#'
#' # what a strongly penalized fit contracts to
#' operator_null(deriv_operator(3))
#'
#' # `order = 2` on a smoother is this operator, and gives the same block
#' set.seed(1)
#' x <- sort(runif(100))
#' a <- smoother_build(bspline_smooth(k = 8, order = 2), x)
#' b <- smoother_build(bspline_smooth(k = 8, order = deriv_operator(2)), x)
#' identical(a$X, b$X) && identical(a$S, b$S)
#' @export
deriv_operator <- function(m = 2) {
  m <- check_whole(m, "m", 1L)
  LinearOperator(
    weights = rep(0, m),
    operator_name = "derivative",
    operator_params = list(kind = "deriv", m = m)
  )
}


#' The Harmonic Acceleration Operator
#'
#' @description
#' The operator for periodic data,
#' \deqn{L x = D \prod_{i=1}^{h} \left(D^2 + (i\nu)^2\right) x,
#'   \qquad \nu = \frac{2\pi}{T},}
#' whose null space is the constant together with the first \eqn{h}
#' harmonics of the period \eqn{T}. At `harmonics = 1` it is the harmonic
#' acceleration operator \eqn{L x = \nu^2 Dx + D^3 x} of Ramsay and
#' Silverman, which is what `fda::vec2Lfd(c(0, (2*pi/T)^2, 0))` builds.
#'
#' @details
#' # Why a periodic basis wants it
#'
#' The derivative operator asks a fit to contract toward a straight line,
#' and a straight line is not what a cyclic phenomenon simplifies to: it is
#' not even periodic. What a seasonal series contracts to is its fundamental
#' harmonic, a constant level plus one sine and one cosine of the period, and
#' this is the operator that leaves exactly that alone.
#'
#' The difference is visible in the fit rather than only in the algebra. On
#' 300 observations of a truth \eqn{2 + 1.5\sin\nu t + \cos\nu t +
#' 0.35\sin 3\nu t}, at a smoothing parameter large enough to flatten the
#' higher harmonics, the fit under \eqn{D^2} has fallen to a standard
#' deviation of 0.745 and a fundamental amplitude of 1.053 against a true
#' 1.803, while the fit under this operator keeps 1.271 and 1.795 and is a
#' pure sinusoid to 3.9e-08.
#'
#' # The period
#'
#' `period` is \eqn{T}, the length of one cycle, in the units of the
#' covariate: 365 for a day of the year, 12 for a month, 24 for an hour.
#' It is **not** \eqn{\nu = 2\pi/T}, which some of the literature also calls
#' omega. `NULL` leaves it to be resolved from the interval of whichever
#' smoother the operator is given to, which is what
#' [fourier_smooth()] does by default; [operator_resolve()] is what fills it
#' in, and the accessors report an unresolved operator as such rather than
#' guessing.
#'
#' @param period The length of one cycle, a positive number, or `NULL` to
#'   take the width of the smoother's interval.
#' @param harmonics How many harmonics the penalty leaves unpenalized, a
#'   whole number of at least 1. The order of the operator is
#'   `2 * harmonics + 1` and the null space has that dimension.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @seealso [fourier_smooth()], whose default it is; [deriv_operator()] for
#'   the non-periodic one; [operator_null()] for the space it leaves.
#'
#' @references
#' Ramsay, J. O. and Silverman, B. W. (2005). *Functional Data Analysis*,
#' second edition. Springer, chapter 5.
#'
#' @examples
#' harmonic_operator(365)
#'
#' # the null space is the constant and the fundamental
#' operator_null(harmonic_operator(365))
#'
#' # two harmonics left free instead of one
#' harmonic_operator(365, harmonics = 2)
#'
#' # the weights are fda's: c(w0, w1, w2) = c(0, (2*pi/365)^2, 0)
#' round(operator_weights(harmonic_operator(365)), 8)
#'
#' # the period may be left to the smoother's interval
#' operator_weights(harmonic_operator())
#' @export
harmonic_operator <- function(period = NULL, harmonics = 1) {
  harmonics <- check_whole(harmonics, "harmonics", 1L)
  check_period(period)
  op <- LinearOperator(
    weights = if (is.null(period)) {
      rep(NA_real_, 2L * harmonics + 1L)
    } else {
      harmonic_weights(period, harmonics)
    },
    operator_name = "harmonic acceleration",
    operator_params = list(
      kind = "harmonic", period = period, harmonics = harmonics
    )
  )
  op
}


#' The Oscillator Operator
#'
#' @description
#' The periodic factor on its own,
#' \deqn{L x = \prod_{i=1}^{h}\left(D^2 + (i\nu)^2\right) x,
#'   \qquad \nu = \frac{2\pi}{T},}
#' whose null space is the first \eqn{h} harmonics of the period \eqn{T}
#' **without** the constant. At `harmonics = 1` it is the equation of simple
#' harmonic motion, \eqn{L x = \nu^2 x + D^2 x}.
#'
#' @details
#' It exists to be composed. [harmonic_operator()] is this operator with a
#' leading \eqn{D}, which is what puts the constant into the null space, so
#' the two are related by
#' `harmonic_operator(T, h)` being `deriv_operator(1) * oscillator_operator(T, h)`
#' exactly, and a penalty that should leave a **linear trend** and a cycle
#' alone rather than a level and a cycle is
#' `deriv_operator(2) * oscillator_operator(T)`.
#'
#' Composing with [harmonic_operator()] instead would raise the order by one
#' and put an extra power of \eqn{t} in the null space, which is why the
#' factor is offered separately rather than left to be written out by hand.
#'
#' @param period The length of one cycle, a positive number, or `NULL` to
#'   take the width of the smoother's interval.
#' @param harmonics How many harmonics the penalty leaves unpenalized, a
#'   whole number of at least 1. The order of the operator is
#'   `2 * harmonics`.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @seealso [harmonic_operator()], which is this with a leading derivative;
#'   [operator-compose] for what composing does to the null space.
#'
#' @examples
#' oscillator_operator(365)
#'
#' # the null space is the fundamental alone: no constant
#' operator_null(oscillator_operator(365))
#'
#' # a linear trend and a yearly cycle, and nothing else
#' op <- deriv_operator(2) * oscillator_operator(365)
#' operator_null(op)
#'
#' # the harmonic operator is this one with a leading derivative
#' a <- operator_weights(deriv_operator(1) * oscillator_operator(365))
#' b <- operator_weights(harmonic_operator(365))
#' all.equal(a, b)
#' @export
oscillator_operator <- function(period = NULL, harmonics = 1) {
  harmonics <- check_whole(harmonics, "harmonics", 1L)
  check_period(period)
  LinearOperator(
    weights = if (is.null(period)) {
      rep(NA_real_, 2L * harmonics)
    } else {
      oscillator_weights(period, harmonics)
    },
    operator_name = "oscillator",
    operator_params = list(
      kind = "oscillator", period = period, harmonics = harmonics
    )
  )
}


#' Refuse a Period That Is Not One Positive Number
#'
#' @param period The period to check, or `NULL`.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @keywords internal
check_period <- function(period) {
  if (is.null(period)) return(invisible(NULL))
  if (!is.numeric(period) || length(period) != 1L || !is.finite(period) ||
    period <= 0) {
    stop(paste0(
      "'period' must be a single positive number, or NULL to take the",
      " width of\n  the smoother's interval. It is the length of one",
      " cycle in the units of the\n  covariate (365 for a day of the",
      " year), not 2*pi over that length."
    ), call. = FALSE)
  }
  invisible(NULL)
}


#' The Weights of an Oscillator Operator
#'
#' @description
#' Multiplies out \eqn{\prod_{i \le h} (r^2 + (i\nu)^2)} and drops the
#' leading coefficient, which is 1.
#'
#' @param period The length of one cycle.
#' @param harmonics How many harmonics to leave in the null space.
#'
#' @return A numeric vector of `2 * harmonics` entries.
#'
#' @keywords internal
oscillator_weights <- function(period, harmonics) {
  nu <- 2 * pi / period
  p <- 1
  for (i in seq_len(harmonics)) p <- poly_mul(p, c((i * nu)^2, 0, 1))
  utils::head(p, -1L)
}


#' The Weights of a Harmonic Acceleration Operator
#'
#' @description
#' Multiplies out \eqn{r \prod_{i \le h} (r^2 + (i\nu)^2)} and drops the
#' leading coefficient, which is 1.
#'
#' @param period The length of one cycle.
#' @param harmonics How many harmonics to leave in the null space.
#'
#' @return A numeric vector of `2 * harmonics` entries, the weights
#'   \eqn{w_0, \ldots, w_{m-1}}.
#'
#' @keywords internal
harmonic_weights <- function(period, harmonics) {
  nu <- 2 * pi / period
  # the characteristic polynomial in increasing degree, starting from r
  p <- c(0, 1)
  for (i in seq_len(harmonics)) p <- poly_mul(p, c((i * nu)^2, 0, 1))
  utils::head(p, -1L)
}


#' Multiply Two Polynomials Given by Their Coefficients
#'
#' @description
#' Convolves two coefficient vectors written in increasing degree, which is
#' what composing two differential operators does to their characteristic
#' polynomials.
#'
#' @param a,b Numeric vectors of coefficients in increasing degree.
#'
#' @return A numeric vector of `length(a) + length(b) - 1` coefficients, in
#'   increasing degree.
#'
#' @keywords internal
poly_mul <- function(a, b) {
  out <- numeric(length(a) + length(b) - 1L)
  for (i in seq_along(a)) {
    j <- seq.int(i, i + length(b) - 1L)
    out[j] <- out[j] + a[[i]] * b
  }
  out
}


#' A Linear Differential Operator From Its Weights
#'
#' @description
#' The general form: `linear_operator(w)` is
#' \eqn{L x = w_1 x + w_2 Dx + \cdots + w_m D^{m-1} x + D^m x}, with the
#' leading coefficient 1 and `w` given in increasing order of
#' differentiation. It is the vector `fda::vec2Lfd()` takes.
#'
#' @details
#' The null space is read from the roots of the characteristic polynomial
#' \eqn{r^m + w_{m-1} r^{m-1} + \cdots + w_0}: a real root \eqn{a} of
#' multiplicity \eqn{\mu} contributes \eqn{t^i e^{at}} for \eqn{i < \mu}, and
#' a complex pair \eqn{a \pm bi} contributes \eqn{t^i e^{at}\cos(bt)} and
#' \eqn{t^i e^{at}\sin(bt)}. [operator_null()] reports it.
#'
#' Use [deriv_operator()] and [harmonic_operator()] where they apply: they
#' say what the operator is for, and they carry a period that can be
#' resolved from a smoother's interval, which a bare weight vector cannot.
#'
#' @param w The weights \eqn{w_0, \ldots, w_{m-1}}, a numeric vector of at
#'   least one finite entry. Its length is the order of the operator.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @seealso [deriv_operator()] and [harmonic_operator()] for the named
#'   instances, [operator_null()] for the null space.
#'
#' @examples
#' # the harmonic acceleration operator written out, as in fda
#' nu2 <- (2 * pi / 365)^2
#' linear_operator(c(0, nu2, 0))
#'
#' # and it is the same operator harmonic_operator() builds
#' all.equal(operator_weights(linear_operator(c(0, nu2, 0))),
#'           operator_weights(harmonic_operator(365)))
#'
#' # exponential growth left unpenalized: L x = Dx - r x
#' operator_null(linear_operator(-0.5))
#' @export
linear_operator <- function(w) {
  if (!is.numeric(w) || length(w) < 1L || !all(is.finite(w))) {
    stop(paste0(
      "'w' must be a numeric vector of at least one finite entry, the",
      " weights\n  w[0], ..., w[m-1] in increasing order of",
      " differentiation. The leading\n  coefficient is 1 and is not given."
    ), call. = FALSE)
  }
  LinearOperator(
    weights = as.numeric(w),
    operator_name = "linear",
    operator_params = list(kind = "linear")
  )
}


#' Compose Two Differential Operators
#' @name operator-compose
#'
#' @description
#' `L1 * L2` is the operator that applies one after the other. Composition
#' multiplies the characteristic polynomials, so the order of the product is
#' the sum of the orders and its null space is the union of the two null
#' spaces.
#'
#' @details
#' The union is what makes composition worth having: an operator whose
#' penalty leaves both a linear trend and a seasonal cycle alone is
#' `deriv_operator(2) * oscillator_operator(365)`, and there is no reason to
#' write its four weights by hand. Composition commutes, these operators
#' having constant coefficients.
#'
#' ⚠️ Compose with [oscillator_operator()] and not with
#' [harmonic_operator()], unless the extra power of \eqn{t} is wanted. The
#' harmonic operator already carries a leading \eqn{D}, so
#' `deriv_operator(2) * harmonic_operator(365)` is \eqn{D^3(D^2 + \nu^2)},
#' of order 5, whose null space is \eqn{1, t, t^2, \sin\nu t, \cos\nu t}. The
#' arithmetic is right either way and the dimension of the null space is
#' always the order; which of the two is meant is the thing to read off
#' [operator_null()] before fitting.
#'
#' An operator whose period has not been resolved cannot be composed: the
#' weights of the product depend on the period, and resolving afterwards
#' would have to know which factor it came from. Give the period, or compose
#' after [operator_resolve()].
#'
#' @param e1,e2 Two [LinearOperator] objects.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @seealso [operator_null()], which reports the union.
#'
#' @examples
#' # a linear trend and a yearly cycle, both left unpenalized
#' op <- deriv_operator(2) * oscillator_operator(365)
#' op
#' operator_null(op)
#'
#' # the order is the sum, the null space the union
#' c(operator_order(deriv_operator(2)), operator_order(oscillator_operator(365)),
#'   operator_order(op))
#'
#' # composing with the harmonic operator instead raises the order by one
#' # and adds a power of t: read the null space, not the name
#' operator_null(deriv_operator(2) * harmonic_operator(365))$label
#' @keywords internal
S7::method(`*`, list(LinearOperator, LinearOperator)) <- function(e1, e2) {
  for (op in list(e1, e2)) {
    if (!operator_resolved(op)) {
      stop(paste0(
        "an operator whose period is NULL cannot be composed: the weights",
        " of the\n  product depend on it. Give 'period', or compose after",
        " operator_resolve()."
      ), call. = FALSE)
    }
  }
  p <- poly_mul(c(e1@weights, 1), c(e2@weights, 1))
  LinearOperator(
    weights = utils::head(p, -1L),
    operator_name = paste(e1@operator_name, "x", e2@operator_name),
    operator_params = list(kind = "product")
  )
}


#' The Order of a Differential Operator
#'
#' @description
#' The highest derivative the operator takes, which is the dimension of its
#' null space and the length of its weight vector.
#'
#' @param op A [LinearOperator].
#'
#' @return A single integer.
#'
#' @seealso [operator_weights()], [operator_null()].
#'
#' @examples
#' c(deriv = operator_order(deriv_operator(2)),
#'   harmonic = operator_order(harmonic_operator(365)),
#'   two_harmonics = operator_order(harmonic_operator(365, harmonics = 2)))
#' @export
operator_order <- function(op) {
  check_operator(op)
  length(op@weights)
}


#' The Weights of a Differential Operator
#'
#' @description
#' The coefficients \eqn{w_0, \ldots, w_{m-1}} of
#' \eqn{L x = \sum_j w_j D^j x + D^m x}, in increasing order of
#' differentiation. The leading coefficient is 1 and is not among them.
#'
#' @param op A [LinearOperator].
#'
#' @return A numeric vector of `operator_order(op)` entries, or `NA` of that
#'   length where the operator's period has not been resolved.
#'
#' @seealso [operator_resolve()], which fills in a period.
#'
#' @examples
#' round(operator_weights(harmonic_operator(365)), 8)
#' operator_weights(deriv_operator(3))
#' @export
operator_weights <- function(op) {
  check_operator(op)
  op@weights
}


#' Has This Operator's Period Been Resolved?
#'
#' @description
#' Reports whether the weights are numbers rather than the `NA` placeholder
#' [harmonic_operator()] leaves when it is given no period.
#'
#' @param op A [LinearOperator].
#'
#' @return A single logical.
#'
#' @keywords internal
operator_resolved <- function(op) !anyNA(op@weights)


#' Fill In an Operator's Period From an Interval
#'
#' @description
#' Returns the operator with its period set to the width of `[lower, upper]`
#' where it was left `NULL`, and unchanged where it was given. It is what a
#' smoother calls before building its penalty, so that
#' `fourier_smooth(lower = 0, upper = 365)` penalizes on a cycle of 365
#' without the period having to be written twice.
#'
#' @param op A [LinearOperator].
#' @param lower,upper The interval, two numbers.
#'
#' @return An S7 object of class [LinearOperator], with numeric weights.
#'
#' @seealso [harmonic_operator()], whose `period = NULL` this resolves.
#'
#' @examples
#' operator_weights(harmonic_operator())
#' round(operator_weights(operator_resolve(harmonic_operator(), 0, 365)), 8)
#'
#' # an operator given a period keeps it, whatever the interval is
#' round(operator_weights(operator_resolve(harmonic_operator(12), 0, 365)), 8)
#' @export
operator_resolve <- function(op, lower, upper) {
  check_operator(op)
  if (operator_resolved(op)) return(op)
  p <- op@operator_params
  if (!p$kind %in% c("harmonic", "oscillator")) {
    stop("only a periodic operator carries a period to resolve.",
      call. = FALSE
    )
  }
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != 1L ||
    length(upper) != 1L || !is.finite(lower) || !is.finite(upper) ||
    upper <= lower) {
    stop("'lower' and 'upper' must be two numbers with 'upper' the larger.",
      call. = FALSE
    )
  }
  f <- if (identical(p$kind, "harmonic")) {
    harmonic_operator
  } else {
    oscillator_operator
  }
  f(period = upper - lower, harmonics = p$harmonics)
}


#' The Null Space of a Differential Operator
#'
#' @description
#' The functions a maximal penalty leaves untouched, that is the solutions of
#' \eqn{L x = 0}. They are what a strongly penalized fit contracts to, so
#' this is the one thing to look at before choosing an operator.
#'
#' @details
#' The null space is read from the roots of the characteristic polynomial
#' \eqn{r^m + w_{m-1}r^{m-1} + \cdots + w_0}, computed by [base::polyroot()] and
#' grouped by multiplicity. A real root \eqn{a} of multiplicity \eqn{\mu}
#' contributes \eqn{t^i e^{at}} for \eqn{i < \mu}, and a complex pair
#' \eqn{a \pm bi} of multiplicity \eqn{\mu} contributes both
#' \eqn{t^i e^{at}\cos(bt)} and \eqn{t^i e^{at}\sin(bt)}.
#'
#' The null space is a property of the **operator** and is computed from it,
#' never from the rank of an assembled penalty matrix. The two are different
#' questions and they give different answers: the same two-harmonic operator
#' has a penalty of null dimension 5 on a Fourier basis and 3 on a cubic
#' B-spline at a relative tolerance of 1e-10, because a spline represents a
#' sine only approximately, while the operator's own null space is
#' five-dimensional in both cases. Which of those functions a **basis** can
#' carry is a separate question, and [smoother_build()] is where it is asked.
#'
#' @param op A [LinearOperator].
#' @param tol The relative tolerance at which two roots count as one, a
#'   positive number.
#'
#' @return A data frame of one row per function in the null space, with
#'   columns `label` (how the function reads), `rate` (the real part of the
#'   root), `freq` (the imaginary part) and `degree` (the power of \eqn{t}
#'   multiplying it).
#'
#' @seealso [operator_null_design()] for those functions evaluated,
#'   [deriv_operator()] and [harmonic_operator()] for the two named
#'   instances.
#'
#' @examples
#' # a straight line
#' operator_null(deriv_operator(2))
#'
#' # a constant plus the fundamental of a yearly cycle
#' operator_null(harmonic_operator(365))
#'
#' # both at once
#' operator_null(deriv_operator(2) * harmonic_operator(365))
#'
#' # the dimension is always the order
#' nrow(operator_null(harmonic_operator(365, harmonics = 3)))
#' @export
operator_null <- function(op, tol = 1e-6) {
  check_operator(op)
  if (!operator_resolved(op)) {
    stop(paste0(
      "the operator's period is NULL, so its null space is not determined",
      " yet.\n  Give 'period', or call operator_resolve() with the",
      " interval."
    ), call. = FALSE)
  }
  w <- op@weights
  m <- length(w)
  if (identical(op@operator_params$kind, "deriv")) {
    # EXACT, and not a root-finding problem: the characteristic polynomial
    # is r^m, so the null space is the polynomials of degree below m.
    return(null_frame(rate = rep(0, m), freq = rep(0, m),
                      degree = seq_len(m) - 1L))
  }
  rt <- polyroot(c(w, 1))
  cl <- cluster_roots(rt, tol)
  rate <- numeric(0)
  freq <- numeric(0)
  deg <- integer(0)
  part <- character(0)
  for (g in cl) {
    a <- Re(g$root)
    b <- Im(g$root)
    # a conjugate pair is taken once, at the positive side, and emits both
    # of its real functions there
    if (b < -tol) next
    real_root <- abs(b) <= tol
    for (i in seq_len(g$mult) - 1L) {
      if (real_root) {
        rate <- c(rate, a)
        freq <- c(freq, 0)
        deg <- c(deg, i)
        part <- c(part, "")
      } else {
        rate <- c(rate, a, a)
        freq <- c(freq, b, b)
        deg <- c(deg, i, i)
        part <- c(part, "sin", "cos")
      }
    }
  }
  null_frame(rate, freq, deg, part)
}


#' Group Complex Roots by Multiplicity
#'
#' @description
#' Clusters the output of [base::polyroot()] so that roots within `tol` of one
#' another, relative to their own size, are read as one root of higher
#' multiplicity.
#'
#' @param rt The complex roots.
#' @param tol The relative tolerance.
#'
#' @return A list of `root` and `mult` pairs, ordered by the real part and
#'   then the imaginary part.
#'
#' @keywords internal
cluster_roots <- function(rt, tol) {
  out <- list()
  for (r in rt) {
    hit <- FALSE
    for (i in seq_along(out)) {
      if (Mod(r - out[[i]]$root) <= tol * max(1, Mod(r))) {
        out[[i]]$mult <- out[[i]]$mult + 1L
        hit <- TRUE
        break
      }
    }
    if (!hit) out[[length(out) + 1L]] <- list(root = r, mult = 1L)
  }
  ord <- order(vapply(out, function(e) Re(e$root), numeric(1)),
    vapply(out, function(e) Im(e$root), numeric(1))
  )
  out[ord]
}


#' Assemble the Null-Space Table
#'
#' @description
#' Builds the data frame [operator_null()] returns, writing the label of each
#' function from its rate, frequency and power of \eqn{t}.
#'
#' @param rate,freq,degree The three descriptions, of equal length.
#' @param part Which half of a conjugate pair the row is, `"sin"` or
#'   `"cos"`, and `""` for a row that comes from a real root. It is carried
#'   as a column rather than deduced from the row's position: a real root
#'   ahead of a pair shifts every position after it, and a pair read off the
#'   parity would then come out as two sines.
#'
#' @return A data frame of five columns.
#'
#' @keywords internal
null_frame <- function(rate, freq, degree, part = rep("", length(rate))) {
  lab <- vapply(seq_along(rate), function(i) {
    a <- rate[[i]]
    b <- freq[[i]]
    d <- degree[[i]]
    parts <- character(0)
    if (d == 1L) parts <- c(parts, "t")
    if (d > 1L) parts <- c(parts, sprintf("t^%d", d))
    if (abs(a) > 1e-12) parts <- c(parts, sprintf("exp(%g t)", a))
    if (nzchar(part[[i]])) {
      parts <- c(parts, sprintf("%s(%g t)", part[[i]], b))
    }
    if (length(parts) == 0L) "1" else paste(parts, collapse = " ")
  }, character(1))
  data.frame(
    label = lab, rate = rate, freq = freq, degree = as.integer(degree),
    part = part, stringsAsFactors = FALSE
  )
}


#' The Null-Space Functions Evaluated
#'
#' @description
#' The design matrix of the functions [operator_null()] lists, evaluated at
#' `x`: one column per function, in the order that table gives them.
#'
#' @details
#' These are the mathematical functions, unscaled. A smoother that restores
#' them as free columns centers and scales them first and records what it
#' did, so that a prediction reapplies the same columns rather than
#' recomputing them from new data; [smoother_span()] is where that happens.
#'
#' A rate far from zero over a wide interval overflows: \eqn{e^{at}} is what
#' it is, and an operator with a large real root is a statement about growth
#' that a design matrix cannot hold. The two operators this package builds
#' for itself, [deriv_operator()] and [harmonic_operator()], have purely
#' imaginary roots and no such difficulty.
#'
#' @param op A [LinearOperator].
#' @param x The points to evaluate at, a numeric vector.
#'
#' @return A numeric matrix of `length(x)` rows and `operator_order(op)`
#'   columns, with the labels of [operator_null()] as column names.
#'
#' @seealso [operator_null()] for what the columns are.
#'
#' @examples
#' x <- seq(0, 365, length.out = 5)
#' round(operator_null_design(harmonic_operator(365), x), 4)
#'
#' # L applied to its own null space is zero, which is what makes it the
#' # null space; here to the accuracy of a central difference
#' round(colSums(abs(operator_null_design(deriv_operator(2), 1:5))), 6)
#' @export
operator_null_design <- function(op, x) {
  nl <- operator_null(op)
  if (!is.numeric(x)) stop("'x' must be numeric.", call. = FALSE)
  out <- vapply(seq_len(nrow(nl)), function(i) {
    v <- exp(nl$rate[[i]] * x) * x^nl$degree[[i]]
    if (nzchar(nl$part[[i]])) {
      v <- v * if (identical(nl$part[[i]], "cos")) {
        cos(nl$freq[[i]] * x)
      } else {
        sin(nl$freq[[i]] * x)
      }
    }
    v
  }, numeric(length(x)))
  out <- matrix(out, length(x), nrow(nl), dimnames = list(NULL, nl$label))
  out
}


#' Apply an Operator to a Basis
#'
#' @description
#' Evaluates \eqn{Lb(x)}, the operator applied to every function of a basis:
#' the matrix whose Gram matrix is the roughness the penalty measures.
#'
#' @details
#' It is one weighted sum of derivative evaluations,
#' \eqn{\sum_j w_j D^j b + D^m b}, and every term of it comes from
#' [basis_deriv()], so a basis that answers its derivatives answers this.
#' Terms whose weight is exactly zero are skipped, which is why the pure
#' derivative operator costs one call and not \eqn{m + 1}.
#'
#' @param basis A [basis].
#' @param x The points to evaluate at.
#' @param op A [LinearOperator], with its period resolved.
#'
#' @return A numeric matrix of `length(x)` rows and `basis@dimension`
#'   columns.
#'
#' @keywords internal
operator_eval <- function(basis, x, op) {
  w <- operator_weights(op)
  out <- basis_deriv(basis, x, order = length(w))
  for (j in seq_along(w)) {
    if (w[[j]] != 0) out <- out + w[[j]] * basis_deriv(basis, x, order = j - 1L)
  }
  out
}


#' Is This a Differential Operator?
#'
#' @description
#' Reports whether an object is a [LinearOperator], which is how the `order`
#' argument of a smoother and of [basis_gram()] tells the two spellings
#' apart.
#'
#' @param x Any object.
#'
#' @return A single logical.
#'
#' @keywords internal
is_operator <- function(x) S7::S7_inherits(x, LinearOperator)


#' Refuse Anything That Is Not an Operator
#'
#' @param op The object to check.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @keywords internal
check_operator <- function(op) {
  if (!is_operator(op)) {
    stop(paste0(
      "expected a linear differential operator, from deriv_operator(),",
      " harmonic_operator()\n  or linear_operator()."
    ), call. = FALSE)
  }
  invisible(NULL)
}


#' Read an Order Argument as an Operator
#'
#' @description
#' Normalizes the `order` argument every smoother family takes: a whole
#' number `m` becomes `deriv_operator(m)` and an operator is returned
#' unchanged. It is what makes `order = 2` the shorthand rather than a
#' second way of saying the same thing.
#'
#' @param order A whole number or a [LinearOperator].
#' @param nm The argument's name, for the error message.
#'
#' @return An S7 object of class [LinearOperator].
#'
#' @keywords internal
as_operator <- function(order, nm = "order") {
  if (is_operator(order)) return(order)
  if (is.numeric(order)) return(deriv_operator(check_whole(order, nm, 1L)))
  stop(sprintf(paste0(
    "'%s' must be a whole number of at least 1, or a linear differential",
    " operator\n  from deriv_operator(), harmonic_operator() or",
    " linear_operator()."
  ), nm), call. = FALSE)
}


#' Is This the Plain Derivative Operator?
#'
#' @description
#' Reports whether an operator is \eqn{D^m}, which is the case every family
#' had before operators existed and the one whose construction must not move.
#'
#' @param op A [LinearOperator].
#'
#' @return A single logical.
#'
#' @keywords internal
is_deriv_operator <- function(op) {
  identical(op@operator_params$kind, "deriv")
}


#' How a Differential Operator Prints
#' @name print.LinearOperator
#'
#' @description
#' Three lines: the order, the operator written out, and the functions its
#' null space holds.
#'
#' @param x A [LinearOperator].
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' harmonic_operator(365)
#' deriv_operator(2) * harmonic_operator(365)
#' @keywords internal
S7::method(print, LinearOperator) <- function(x, ...) {
  cat(sprintf("<linear differential operator>  order %d\n",
              operator_order(x)))
  cat(sprintf("  L x = %s\n", operator_formula(x)))
  if (operator_resolved(x)) {
    cat(sprintf("  null space: %s\n",
                paste(operator_null(x)$label, collapse = ", ")))
  } else {
    cat(sprintf("  period: the smoother's interval (%d harmonic%s)\n",
      x@operator_params$harmonics,
      if (x@operator_params$harmonics == 1L) "" else "s"
    ))
  }
  invisible(x)
}


#' The Operator Written Out
#'
#' @description
#' Assembles the string `print()` shows, dropping the terms whose weight is
#' zero and writing the leading term without a coefficient.
#'
#' @param op A [LinearOperator].
#'
#' @return A single string.
#'
#' @keywords internal
operator_formula <- function(op) {
  w <- op@weights
  m <- length(w)
  if (!operator_resolved(op)) {
    return(sprintf("(%d weights, to be resolved) + D^%d x", m, m))
  }
  d <- function(j) {
    if (j == 0L) "x" else if (j == 1L) "Dx" else sprintf("D^%d x", j)
  }
  parts <- character(0)
  for (j in seq_along(w)) {
    if (w[[j]] != 0) {
      parts <- c(parts, sprintf("%s %s", format(w[[j]], digits = 4),
                                d(j - 1L)))
    }
  }
  paste(c(parts, d(m)), collapse = " + ")
}
