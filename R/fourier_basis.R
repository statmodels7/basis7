#' @include numerical_fallbacks.R
NULL


#' Fourier Basis
#'
#' @description
#' The S7 class of Fourier bases, the objects [fourier_basis()] returns. It
#' adds no property to [basis] and exists as the class the trigonometric
#' methods dispatch on. A Fourier basis holds a constant function and pairs of
#' sines and cosines of increasing frequency, and answers every derivative and
#' its integral from one identity.
#'
#' @details
#' # One identity for every order
#'
#' Writing \eqn{z = 2\pi (x - \ell)/\omega} for the phase and \eqn{\omega}
#' for the period,
#' \deqn{\frac{\mathrm{d}^{k}}{\mathrm{d}x^{k}} \sin(j z)
#'       = \left(\frac{2\pi j}{\omega}\right)^{k}
#'         \sin\!\left(j z + \frac{k\pi}{2}\right),}
#' and the same for the cosine. Differentiating a sinusoid shifts its phase by
#' a quarter turn and multiplies it by its frequency, so no order is a special
#' case and the fourth derivative costs exactly what the first does.
#'
#' The identity holds for negative \eqn{k}, which is where the antiderivative
#' comes from: [basis_int.FourierBasis()] evaluates it at \eqn{k = -1} and
#' subtracts the value at the lower endpoint.
#'
#' # Its `basis_params`
#'
#' Three entries. `omega` is the period, `n_pairs` is
#' `(dimension - 1) %/% 2`, and `full_period` records whether `omega` equals
#' the width of the interval. The last decides which route
#' [basis_gram.FourierBasis()] takes: a closed diagonal matrix when it is
#' `TRUE`, a quadrature when it is not.
#'
#' @inheritParams basis
#'
#' @return An object of class `FourierBasis`, inheriting from [basis], with
#'   the same five properties and `basis_params` holding `omega`, `n_pairs`
#'   and `full_period`. Call [fourier_basis()] instead of this class directly;
#'   it rejects an even dimension, defaults the period and records all three.
#'
#' @seealso [fourier_basis()], the constructor;
#'   [basis_eval.FourierBasis()], [basis_deriv.FourierBasis()],
#'   [basis_int.FourierBasis()] and [basis_gram.FourierBasis()] for the
#'   methods registered on it.
#'
#' @examples
#' f <- fourier_basis(dimension = 5)
#' S7::S7_inherits(f, FourierBasis)
#' f@basis_params
#'
#' # Differentiating shifts the phase by a quarter turn: the first derivative
#' # of sin1 at the start of the period is its own frequency, 2*pi.
#' basis_deriv(f, 0, order = 1)
#' 2 * pi
#'
#' @export
FourierBasis <- S7::new_class("FourierBasis", parent = basis)


#' Construct a Fourier Basis
#'
#' @description
#' Returns a basis of a constant function and \eqn{(K-1)/2} sine-cosine pairs
#' of increasing frequency on \eqn{[\ell, u]}. It is the basis to reach for
#' when the function being modeled is periodic or nearly so: over a whole
#' period the functions are mutually orthogonal, and stay orthogonal after
#' differentiation, so the Gram matrix is diagonal at every order and a
#' roughness penalty is one vector of numbers.
#'
#' @details
#' # What the columns are
#'
#' With \eqn{z = 2\pi (x - \ell)/\omega}, the columns are \eqn{1},
#' \eqn{\sin z}, \eqn{\cos z}, \eqn{\sin 2z}, \eqn{\cos 2z}, and so on
#' up to `n_pairs` frequencies, named `const`, `sin1`, `cos1`, `sin2`, `cos2`.
#'
#' # Why the dimension must be odd
#'
#' A sine without its cosine represents a wave at one phase but not at the
#' next, so the fitted function would depend on where the interval was cut.
#' An even dimension throws, naming the two odd numbers on either side. It is
#' not adjusted silently: growing it would return a basis of a size the caller
#' did not ask for, and the constructor is the only place the mismatch can be
#' caught.
#'
#' # The period, and what changes when it is not the interval
#'
#' `omega` defaults to `upper - lower`, so the interval is exactly one period
#' and the functions are orthogonal on it. The order-\eqn{d} Gram matrix is
#' then diagonal with entries \eqn{\omega} for the constant at order 0, zero
#' for it above, and \eqn{(\omega/2)(2\pi j/\omega)^{2d}} for both members
#' of pair \eqn{j}, written in closed form.
#'
#' Any other positive period is accepted and every generic still answers, but
#' the interval is no longer a whole number of periods, the orthogonality
#' fails, and the Gram matrix is computed by composite Gauss-Legendre instead.
#' `basis_params$full_period` records which case the object is in, and
#' [basis_is_numerical()] reports `basis_gram` as `TRUE` there: the family
#' says so itself through [basis_numerical_route.FourierBasis()], where
#' reading which class the method is registered on would answer `FALSE`,
#' the owner being `FourierBasis` either way.
#'
#' # Periodic by construction
#'
#' Every column except the constant takes the same value at both ends of a
#' full period, so a fitted curve joins up. That is the property to want here,
#' and the reason not to reach for [bspline_basis()], whose ends are free.
#'
#' @param lower,upper The endpoints of the interval, each a single finite
#'   number with `lower < upper`. Default \eqn{[0, 1]}. Evaluating outside
#'   throws.
#' @param dimension The number of basis functions, a single **odd** whole
#'   number of at least 1, default `5`. `1` is the constant alone, `5` is the
#'   constant and two pairs. An even value throws.
#' @param omega The period, a single positive finite number. `NULL`, the
#'   default, uses `upper - lower`, the value that makes the basis orthogonal
#'   on its interval. Anything not a single positive number throws.
#'
#' @return An object of class [FourierBasis], with `basis_name` `"fourier"`,
#'   `basis_params` holding `omega`, `n_pairs` and `full_period`, and column
#'   names `const`, `sin1`, `cos1`, and so on.
#'
#' @seealso [bspline_basis()] for a local basis with free ends and
#'   [poly_basis()] for a global polynomial one; [basis_gram()] for the
#'   roughness penalty this basis diagonalizes.
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#' b
#' basis_colnames(b)
#'
#' # Orthogonal on a whole period, so the Gram matrix is diagonal: omega for
#' # the constant, omega/2 for every sine and cosine.
#' round(basis_gram(b), 10)
#'
#' # It stays diagonal after differentiating, with entries (2 pi j)^(2d) / 2.
#' round(diag(basis_gram(b, order = 2)), 4)
#' c(0, rep((2 * pi * (1:2))^4 / 2, each = 2))
#'
#' # Periodic: the ends of a full period agree.
#' max(abs(basis_eval(b, 0) - basis_eval(b, 1)))
#'
#' # An even dimension would leave half a pair, and is rejected.
#' try(fourier_basis(dimension = 4))
#'
#' # A period that is not the interval width gives up the orthogonality, and
#' # the Gram matrix is then a quadrature and no longer diagonal.
#' f <- fourier_basis(dimension = 5, omega = 0.7)
#' f@basis_params$full_period
#' round(basis_gram(f), 4)
#'
#' @export
fourier_basis <- function(lower = 0, upper = 1, dimension = 5, omega = NULL) {
  dimension <- check_basis_args(lower, upper, dimension)

  if (dimension %% 2L == 0L) {
    stop(sprintf(
      paste0(
        "'dimension' must be odd: a Fourier basis holds a constant plus ",
        "complete sine-cosine pairs, so %d would leave half a pair. Use %d ",
        "or %d."
      ),
      dimension, dimension - 1L, dimension + 1L
    ), call. = FALSE)
  }

  if (is.null(omega)) {
    omega <- upper - lower
  } else if (!is.numeric(omega) || length(omega) != 1L || !is.finite(omega) ||
    omega <= 0) {
    stop("'omega' must be a single positive number.", call. = FALSE)
  }

  FourierBasis(
    basis_name = "fourier",
    dimension = dimension,
    lower = lower,
    upper = upper,
    basis_params = list(
      omega = omega,
      n_pairs = (dimension - 1L) %/% 2L,
      full_period = isTRUE(all.equal(omega, upper - lower))
    )
  )
}


#' Column Names of a Fourier Basis
#'
#' @name basis_colnames.FourierBasis
#'
#' @description
#' Names the columns `const`, `sin1`, `cos1`, `sin2`, `cos2`, and so on: the
#' constant first, then the sine and cosine of each frequency in the order the
#' matrix holds them. A coefficient's name therefore says which harmonic it
#' belongs to, where the default `fo1 ... fo5` of [basis_colnames.basis()]
#' would leave the reader counting.
#'
#' @details
#' At `n_pairs == 0`, which is `dimension = 1`, the answer is the single name
#' `"const"` and the general branch is skipped. Falling through it would give
#' `"sin"` and `"cos"` with no number, `paste0()` recycling a zero-length
#' argument to the empty string.
#'
#' @param basis A [FourierBasis] object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A character vector of length `basis@dimension`, `"const"` first.
#'
#' @seealso [basis_colnames()] for the generic.
#'
#' @keywords internal
S7::method(basis_colnames, FourierBasis) <- function(basis, ...) {
  n_pairs <- basis@basis_params$n_pairs
  # With no pairs the basis is the constant alone. Falling through would give
  # "sin" and "cos" without a number, because paste0() recycles a zero-length
  # argument to the empty string rather than to nothing.
  if (n_pairs == 0L) return("const")
  j <- seq_len(n_pairs)
  c("const", as.character(rbind(paste0("sin", j), paste0("cos", j))))
}


#' Evaluate a Fourier Basis
#'
#' @name basis_eval.FourierBasis
#'
#' @description
#' Evaluates the constant and the sine-cosine pairs at the given points, from
#' the phase-shift identity of [FourierBasis] taken at order zero, where it is
#' the sinusoid itself. The first column is `1` at every point; column
#' \eqn{2j} is \eqn{\sin(j z)} and column \eqn{2j + 1} is \eqn{\cos(j z)},
#' with \eqn{z = 2\pi (x - \ell)/\omega}.
#'
#' @details
#' Cost is two `sin()` and `cos()` calls per frequency per point, with no
#' recurrence and no accumulation, so a high frequency is evaluated as
#' accurately as a low one. Missing points give missing rows.
#'
#' @param basis A [FourierBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, with column names `const`, `sin1`, `cos1`, and so on.
#'
#' @seealso [fourier_trig()], which builds the trigonometric columns;
#'   [basis_eval()] for the generic.
#'
#' @keywords internal
S7::method(basis_eval, FourierBasis) <- function(basis, x, ...) {
  out <- cbind(rep(1, length(x)), fourier_trig(basis, x, 0L))
  out[is.na(x), ] <- NA_real_
  name_columns(out, basis)
}


#' Derivatives of a Fourier Basis
#'
#' @name basis_deriv.FourierBasis
#'
#' @description
#' Returns the `order`-th derivative of every column, exactly and at any
#' order, from the phase-shift identity of [FourierBasis]: differentiating
#' \eqn{\sin(jz)} shifts its phase by \eqn{k\pi/2} and multiplies it by
#' \eqn{(2\pi j/\omega)^{k}}. The constant column is zero at every order
#' above 0.
#'
#' @details
#' Reaching order \eqn{k} costs the same as reaching order 1: the shift and
#' the scale are both computed directly from \eqn{k}, with no recursion over
#' the orders below it. A Fourier basis therefore has no order at which its
#' derivatives stop being available, in contrast with a spline, whose
#' derivatives run out at its degree.
#'
#' The scale grows as \eqn{j^{k}}, so a high frequency differentiated many
#' times is a large number: at `dimension = 21` and `order = 4` the largest
#' entry is \eqn{(20\pi)^4}, about 1.6e+06 on the unit interval. That is the
#' value, not a loss of accuracy.
#'
#' @param basis A [FourierBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param order The derivative order, a single non-negative whole number,
#'   default `1`. Any order is available.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, first column zero, with column names `const`, `sin1`, and so on.
#'
#' @seealso [fourier_trig()], which applies the identity;
#'   [basis_deriv()] for the generic.
#'
#' @keywords internal
S7::method(basis_deriv, FourierBasis) <- function(basis, x, order = 1L, ...) {
  out <- cbind(rep(0, length(x)), fourier_trig(basis, x, order))
  out[is.na(x), ] <- NA_real_
  name_columns(out, basis)
}


#' Integral of a Fourier Basis
#'
#' @name basis_int.FourierBasis
#'
#' @description
#' Returns \eqn{\int_{\ell}^{x} \varphi_j(t)\,\mathrm{d}t} for every column
#' in closed form, with no quadrature. The constant integrates to
#' \eqn{x - \ell} and the sinusoids come from the phase-shift identity of
#' [FourierBasis] taken at order \eqn{-1}.
#'
#' @details
#' The identity at \eqn{k = -1} gives *an* antiderivative, which is not the
#' one [basis_int()] promises. The two differ by a constant that is not the
#' same in every column: at the lower endpoint the sine columns of the raw
#' antiderivative are \eqn{-\omega/(2\pi j)} while the cosine columns are
#' already zero. Subtracting the row at the lower endpoint corrects every
#' column at once and makes the anchoring exact.
#'
#' Over a full period every sinusoid integrates to zero, so the row at
#' `basis@upper` is \eqn{(\omega, 0, 0, \ldots)}, which is the area under a
#' fitted curve being the constant's coefficient times the period.
#'
#' @param basis A [FourierBasis] object.
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension`
#'   columns, exactly zero in the row at `basis@lower`.
#'
#' @seealso [basis_int()] for the generic and the anchoring convention;
#'   [fourier_trig()], which supplies the raw antiderivative.
#'
#' @keywords internal
S7::method(basis_int, FourierBasis) <- function(basis, x, ...) {
  anti <- fourier_trig(basis, x, -1L)
  at_lower <- fourier_trig(basis, basis@lower, -1L)
  trig <- sweep(anti, 2L, at_lower, "-")
  out <- cbind(x - basis@lower, trig)
  out[is.na(x), ] <- NA_real_
  name_columns(out, basis)
}


#' Gram Matrix of a Fourier Basis
#'
#' @name basis_gram.FourierBasis
#'
#' @description
#' Returns the inner products of the `order`-th derivatives. On a whole period
#' the matrix is diagonal and written down in closed form; on any other period
#' the orthogonality fails and the method computes a quadrature instead.
#'
#' @details
#' # The closed form
#'
#' Over a whole period the trigonometric functions are mutually orthogonal,
#' and stay orthogonal after differentiation, a derivative only shifting the
#' phase and rescaling. The order-\eqn{d} matrix is therefore diagonal, with
#'
#' \deqn{G_{00} = \omega \ (d = 0), \quad 0 \ (d \ge 1), \qquad
#'   G_{jj} = \frac{\omega}{2}\left(\frac{2\pi j}{\omega}\right)^{2d}}
#'
#' for both members of pair \eqn{j}. At `dimension = 7` on \eqn{[0, 1]} the
#' order-2 diagonal is `0, 779.27, 779.27, 12468.4, 12468.4, 63121.1,
#' 63121.1`, matching \eqn{(2\pi j)^4/2} exactly.
#'
#' # The other period
#'
#' When `basis_params$full_period` is `FALSE` the method delegates to
#' [numerical_gram()], composite Gauss-Legendre over the interval, and the
#' result is a full matrix.
#'
#' [basis_is_numerical()] reports `basis_gram` as `TRUE` on such a basis,
#' the family answering through [basis_numerical_route.FourierBasis()]
#' rather than through the class the method is registered on, which is
#' `FourierBasis` in both branches. [check_basis()] reads the same predicate,
#' so it holds this matrix to the tolerance a quadrature deserves.
#'
#' @param basis A [FourierBasis] object.
#' @param order The derivative order, a single non-negative whole number,
#'   default `0`.
#' @param at,weight Handled in the body of [basis_gram()] before dispatch, so
#'   they never arrive here. Named only because S7 requires a method's formals
#'   to contain the generic's.
#' @param ... Passed to [numerical_gram()] on the non-full-period branch,
#'   where `panels` and `nodes` control the quadrature. Ignored on the closed
#'   branch.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with column names `const`, `sin1`, and so on. Diagonal when the interval
#'   is one whole period, full otherwise; singular for any `order >= 1`, the
#'   constant differentiating away.
#'
#' @seealso [basis_gram()] for the generic and the alternative measures;
#'   [numerical_gram()] for the branch taken on another period.
#'
#' @keywords internal
S7::method(basis_gram, FourierBasis) <- function(basis, order = 0L, at = NULL,
                                                 weight = NULL, ...) {
  p <- basis@basis_params
  if (!p$full_period) {
    return(numerical_gram(basis, order, ...))
  }

  omega <- p$omega
  j <- seq_len(p$n_pairs)
  scale <- (2 * pi * j / omega)^(2 * order)
  const <- if (order == 0L) omega else 0
  d <- c(const, rep(omega / 2 * scale, each = 2L))

  nm <- basis_colnames(basis)
  matrix(diag(d, nrow = basis@dimension), basis@dimension, basis@dimension,
    dimnames = list(nm, nm)
  )
}


#' The Trigonometric Columns of a Fourier Basis
#'
#' @description
#' Builds the sine and cosine columns at derivative order `d` from the
#' phase-shift identity
#' \deqn{\frac{\mathrm{d}^{d}}{\mathrm{d}x^{d}} \sin(jz)
#'   = \left(\frac{2\pi j}{\omega}\right)^{d}
#'     \sin\!\left(jz + \frac{d\pi}{2}\right),}
#' the one routine behind [basis_eval.FourierBasis()] at `d = 0`,
#' [basis_deriv.FourierBasis()] above it and [basis_int.FourierBasis()] at
#' `d = -1`. The constant column is not included; each caller prepends its
#' own.
#'
#' @details
#' At `d = -1` the identity gives an antiderivative whose constant is not the
#' one [basis_int()] promises; the caller subtracts the row at the lower
#' endpoint. Columns are interleaved sine-then-cosine per frequency, matching
#' [basis_colnames.FourierBasis()].
#'
#' At `n_pairs == 0` the result is a `length(x)` by 0 matrix, and the loop is
#' skipped.
#'
#' @param basis A [FourierBasis] object.
#' @param x A numeric vector of evaluation points. Not range-checked here.
#' @param d The order, a single whole number that may be negative. `-1` gives
#'   an antiderivative, `0` the functions themselves.
#'
#' @return A numeric matrix with `length(x)` rows and `basis@dimension - 1`
#'   columns, no dimnames, sine and cosine interleaved by frequency.
#'
#' @seealso [basis_eval.FourierBasis()], [basis_deriv.FourierBasis()] and
#'   [basis_int.FourierBasis()], its three callers.
#'
#' @keywords internal
fourier_trig <- function(basis, x, d) {
  p <- basis@basis_params
  n_pairs <- p$n_pairs
  out <- matrix(0, length(x), 2L * n_pairs)
  if (n_pairs == 0L) return(out)

  z <- 2 * pi * (x - basis@lower) / p$omega
  shift <- d * pi / 2
  for (j in seq_len(n_pairs)) {
    scale <- (2 * pi * j / p$omega)^d
    out[, 2L * j - 1L] <- sin(j * z + shift) * scale
    out[, 2L * j] <- cos(j * z + shift) * scale
  }
  out
}


#' @title The Gram Route of a Fourier Basis
#' @name basis_numerical_route.FourierBasis
#' @description
#' Reports `basis_gram` as `TRUE` when `basis_params$full_period` is `FALSE`,
#' where the owner test would read `FourierBasis` and answer `FALSE`.
#' [basis_gram.FourierBasis()] delegates to [numerical_gram()] in that case, so
#' the matrix is a composite Gauss-Legendre quadrature and carries its error.
#' The evaluation, the derivatives and the anchored integral are closed form at
#' any period and are left as the owner test finds them.
#'
#' What this buys is that [check_basis()] holds the Gram matrix to the tolerance
#' a quadrature deserves rather than the one meant for a closed form, and that
#' [print.basis()] names the route in use.
#'
#' @param basis A [FourierBasis] object.
#' @param ... Unused, and accepted so the signature matches the generic's.
#'
#' @return The named logical vector [basis_numerical_route()] describes, with
#'   `basis_gram` `TRUE` for a basis whose period is not the interval width.
#'
#' @seealso [basis_gram.FourierBasis()] for the two branches, and
#'   [basis_numerical_route.basis()] for the default this starts from.
#'
#' @keywords internal
S7::method(basis_numerical_route, FourierBasis) <- function(basis, ...) {
  out <- route_by_owner(basis)
  out[["basis_gram"]] <- out[["basis_gram"]] || !basis@basis_params$full_period
  out
}
