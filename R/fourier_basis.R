#' @include numerical_fallbacks.R
NULL


#' Fourier Basis
#'
#' @description
#' The S7 class of Fourier bases. Constructed by \code{\link{fourier_basis}}.
#'
#' @details
#' The basis holds a constant function and pairs of sines and cosines of
#' increasing frequency. Every derivative and the antiderivative come from one
#' identity, so no order is special:
#' \deqn{\frac{\mathrm{d}^{k}}{\mathrm{d}x^{k}} \sin(j z)
#'       = \left(\frac{2\pi j}{\omega}\right)^{k} \sin\!\left(j z + \frac{k\pi}{2}\right),
#'       \qquad z = \frac{2\pi (x - \ell)}{\omega},}
#' and likewise for the cosine. The identity holds for negative \eqn{k} as
#' well, which is where the antiderivative comes from.
#'
#' @inheritParams basis
#'
#' @return An object of class \code{FourierBasis}. Use
#'   \code{\link{fourier_basis}} rather than calling the class directly, so
#'   that the dimension is checked and the period recorded.
#'
#' @seealso \code{\link{fourier_basis}}
#'
#' @examples
#' f <- fourier_basis(dimension = 5)
#' S7::S7_inherits(f, FourierBasis)
#'
#' @export
FourierBasis <- S7::new_class("FourierBasis", parent = basis)


#' Construct a Fourier Basis
#'
#' @description
#' A basis of a constant function and \eqn{(K-1)/2} sine-cosine pairs of
#' increasing frequency on \eqn{[\ell, u]}.
#'
#' @details
#' The dimension must be odd. A sine without its cosine is a basis that can
#' represent a wave at one phase and not at another, so a dimension that would
#' leave a half pair is rejected rather than adjusted: growing it silently would
#' return a basis of a size the caller did not ask for, and the constructor is
#' the only place the inconsistency can be caught.
#'
#' The period \code{omega} defaults to the width of the interval, which is what
#' makes the basis functions orthogonal on it. A different period is accepted
#' and everything still works, but the interval is then no longer a whole
#' number of periods, so the Gram matrix stops being diagonal and is computed
#' by quadrature instead of in closed form.
#'
#' @param lower,upper The endpoints of the interval.
#' @param dimension The number of basis functions, an odd positive integer.
#' @param omega The period. Defaults to \code{upper - lower}.
#'
#' @return An object of class \code{\link{FourierBasis}}.
#'
#' @seealso \code{\link{bspline_basis}}, \code{\link{check_basis}}
#'
#' @examples
#' b <- fourier_basis(dimension = 5)
#' b
#' basis_colnames(b)
#'
#' # orthogonal on a whole period, so the Gram matrix is diagonal
#' round(basis_gram(b), 10)
#'
#' # an even dimension would leave half a pair, and is rejected
#' try(fourier_basis(dimension = 4))
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
#' @description
#' The constant followed by the sine and cosine of each frequency.
#' @param basis A \code{\link{FourierBasis}} object.
#' @param ... Unused.
#' @return A character vector of length \code{basis@dimension}.
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
#' @description
#' The constant and the sine-cosine pairs, from the shift identity of
#' \code{\link{FourierBasis}} at order zero.
#' @param basis A \code{\link{FourierBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_eval, FourierBasis) <- function(basis, x, ...) {
  out <- cbind(rep(1, length(x)), fourier_trig(basis, x, 0L))
  out[is.na(x), ] <- NA_real_
  name_columns(out, basis)
}


#' Derivatives of a Fourier Basis
#'
#' @name basis_deriv.FourierBasis
#' @description
#' Exact derivatives of any order, from the shift identity of
#' \code{\link{FourierBasis}}. The constant differentiates to zero.
#' @param basis A \code{\link{FourierBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param order The derivative order.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_deriv, FourierBasis) <- function(basis, x, order = 1L, ...) {
  out <- cbind(rep(0, length(x)), fourier_trig(basis, x, order))
  out[is.na(x), ] <- NA_real_
  name_columns(out, basis)
}


#' Integral of a Fourier Basis
#'
#' @name basis_int.FourierBasis
#' @description
#' The exact integral from the lower endpoint.
#' @details
#' The shift identity at order \eqn{-1} gives an antiderivative, not the one
#' this package's convention asks for, and the two differ by a constant that is
#' not the same for every column: at the lower endpoint the sine columns of the
#' raw antiderivative are \eqn{-\omega/(2\pi j)} while the cosine columns are
#' already zero. Subtracting the value at the lower endpoint fixes every column
#' at once and makes the convention hold exactly.
#' @param basis A \code{\link{FourierBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param ... Unused.
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension} columns.
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
#' @description
#' The inner products of the basis functions, diagonal in closed form when the
#' interval is a whole period.
#' @details
#' Over a whole period the trigonometric functions are mutually orthogonal, and
#' remain so after differentiation, since a derivative only shifts the phase and
#' rescales. The order-\eqn{d} Gram matrix is therefore diagonal, with entries
#' \eqn{\omega} for the constant at order zero, zero for it at every higher
#' order, and \eqn{(\omega/2)(2\pi j/\omega)^{2d}} for both members of pair
#' \eqn{j}.
#'
#' When \code{omega} is not the width of the interval the orthogonality fails,
#' and the numerical method of the \code{\link{basis}} class is used instead.
#' @param basis A \code{\link{FourierBasis}} object.
#' @param order The derivative order.
#' @param ... Passed to the numerical method when the period is not the
#'   interval width.
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
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
#' The sine and cosine columns at derivative order \code{d}, from the shift
#' identity. Order \eqn{-1} gives an antiderivative.
#'
#' @param basis A \code{\link{FourierBasis}} object.
#' @param x A numeric vector of evaluation points.
#' @param d The order, which may be negative.
#'
#' @return A numeric matrix with \code{length(x)} rows and
#'   \code{basis@dimension - 1} columns.
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
