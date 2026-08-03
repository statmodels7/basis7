#' @include generics.R
NULL


#' Print a Basis
#'
#' @name print.basis
#' @description
#' Reports the family, the number of functions, the interval, any parameters
#' the family carries, and which of the derived quantities are computed
#' numerically.
#' @param x An object inheriting from class \code{basis}.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' bspline_basis(dimension = 6)
#' fourier_basis(dimension = 5)
S7::method(print, basis) <- function(x, ...) {
  cat("Basis: ", x@basis_name, "\n", sep = "")
  cat("Functions: ", x@dimension,
    "   Interval: [", format(x@lower), ", ", format(x@upper), "]\n",
    sep = ""
  )

  p <- x@basis_params
  if (length(p)) {
    shown <- vapply(p, function(v) {
      if (is.numeric(v) && length(v) > 4L) {
        sprintf("<%d values>", length(v))
      } else if (is.numeric(v)) {
        paste(format(v, digits = 4), collapse = ", ")
      } else {
        paste(format(v), collapse = ", ")
      }
    }, character(1))
    cat("Parameters:\n")
    cat(paste0("  ", format(names(shown)), "  ", shown, collapse = "\n"), "\n",
      sep = ""
    )
  }

  num <- basis_is_numerical(x)
  cat("Numerical: ",
    if (any(num)) paste(names(num)[num], collapse = ", ") else "none",
    "\n",
    sep = ""
  )
  invisible(x)
}


#' Plot a Basis
#'
#' @name plot.basis
#' @description
#' Draws every basis function over the interval, or its derivative or integral.
#' @param x An object inheriting from class \code{basis}.
#' @param order What to draw: \code{0} for the basis functions, a positive
#'   integer for that derivative, \code{-1} for the integral from the lower
#'   endpoint.
#' @param n The number of points at which to evaluate.
#' @param ... Passed to \code{\link[graphics]{matplot}}.
#' @return \code{x}, invisibly.
#' @examples
#' b <- bspline_basis(dimension = 6)
#' plot(b)
#' plot(b, order = 1)
#' plot(b, order = -1)
S7::method(plot, basis) <- function(x, order = 0L, n = 200L, ...) {
  if (!is.numeric(order) || length(order) != 1L || order < -1 ||
    order != round(order)) {
    stop(
      "'order' must be -1 (integral), 0 (the basis), or a positive integer.",
      call. = FALSE
    )
  }
  order <- as.integer(order)

  grid <- seq(x@lower, x@upper, length.out = n)
  y <- if (order == -1L) {
    basis_int(x, grid)
  } else {
    basis_deriv(x, grid, order = order)
  }

  ylab <- switch(as.character(min(order, 2L)),
    "-1" = expression(integral(B(t) * dt)),
    "0" = expression(B(x)),
    "1" = expression(B * minute * (x)),
    bquote(B^(.(order)) * (x))
  )

  old <- graphics::par(mar = c(4.5, 5, 3, 1))
  on.exit(graphics::par(old))
  graphics::matplot(grid, y,
    type = "l", lty = 1, xlab = "x", ylab = ylab,
    main = x@basis_name, ...
  )
  invisible(x)
}
