#' @include generics.R ggplot.R
NULL


#' Print a Basis
#'
#' @name print.basis
#'
#' @description
#' Prints a summary of a basis object: the family from which it comes, how many
#' functions it holds and over how many variables, the interval over which each
#' variable runs, the parameters that the family carries, and which of the three
#' derived quantities are computed numerically instead of from a formula. The
#' summary has four lines, plus a header and one line per parameter when the
#' basis carries parameters.
#'
#' @details
#' # The lines
#'
#' `Basis:` is `@basis_name`, which records how the object was built. A
#' wrapper or a product names its parents, so an orthonormalized B-spline
#' reads `orthonorm(bspline)` and a two-way product reads
#' `tensor(bspline, fourier)`.
#'
#' `Functions:` is `@dimension`, the number of columns that [basis_eval()]
#' returns, and `Variables:` is [basis_nvar()], the number of columns of `x`
#' that it expects. Every family here has one variable except
#' [tensor_basis()].
#'
#' `Domain:` is `@lower` and `@upper` as a closed interval per variable,
#' separated by `x` for a product. It is the interval against which every
#' generic checks its evaluation points. A point outside signals an error that
#' states how many of the points were outside and what the interval is.
#'
#' `Parameters:` appears only when `@basis_params` is non-empty, and prints
#' one line per entry. A numeric entry of more than four values is abbreviated
#' to its length, so a B-spline over many knots reads `<8 values>` and leaves
#' the knots to `b@basis_params$knots`.
#'
#' `Numerical:` reads [basis_is_numerical()] and names those of `basis_deriv`,
#' `basis_int` and `basis_gram` that are computed by the numerical route:
#' finite differences for the derivatives, quadrature for the integral and the
#' Gram matrix. The shipped families report `none`, except a Fourier
#' basis whose interval is not a whole period, whose Gram matrix is computed
#' by quadrature. A basis defined from its evaluation alone reports all
#' three.
#'
#' @param x A basis object, of any class inheriting from [basis].
#' @param ... Unused, and accepted so that the signature matches the `print()`
#'   generic's.
#'
#' @return `x`, invisibly.
#'
#' @seealso [basis_is_numerical()] for the last line as a named logical vector,
#'   [plot.basis()] to see the functions themselves, and [check_basis()] for a
#'   verification of the components that this summary only names.
#'
#' @examples
#' # The parameter block differs by family: the degree and the knots for a
#' # B-spline, the period, the number of pairs and a full-period flag for
#' # Fourier, the degree alone for Legendre.
#' bspline_basis(dimension = 6)
#' fourier_basis(dimension = 5)
#' poly_basis(dimension = 4)
#'
#' # A product names its margins and reports one interval per variable.
#' tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#'
#' # More than four values in a parameter are abbreviated to a count. The
#' # values themselves stay reachable on the object.
#' b <- bspline_basis(dimension = 12)
#' b
#' b@basis_params$knots
#'
#' # A basis defined from its evaluation alone reports all three quantities
#' # as numerical, where every shipped family reports none.
#' Bumps <- S7::new_class("Bumps", parent = basis)
#' S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
#'   out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
#'                           "-")^2 / 0.12^2)
#'   colnames(out) <- basis_colnames(basis)
#'   out
#' }
#' Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)
S7::method(print, basis) <- function(x, ...) {
  cat("Basis: ", x@basis_name, "\n", sep = "")
  cat("Functions: ", x@dimension, "   Variables: ", basis_nvar(x), "\n",
    sep = ""
  )
  cat("Domain: ",
    paste(sprintf("[%s, %s]", format(x@lower), format(x@upper)),
      collapse = " x "
    ),
    "\n",
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
#'
#' @description
#' `autoplot()` returns a graph of all `@dimension` functions of a basis on
#' one panel, over an equally spaced grid covering the whole interval, and
#' `plot()` prints that graph. `order` selects what is drawn: the functions
#' themselves, a derivative of any order, or the integral anchored at the
#' lower endpoint. One line is drawn per basis function, each in its own
#' colour. The graph has no legend, and its title names the basis and the
#' number of its functions.
#'
#' @details
#' # What is drawn
#'
#' The grid is `n` equally spaced points from `@lower` to `@upper`, endpoints
#' included, and the curves are whichever of [basis_eval()], [basis_deriv()]
#' or [basis_int()] is the one that `order` names. The vertical axis is labeled
#' to match:
#' \eqn{B(x)} at order `0`, \eqn{B'(x)} at order `1`, \eqn{B^{(k)}(x)} above
#' that, and \eqn{\int B(t)\,\mathrm{d}t} for the integral.
#'
#' The graph is drawn with \pkg{ggplot2}, which is suggested; the methods
#' signal an error that names it when it is not installed. The title, the
#' labels of the axes and the colours are changed on the returned object, for
#' instance with `+ ggplot2::labs(title = "my title")`.
#'
#' # Only one variable
#'
#' A basis of several variables signals an error. A product of two bases is
#' a surface over a rectangle; a margin can be plotted instead, which is
#' `tb@marginals[[1]]` for a [tensor_basis()].
#'
#' # Orders above the smoothness of the family
#'
#' An order above the smoothness of the family is accepted and draws a flat
#' line at zero, because the fourth derivative of a piecewise cubic is zero
#' away from the knots.
#'
#' @param x,object A basis object of one variable, of any class inheriting
#'   from [basis]. A basis of several variables signals an error.
#' @param order What to draw. `0`, the default, draws the basis functions; a
#'   positive whole number draws that derivative; `-1` draws the integral from
#'   the lower endpoint. Any other value signals an error, including a
#'   fraction, a value below `-1` and a vector of length other than one.
#' @param n The number of grid points, default `200`. A basis with many
#'   knots or a high frequency needs more, because 200 points leave a curve
#'   visibly polygonal there.
#' @param ... Unused; an argument passed here signals an error.
#'
#' @return `autoplot()` returns a ggplot2 object; `plot()` prints the same
#'   object and returns it invisibly.
#'
#' @seealso [basis_eval()], [basis_deriv()] and [basis_int()] for the numbers
#'   behind the three cases of `order`, and [print.basis()] for the object's
#'   summary.
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' b <- bspline_basis(dimension = 6)
#'
#' # The three graphs show the six cubic B-splines, their first derivative,
#' # and their integrals.
#' plot(b)
#' plot(b, order = 1)
#' plot(b, order = -1)
#'
#' # Every curve in the third graph starts at zero, because the integral is
#' # anchored at the lower endpoint.
#' basis_int(b, 0)
#'
#' # The returned graph is modified like any other ggplot2 object.
#' p <- plot(b)
#' p + ggplot2::labs(title = "six cubic B-splines on [0, 1]")
#'
#' # A Fourier basis at a high frequency needs a finer grid than the default.
#' plot(fourier_basis(dimension = 21), n = 1000)
#' @export
plot.basis <- function(x, order = 0L, n = 200L, ...) {
  p <- autoplot_basis(x, order = order, n = n, ...)
  print(p)
  invisible(p)
}
S7::method(plot, basis) <- plot.basis

#' @rdname plot.basis
autoplot_basis <- function(object, order = 0L, n = 200L, ...) {
  if (...length()) {
    stop("The graph of a basis takes no further arguments; modify the ",
         "returned object instead.", call. = FALSE)
  }
  x <- object
  if (basis_nvar(x) > 1L) {
    stop(
      "plot() draws a basis of one variable. A product of several has no ",
      "single picture: draw its marginals, or a surface of one function.",
      call. = FALSE
    )
  }
  if (!is.numeric(order) || length(order) != 1L || order < -1 ||
    order != round(order)) {
    stop(
      "'order' must be -1 (integral), 0 (the basis), or a positive integer.",
      call. = FALSE
    )
  }
  need_ggplot()
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
    as.expression(bquote(B^(.(order)) * (x)))
  )

  fn <- basis_colnames(x)
  d <- data.frame(x = rep(grid, ncol(y)), y = as.vector(y),
                  fn = factor(rep(fn, each = n), levels = fn))
  .data <- ggplot2::.data
  ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y,
                                  colour = .data$fn)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(title = sprintf("%s basis, %d functions", x@basis_name,
                                  ncol(y)),
                  x = "x", y = ylab, colour = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

S7::method(autoplot, basis) <- autoplot_basis
