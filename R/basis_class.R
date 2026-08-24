#' Basis Expansion
#'
#' @description
#' The abstract parent of every basis in the package, and the class to inherit
#' from when writing one of your own. A basis is a finite collection of
#' functions on an interval; the object carries the interval, how many
#' functions there are, and whatever else the family needs, and the generics
#' evaluate that collection, differentiate it, integrate it and take its inner
#' products. The class is abstract, so `basis(...)` throws; construct one of
#' the concrete families or a subclass.
#'
#' @details
#' # What a basis is
#'
#' A basis of dimension \eqn{d} on \eqn{[a, b]} is a collection
#' \eqn{\varphi_1, \dots, \varphi_d}, and the object exists so that a
#' function may be written as a linear combination of them,
#'
#' \deqn{f(x) = \sum_{j=1}^{d} \beta_j \varphi_j(x) = B(x)\beta,
#'   \qquad B(x)_{ij} = \varphi_j(x_i),}
#'
#' with \eqn{B(x)} the \eqn{n \times d} design matrix [basis_eval()]
#' returns. Fitting \eqn{f} is a linear problem in \eqn{\beta} whatever the
#' family, so the derivative, the anchored integral and the Gram matrix are
#' all properties of the basis alone and can be computed once, before any
#' data arrive.
#'
#' # Writing a subclass
#'
#' A concrete basis is a subclass of this one. It must implement
#' [basis_eval()]; [basis_deriv()], [basis_int()] and [basis_gram()] have
#' numerical methods registered on this class, so a subclass supplying its
#' evaluation alone answers all four generics immediately. A closed form
#' registered later takes over through dispatch, with no change to calling
#' code, and [basis_is_numerical()] reports which of the three are still on
#' the fallback. The vignette `vignette("defining-a-basis")` works one
#' through.
#'
#' # Bases here are complete
#'
#' A B-spline basis carries all its functions, so its rows sum to one and it
#' spans the constant. Restricting a basis, for identifiability or to separate
#' a linear part from a nonlinear one, is a linear transformation of it:
#' [constrain_basis()], [orthonorm_basis()] and [dr_basis()] each return a
#' [TransformedBasis] that is itself a basis, so nothing downstream has to
#' know a restriction happened.
#'
#' # One variable or several
#'
#' A basis lives on an interval, or, when it is a product of several, on a
#' box. `@lower` and `@upper` then hold one endpoint per variable and
#' [basis_nvar()] reports how many; [basis_eval()] takes a matrix of that
#' many columns instead of a vector. Nothing else changes, and a basis of one
#' variable is the case \eqn{d = 1} of the same object.
#'
#' # What the validator enforces
#'
#' `@basis_name` must be one string, `@dimension` one integer of at least 1,
#' and `@lower` and `@upper` must be finite, of equal length and strictly
#' ordered within each variable. A double `@dimension` is rejected by S7's own
#' property check before the validator runs, so pass `6L` or an
#' `as.integer()`; the constructors take a plain `6` and convert it.
#'
#' @param basis_name A single string naming the family, printed by
#'   [print.basis()] and used by wrappers to build their own name. Not read by
#'   any computation.
#' @param dimension The number of functions in the basis, a single integer of
#'   at least 1. Must be of storage mode integer.
#' @param lower,upper The endpoints of the interval the basis lives on, or one
#'   endpoint per variable for a basis of several. Both finite, of the same
#'   length, and `lower[j] < upper[j]` for every `j`. Their length is what
#'   [basis_nvar()] reports.
#' @param basis_params A named list of whatever else the subclass needs: the
#'   knots and degree of a B-spline, the frequency of a Fourier basis, the
#'   marginal dimensions of a product. [print.basis()] shows it, abbreviating
#'   any numeric entry of more than four values. Defaults to an empty list.
#'
#' @return This class is abstract and cannot be constructed. A subclass
#'   constructed from it is an S7 object with properties `basis_name`
#'   (character), `dimension` (integer), `lower` and `upper` (numeric, one
#'   entry per variable) and `basis_params` (list).
#'
#' @seealso [bspline_basis()], [fourier_basis()] and [poly_basis()] for the
#'   three concrete families; [basis_eval()], [basis_deriv()], [basis_int()]
#'   and [basis_gram()] for the generics every basis answers; [check_basis()]
#'   to verify a subclass of your own.
#'
#' @examples
#' # The class is abstract, so it is a subclass that gets constructed.
#' b <- bspline_basis(dimension = 6)
#' b@dimension
#' c(b@lower, b@upper)
#' b@basis_params$degree
#'
#' # Its five properties are the same five on every basis in the package.
#' S7::prop_names(b)
#' S7::prop_names(fourier_basis(dimension = 5))
#'
#' # Inheriting from it and writing one method gives a working basis.
#' Bumps <- S7::new_class("Bumps", parent = basis)
#' S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
#'   exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension), "-")^2 / 0.12^2)
#' }
#' bump <- Bumps(basis_name = "bumps", dimension = 4L, lower = 0, upper = 1)
#'
#' # The three unwritten generics answer anyway, from finite differences.
#' basis_is_numerical(bump)
#' round(basis_deriv(bump, c(0.25, 0.75), order = 1), 4)
#'
#' # An interval the wrong way round is caught at construction.
#' try(Bumps(basis_name = "b", dimension = 4L, lower = 1, upper = 0))
#'
#' @export
basis <- S7::new_class(
  "basis",
  abstract = TRUE,
  properties = list(
    basis_name = S7::class_character,
    dimension = S7::class_integer,
    lower = S7::class_numeric,
    upper = S7::class_numeric,
    basis_params = S7::class_list
  ),
  validator = function(self) {
    if (length(self@basis_name) != 1L) {
      return("@basis_name must be a single string")
    }
    if (length(self@dimension) != 1L || is.na(self@dimension) ||
      self@dimension < 1L) {
      return("@dimension must be a single positive integer")
    }
    if (length(self@lower) < 1L || length(self@lower) != length(self@upper) ||
      !all(is.finite(self@lower)) || !all(is.finite(self@upper))) {
      return("@lower and @upper must be finite and of the same length")
    }
    if (any(self@lower >= self@upper)) {
      return("every @lower must be strictly less than its @upper")
    }
    NULL
  }
)


#' How Many Variables a Basis Takes
#'
#' @description
#' Returns the number of variables a basis is a function of: `1` for the three
#' shipped families and for any wrapper over them, and the number of margins
#' for a [tensor_basis()]. It is the number of columns [basis_eval()] expects
#' its `x` to have, so it answers what shape of input the object takes.
#'
#' @details
#' The count is `length(basis@lower)`, the endpoints carrying one entry per
#' variable. A basis therefore declares its input dimension by construction:
#' there is no separate property that could disagree with the interval, and
#' the validator already requires `@lower` and `@upper` to be of the same
#' length.
#'
#' At `basis_nvar(b) == 1` the evaluation points are a plain numeric vector
#' and the generics return an `n` by `@dimension` matrix. Above one they are a
#' matrix of `basis_nvar(b)` columns, one per variable; a plain vector is
#' reshaped by row, so `c(0.1, 0.2, 0.5, 0.6)` on a two-variable basis is the
#' two points `(0.1, 0.2)` and `(0.5, 0.6)`. `plot()` refuses a basis of more
#' than one variable, having no single picture to draw.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#'
#' @return A single positive integer, of length one and never `NA`.
#'
#' @seealso [basis_eval()], whose input shape this describes;
#'   [tensor_basis()], the only shipped family that answers more than `1`;
#'   [print.basis()], which shows it on the `Variables:` line.
#'
#' @examples
#' # One variable for every family and every wrapper over one.
#' basis_nvar(bspline_basis(dimension = 5))
#' basis_nvar(orthonorm_basis(bspline_basis(dimension = 5)))
#'
#' # A product answers with its number of margins.
#' tb <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#' basis_nvar(tb)
#'
#' # Which is the number of columns basis_eval() reads: four numbers are two
#' # points on a two-variable basis, one row each.
#' dim(basis_eval(tb, c(0.1, 0.2, 0.5, 0.6)))
#' dim(basis_eval(tb, cbind(c(0.1, 0.5), c(0.2, 0.6))))
#'
#' @export
basis_nvar <- function(basis) {
  length(basis@lower)
}


#' Validate the Arguments Every Basis Constructor Takes
#'
#' @description
#' Checks the three arguments common to [bspline_basis()], [fourier_basis()]
#' and [poly_basis()], and returns the dimension coerced to integer, which is
#' the storage mode the class validator requires. Every constructor calls it
#' first, before anything family-specific, so the three report the same errors
#' in the same words.
#'
#' @details
#' Three conditions, each throwing with `call. = FALSE`:
#'
#' - `lower` and `upper` must each be a single finite number, else
#'   `'lower' and 'upper' must be single finite numbers.`
#' - `lower` must be strictly below `upper`, else
#'   `'lower' must be strictly less than 'upper'.`
#' - `dimension` must be a single finite number, at least 1 and whole, else
#'   `'dimension' must be a single positive integer.`
#'
#' A whole-valued double passes and is converted, so a caller writing `6`
#' rather than `6L` is served. The multivariate case is not reachable here:
#' [tensor_basis()] builds its endpoints from margins already validated.
#'
#' @param lower,upper The endpoints of the interval, each a single finite
#'   number with `lower < upper`.
#' @param dimension The number of basis functions, a single whole number of at
#'   least 1, integer or double.
#'
#' @return `dimension` as an integer of length one.
#'
#' @seealso [check_eval_points()], the same role for the points a basis is
#'   evaluated at.
#'
#' @keywords internal
check_basis_args <- function(lower, upper, dimension) {
  if (!is.numeric(lower) || length(lower) != 1L || !is.finite(lower) ||
    !is.numeric(upper) || length(upper) != 1L || !is.finite(upper)) {
    stop("'lower' and 'upper' must be single finite numbers.", call. = FALSE)
  }
  if (lower >= upper) {
    stop("'lower' must be strictly less than 'upper'.", call. = FALSE)
  }
  if (!is.numeric(dimension) || length(dimension) != 1L ||
    !is.finite(dimension) || dimension < 1 || dimension != round(dimension)) {
    stop("'dimension' must be a single positive integer.", call. = FALSE)
  }
  as.integer(dimension)
}


#' Validate Evaluation Points Against a Basis
#'
#' @description
#' Checks that `x` is numeric and lies inside the basis interval, reshapes it
#' for a basis of several variables, and returns it with near-endpoint values
#' clamped exactly onto the endpoints. Every method that evaluates a basis
#' calls it first, so the three families and the fallbacks agree on what a
#' legal evaluation point is. Missing values pass through and become a missing
#' row of the result.
#'
#' @details
#' A basis is defined on its interval and nowhere else, so a point outside is
#' rejected with an error naming how many points were outside and what the
#' interval is. No extrapolation rule is applied: which one a model wants, if
#' any, belongs to the layer that owns the meaning of the covariate, and a
#' silent answer here would take that decision away from it.
#'
#' The comparison carries a tolerance of `1e-8` times the width of the
#' interval, so a point that is an endpoint up to rounding is accepted and
#' then set to the endpoint exactly. On \eqn{[0, 1]} that admits `1 + 5e-9`
#' and refuses `1 + 5e-8`; on \eqn{[0, 1000]} it admits `1000 + 5e-6`. The
#' clamp matters because a B-spline evaluated a rounding step past its last
#' knot is zero in every column.
#'
#' For a basis of several variables `x` is a matrix of one column per
#' variable, each column checked against its own endpoints. A plain vector is
#' accepted there and reshaped **by row**, so on two variables
#' `c(0.1, 0.2, 0.5, 0.6)` is the two points `(0.1, 0.2)` and `(0.5, 0.6)`. A
#' matrix of the wrong number of columns throws, naming the number wanted.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x A numeric vector of evaluation points for a basis of one variable,
#'   or a matrix of [basis_nvar()] columns for a basis of several. A
#'   one-column matrix is accepted for a basis of one variable and flattened;
#'   a vector is accepted for a basis of several and taken by row. `NA` is
#'   allowed anywhere and is neither range-checked nor clamped.
#'
#' @return `x` with near-endpoint values clamped onto the endpoints: a numeric
#'   vector for a basis of one variable, a numeric matrix of [basis_nvar()]
#'   columns otherwise.
#'
#' @seealso [clamp_to_range()], which does the per-variable work;
#'   [check_basis_args()], the same role for a constructor's arguments.
#'
#' @keywords internal
check_eval_points <- function(basis, x) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  d <- basis_nvar(basis)

  if (d == 1L) {
    if (is.matrix(x)) {
      if (ncol(x) != 1L) {
        stop("'x' must be a vector for a basis of one variable.", call. = FALSE)
      }
      x <- as.numeric(x)
    }
    return(clamp_to_range(x, basis@lower, basis@upper, "the basis interval"))
  }

  if (!is.matrix(x)) x <- matrix(x, ncol = d, byrow = TRUE)
  if (ncol(x) != d) {
    stop(sprintf(
      "'x' must have %d columns, one per variable of the basis.", d
    ), call. = FALSE)
  }
  for (j in seq_len(d)) {
    x[, j] <- clamp_to_range(
      x[, j], basis@lower[j], basis@upper[j],
      sprintf("the range of variable %d", j)
    )
  }
  x
}


#' Reject Points Outside a Range, and Clamp Those On Its Edge
#'
#' @description
#' Throws if any entry of `z` lies outside `[lo, hi]` by more than `1e-8`
#' times the width of the range, and otherwise returns `z` with anything
#' inside that tolerance moved onto the nearer endpoint. The one place the
#' range rule of the package is written; [check_eval_points()] calls it once
#' per variable.
#'
#' @details
#' The tolerance is relative to the width, `1e-8 * (hi - lo)`, so the rule
#' means the same on \eqn{[0, 1]} and on \eqn{[0, 1000]}. Entries that are
#' `NA` are excluded from both the test and the clamp, and travel through.
#'
#' The error message names the count and the range, as in
#' `2 of 5 evaluation points fall outside the basis interval [0, 1].`, and is
#' raised with `call. = FALSE`, so the phrase `what` supplies is what a user
#' sees in place of this function's name.
#'
#' @param z A numeric vector. `NA` entries are ignored.
#' @param lo,hi The endpoints, single finite numbers with `lo < hi`. Not
#'   validated here; the class validator has already required it.
#' @param what A phrase naming the range, spliced into the error message
#'   between "fall outside" and the interval. Callers pass
#'   `"the basis interval"` or `"the range of variable 2"`.
#'
#' @return `z`, of the same length, with entries within tolerance of an
#'   endpoint set to that endpoint exactly.
#'
#' @seealso [check_eval_points()], its only caller.
#'
#' @keywords internal
clamp_to_range <- function(z, lo, hi, what) {
  tol <- 1e-8 * (hi - lo)
  bad <- !is.na(z) & (z < lo - tol | z > hi + tol)
  if (any(bad)) {
    stop(sprintf(
      "%d of %d evaluation points fall outside %s [%s, %s].",
      sum(bad), length(z), what, format(lo), format(hi)
    ), call. = FALSE)
  }
  z[!is.na(z) & z < lo] <- lo
  z[!is.na(z) & z > hi] <- hi
  z
}


#' Name the Columns of a Basis Matrix
#'
#' @description
#' Sets `colnames(m)` to [basis_colnames()] of the basis and returns the
#' matrix. Every method returning an `n` by `@dimension` matrix ends with a
#' call to it, so the evaluation, the derivatives of every order, the integral
#' and the Gram matrix of one basis all carry the same column names in the
#' same order.
#'
#' @details
#' The names come from the basis, so a matrix with the wrong number of
#' columns produces R's own recycling error instead of a silently mislabeled
#' result. Nothing else is checked.
#'
#' @param m A numeric matrix with exactly `basis@dimension` columns.
#' @param basis A basis object, of any class inheriting from [basis].
#'
#' @return `m`, with its column names set and its contents untouched.
#'
#' @seealso [basis_colnames()], which supplies the names.
#'
#' @keywords internal
name_columns <- function(m, basis) {
  colnames(m) <- basis_colnames(basis)
  m
}
