#' @include numerical_fallbacks.R
NULL


#' Tensor Product Basis
#'
#' @description
#' The S7 class of tensor product bases: all products of the functions of
#' several bases, one basis per variable,
#' \deqn{B(x_1, \ldots, x_D) = B_1(x_1) \otimes \cdots \otimes B_D(x_D).}
#' It carries the marginals whole, and every quantity it answers with is built
#' from theirs, so a product of exactly integrated marginals is itself exact
#' at any number of variables. Constructed by [tensor_basis()].
#'
#' @details
#' # Everything follows from the marginals
#'
#' The product separates, so each generic reduces to its marginals':
#'
#' - a partial derivative differentiates one marginal to its own order and
#'   leaves the others alone, `order` being a multi-index;
#' - the integral over the box from the lower corner is the product of the
#'   marginal integrals, and the anchor survives, a product whose every factor
#'   is zero at the corner being zero there;
#' - the Gram matrix is the Kronecker product of the marginal Gram matrices,
#'   which is a product of one-dimensional integrals and never a quadrature
#'   over the box.
#'
#' That last point is what keeps the construction affordable. A quadrature
#' over a box of \eqn{D} variables costs a node count exponential in \eqn{D}
#' and carries an error to match; a Kronecker product of exact marginal
#' matrices is exact, and costs one marginal Gram matrix per variable.
#'
#' # Column order
#'
#' Columns follow [base::kronecker()]: the **last** marginal varies fastest.
#' At two variables of 3 and 2 functions the names are `bs1.P0`, `bs1.P1`,
#' `bs2.P0`, `bs2.P1`, `bs3.P0`, `bs3.P1`.
#'
#' That order matters when coefficients are supplied as an array, R storing an
#' array with its *first* index fastest. [basis_contract()] reverses the
#' dimensions for you; a hand-written `as.numeric(coef)` does not, and pairs
#' every coefficient with the wrong function.
#'
#' # What the validator enforces
#'
#' `@marginals` must be non-empty and hold only bases, and `@dimension` must
#' equal the product of their dimensions.
#'
#' @inheritParams basis
#' @param marginals The list of bases being multiplied, one per variable, each
#'   of one variable itself. Kept whole, so each can still be evaluated and
#'   asked what it is.
#'
#' @return An object of class `TensorBasis`, inheriting from [basis], with the
#'   five properties of a basis plus `marginals`, and `basis_params` holding
#'   `marginal_dimensions`.
#'
#' @seealso [tensor_basis()], the constructor; [basis_contract()], which
#'   evaluates a fit without forming the product; [basis_gram()] for the
#'   Kronecker identity.
#'
#' @examples
#' t2 <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#' c(basis_nvar(t2), t2@dimension)
#'
#' # The last marginal varies fastest, and the names record both margins.
#' basis_colnames(tensor_basis(bspline_basis(dimension = 3, degree = 2),
#'                             poly_basis(dimension = 2)))
#'
#' # A product is exact wherever its margins are.
#' basis_is_numerical(t2)
#'
#' @export
TensorBasis <- S7::new_class(
  "TensorBasis",
  parent = basis,
  properties = list(marginals = S7::class_list),
  validator = function(self) {
    if (!length(self@marginals)) return("@marginals must not be empty")
    if (!all(vapply(self@marginals, function(m) S7::S7_inherits(m, basis),
      logical(1)
    ))) {
      return("every element of @marginals must be a basis")
    }
    if (self@dimension != prod(vapply(self@marginals, function(m) m@dimension,
      numeric(1)
    ))) {
      return("@dimension must be the product of the marginal dimensions")
    }
    NULL
  }
)


#' Construct a Tensor Product Basis
#'
#' @description
#' Multiplies bases, one per variable, into the basis of all products of their
#' functions. This is how a smooth surface of several covariates is built from
#' one-dimensional pieces: the result is a basis like any other, evaluated at
#' a matrix of points instead of a vector, and every generic answers exactly
#' where the margins do.
#'
#' @details
#' # Size
#'
#' The result has \eqn{\prod_j K_j} functions and takes \eqn{D} variables, so
#' the evaluation points become a matrix of \eqn{D} columns. The dimension
#' grows geometrically: four cubic B-splines of eight functions each give
#' 4096 columns, and the design matrix at 20000 observations would be 625 MB.
#'
#' [basis_contract()] exists for that reason. It computes what a fit needs
#' from the marginal evaluations, in blocks for a full coefficient array and
#' without forming anything at all for a factorized one.
#'
#' # Flattening
#'
#' Each marginal must take one variable. A product of products is flattened
#' and never nested, so `@marginals` always holds the original bases and
#' [basis_nvar()] counts them all: multiplying a two-way product by a third
#' basis gives one three-way product with three margins.
#'
#' # The interval
#'
#' `@lower` and `@upper` hold one entry per variable, taken from the margins,
#' so the domain is the box they span and a point outside any margin's
#' interval throws.
#'
#' @param ... The bases to multiply, two or more, or a single list of them.
#'   Each must take one variable; a [TensorBasis] among them is flattened into
#'   its own margins.
#'
#' @return An object of class [TensorBasis] with `basis_name`
#'   `tensor(<names>)`, `basis_params` holding `marginal_dimensions`, and
#'   column names pasting the margins' with dots.
#'
#' @references
#' Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. *Biometrics* **62**, 1025-1036.
#'
#' @seealso [basis_contract()], for evaluating a fit without the design
#'   matrix; [basis_gram()], which is a Kronecker product here; [dr_basis()]
#'   and [constrain_basis()], which apply to a product as to any basis.
#'
#' @examples
#' b <- tensor_basis(bspline_basis(dimension = 4),
#'                   bspline_basis(dimension = 3, degree = 2))
#' b
#'
#' x <- cbind(c(0.2, 0.5), c(0.7, 0.1))
#' round(basis_eval(b, x), 4)
#'
#' # The Gram matrix is the Kronecker product of the marginal ones, exactly.
#' max(abs(basis_gram(b) - kronecker(
#'   basis_gram(bspline_basis(dimension = 4)),
#'   basis_gram(bspline_basis(dimension = 3, degree = 2))
#' )))
#'
#' # The integral is anchored at the lower corner of the box.
#' basis_int(b, cbind(0, 0))
#'
#' # A product of products is flattened, so the margins stay one-dimensional.
#' t3 <- tensor_basis(b, fourier_basis(dimension = 3))
#' c(basis_nvar(t3), length(t3@marginals), t3@dimension)
#'
#' # The dimension grows geometrically. basis_contract() exists for that.
#' vapply(2:5, function(D) {
#'   tensor_basis(rep(list(bspline_basis(dimension = 8)), D))@dimension
#' }, numeric(1))
#'
#' @export
tensor_basis <- function(...) {
  args <- list(...)
  if (length(args) == 1L && is.list(args[[1L]]) &&
    !S7::S7_inherits(args[[1L]], basis)) {
    args <- args[[1L]]
  }
  if (!length(args)) {
    stop("'tensor_basis' needs at least one basis.", call. = FALSE)
  }

  # A product of products is one product: flattening keeps the marginals the
  # bases the caller named, and keeps the column order unambiguous.
  marginals <- unlist(lapply(args, function(b) {
    if (S7::S7_inherits(b, TensorBasis)) b@marginals else list(b)
  }), recursive = FALSE)

  for (m in marginals) {
    if (!S7::S7_inherits(m, basis)) {
      stop("every argument must be a basis.", call. = FALSE)
    }
    if (basis_nvar(m) != 1L) {
      stop("every marginal must take a single variable.", call. = FALSE)
    }
  }
  if (length(marginals) == 1L) return(marginals[[1L]])

  TensorBasis(
    basis_name = paste0(
      "tensor(",
      paste(vapply(marginals, function(m) m@basis_name, character(1)),
        collapse = ", "
      ),
      ")"
    ),
    dimension = as.integer(prod(vapply(marginals, function(m) m@dimension,
      numeric(1)
    ))),
    lower = vapply(marginals, function(m) m@lower, numeric(1)),
    upper = vapply(marginals, function(m) m@upper, numeric(1)),
    basis_params = list(
      marginal_dimensions = vapply(marginals, function(m) m@dimension,
        integer(1)
      )
    ),
    marginals = marginals
  )
}


#' Column Names of a Tensor Product Basis
#'
#' @name basis_colnames.TensorBasis
#'
#' @description
#' Pastes the margins' column names with dots, in the order the columns come
#' out: the last marginal varies fastest, following [base::kronecker()]. A
#' column named `bs2.cos1` is the product of the parent B-spline's second
#' function with the Fourier margin's first cosine, so a coefficient's name
#' says which marginal function it belongs to in each variable.
#'
#' @details
#' The names are built with the same recycling [khatri_rao()] uses, the second
#' factor repeated within each element of the first. `outer()` would give the
#' transpose of this, which reads plausibly and labels every column but the
#' first and last wrongly.
#'
#' @param basis A [TensorBasis] object.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A character vector of length `basis@dimension`.
#'
#' @seealso [basis_colnames()] for the generic; [khatri_rao()], whose ordering
#'   this matches.
#'
#' @keywords internal
S7::method(basis_colnames, TensorBasis) <- function(basis, ...) {
  parts <- lapply(basis@marginals, basis_colnames)
  # The names are built with the same recycling as khatri_rao(): the second
  # factor repeated within each element of the first. outer() would give the
  # transpose of this, which reads plausibly and labels every column but the
  # first and last wrongly.
  Reduce(function(a, b) {
    paste(rep(a, each = length(b)), rep(b, times = length(a)), sep = ".")
  }, parts)
}


#' Evaluate a Tensor Product Basis
#'
#' @name basis_eval.TensorBasis
#'
#' @description
#' Evaluates each margin at its own column of the points and takes the
#' row-wise Kronecker product, so row \eqn{i} of the result is
#' \eqn{B_1(x_{i1}) \otimes \cdots \otimes B_D(x_{iD})}. The matrix has
#' \eqn{\prod_j K_j} columns, so at several variables it is worth avoiding,
#' and [basis_contract()] avoids it.
#'
#' @details
#' Cost is one marginal evaluation per variable plus the products, and the
#' result is `nrow(x)` by `basis@dimension`, which at four margins of eight
#' functions and 20000 points is 625 MB. Nothing is cached.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable, or a vector taken
#'   by row. Each column is checked against its own margin's interval.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `nrow(x)` rows and `basis@dimension` columns,
#'   with column names pasting the margins' with dots.
#'
#' @seealso [tensor_design()] and [khatri_rao()], which do the work;
#'   [basis_contract()], for the value of a fit without this matrix.
#'
#' @keywords internal
S7::method(basis_eval, TensorBasis) <- function(basis, x, ...) {
  name_columns(tensor_design(basis, x, order = NULL, integral = FALSE), basis)
}


#' Partial Derivatives of a Tensor Product Basis
#'
#' @name basis_deriv.TensorBasis
#'
#' @description
#' Returns the mixed partial derivative named by a multi-index, one order per
#' variable: `c(2, 0)` is \eqn{\partial^2/\partial x_1^2} and `c(1, 1)` is
#' \eqn{\partial^2/\partial x_1 \partial x_2}. Exact wherever the margins are,
#' and the one place in the package a mixed partial is available at all, the
#' numerical fallback refusing them.
#'
#' @details
#' The product separates, so the derivative differentiates each margin to its
#' own order and multiplies the results. No cross term appears and no stencil
#' in the plane is needed, which is why the accuracy of a mixed partial here
#' is the accuracy of the margins' own derivatives.
#'
#' An order beyond what a margin carries makes that factor zero, so the whole
#' product is zero: `c(4, 0)` on a cubic B-spline margin is the exact zero
#' matrix whatever the other margins do.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable.
#' @param order An integer vector with one entry per variable, or a single
#'   `0`. A single non-zero order throws, having two readings; see
#'   [check_order()].
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `nrow(x)` rows and `basis@dimension` columns.
#'
#' @seealso [tensor_design()], which does the work; [basis_deriv()] for the
#'   generic and the multi-index rule.
#'
#' @keywords internal
S7::method(basis_deriv, TensorBasis) <- function(basis, x, order = 1L, ...) {
  name_columns(tensor_design(basis, x, order = order, integral = FALSE), basis)
}


#' Integral of a Tensor Product Basis
#'
#' @name basis_int.TensorBasis
#'
#' @description
#' Returns the integral over the box from the lower corner to each point, one
#' iterated integral per variable, in closed form wherever the margins have
#' one. It is the only integral over a box in the package, the numerical
#' fallback taking one variable alone.
#'
#' @details
#' The integrand separates, so the multiple integral is the product of the
#' marginal integrals. The anchoring convention of [basis_int()] survives
#' without a correction: a product in which every factor is zero at the corner
#' is zero at the corner, so the row at `basis@lower` is exactly zero.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable, the upper corners
#'   of the boxes.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric matrix with `nrow(x)` rows and `basis@dimension` columns,
#'   exactly zero in the row at the lower corner.
#'
#' @seealso [tensor_design()], which does the work; [basis_int()] for the
#'   generic and the anchoring convention.
#'
#' @keywords internal
S7::method(basis_int, TensorBasis) <- function(basis, x, ...) {
  name_columns(tensor_design(basis, x, order = NULL, integral = TRUE), basis)
}


#' Gram Matrix of a Tensor Product Basis
#'
#' @name basis_gram.TensorBasis
#'
#' @description
#' Returns the Kronecker product of the marginal Gram matrices, which is the
#' integral over the box exactly, with no quadrature. A product of exactly
#' integrated margins is therefore exact at any number of variables, where a
#' rule over the box would cost a node count exponential in the number of
#' variables and still carry an error.
#'
#' @details
#' The integral over the box of a product of separable functions factorizes
#' into one-dimensional integrals, so
#' \eqn{G = G_1 \otimes \cdots \otimes G_D} with \eqn{G_j} the margin's own
#' matrix at that margin's own order. It costs one marginal Gram matrix per
#' variable and no integration over the box at all.
#'
#' The result is symmetrized as `(G + t(G))/2`, a Kronecker product of
#' symmetric matrices being symmetric only up to the order its entries were
#' formed in. It is singular whenever any margin's is, so an `order` with any
#' non-zero entry gives a singular matrix.
#'
#' @param basis A [TensorBasis] object.
#' @param order The derivative order, an integer vector with one entry per
#'   variable, or a single `0`. Each margin is asked for its own entry.
#' @param at,weight Handled in the body of [basis_gram()] before dispatch, so
#'   they never arrive here. Named only because S7 requires a method's formals
#'   to contain the generic's.
#' @param ... Passed to each margin's method.
#'
#' @return A symmetric numeric matrix of `basis@dimension` rows and columns,
#'   with the product column names on both margins.
#'
#' @seealso [basis_gram()] for the generic and the alternative measures;
#'   [base::kronecker()], whose column order this follows.
#'
#' @keywords internal
S7::method(basis_gram, TensorBasis) <- function(basis, order = 0L, at = NULL,
                                                weight = NULL, ...) {
  order <- check_order(order, basis_nvar(basis))
  g <- Reduce(kronecker, Map(
    function(m, d) basis_gram(m, order = d, ...), basis@marginals, order
  ))
  g <- (g + t(g)) / 2
  nm <- basis_colnames(basis)
  dimnames(g) <- list(nm, nm)
  g
}


#' The Row-Wise Kronecker Product of the Marginal Designs
#'
#' @description
#' Evaluates each margin, differentiated or integrated as asked, and reduces
#' the results with [khatri_rao()] into the product's own design matrix. The
#' one body behind [basis_eval.TensorBasis()], [basis_deriv.TensorBasis()]
#' and [basis_int.TensorBasis()], which differ only in what they ask the
#' margins for.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable, already validated.
#' @param order An integer vector with one entry per variable, or `NULL` for
#'   the functions themselves.
#' @param integral `TRUE` to ask each margin for its anchored integral instead
#'   of its evaluation. Not combined with `order`; each caller sets at most
#'   one.
#'
#' @return A numeric matrix with `nrow(x)` rows and `basis@dimension` columns,
#'   no dimnames; callers add them through [name_columns()].
#'
#' @seealso [marginal_designs()] and [khatri_rao()], its two steps.
#'
#' @keywords internal
tensor_design <- function(basis, x, order = NULL, integral = FALSE) {
  parts <- marginal_designs(basis, x, order, integral)
  Reduce(khatri_rao, parts)
}


#' Evaluate Every Marginal at Its Own Column
#'
#' @description
#' Returns a list of the margins' design matrices, each evaluated at its own
#' column of the points, and each differentiated or integrated as asked. The
#' first step of [tensor_design()].
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable.
#' @param order An integer vector with one entry per variable, or `NULL` for
#'   the functions themselves.
#' @param integral `TRUE` to ask each margin for its anchored integral.
#'
#' @return A list of `basis_nvar(basis)` numeric matrices, the \eqn{j}th with
#'   `nrow(x)` rows and `basis@marginals[[j]]@dimension` columns.
#'
#' @seealso [tensor_design()], its only caller.
#'
#' @keywords internal
marginal_designs <- function(basis, x, order = NULL, integral = FALSE) {
  lapply(seq_along(basis@marginals), function(j) {
    m <- basis@marginals[[j]]
    z <- x[, j]
    if (integral) {
      basis_int(m, z)
    } else if (is.null(order)) {
      basis_eval(m, z)
    } else {
      basis_deriv(m, z, order = order[j])
    }
  })
}


#' Row-Wise Kronecker Product of Two Matrices
#'
#' @description
#' For two matrices with the same number of rows, returns the matrix whose
#' \eqn{i}th row is the Kronecker product of their \eqn{i}th rows. Reducing
#' the margins' design matrices with it gives the product's, one row of
#' \eqn{\prod_j K_j} entries per observation.
#'
#' @details
#' The columns come out with `b` varying fastest, matching
#' [base::kronecker()] and [basis_colnames.TensorBasis()]. That is the whole
#' column-order convention of a tensor basis, and the reason [basis_contract()]
#' reverses an array's dimensions before flattening it.
#'
#' The result has `ncol(a) * ncol(b)` columns, so reducing across several
#' margins grows geometrically; nothing here bounds that, and
#' [basis_contract()] is what avoids paying it.
#'
#' @param a,b Numeric matrices with the same number of rows. The row counts
#'   are not checked; a mismatch gives R's own recycling behavior.
#'
#' @return A numeric matrix with `nrow(a)` rows and `ncol(a) * ncol(b)`
#'   columns, no dimnames.
#'
#' @seealso [tensor_design()], its only caller;
#'   [basis_colnames.TensorBasis()], which names its columns.
#'
#' @keywords internal
khatri_rao <- function(a, b) {
  pa <- ncol(a)
  pb <- ncol(b)
  a[, rep(seq_len(pa), each = pb), drop = FALSE] *
    b[, rep(seq_len(pb), times = pa), drop = FALSE]
}


#' Evaluate a Basis Against Coefficients
#'
#' @description
#' Returns the values of the function a coefficient vector describes,
#' \eqn{B(x)\beta}, computed without forming \eqn{B(x)} in full where that is
#' possible. For an ordinary basis it is the design matrix times the
#' coefficients and there is nothing to save; for a [tensor_basis()] the
#' design matrix has \eqn{\prod_j K_j} columns, and avoiding it is what keeps
#' a model of several variables affordable.
#'
#' @details
#' # Two shapes of coefficient
#'
#' An **array** of dimension \eqn{(K_1, \ldots, K_D)} is the general case.
#' The rows are processed in blocks of `block`, so the peak memory is the
#' block's design in place of the whole one: at four margins of eight functions
#' and 20000 points, 32 MB against 625 MB, measured, at the same speed
#' (0.83 s against 1.00 s) and to an exact agreement.
#'
#' A **list of factor matrices** \eqn{\Gamma_j} of size \eqn{K_j \times F} is
#' the canonical polyadic form, in which the coefficient array is a sum of
#' \eqn{F} outer products. The value is then
#' \eqn{\sum_f \prod_j B_j(x_j)^\top \gamma_{j,f}}, costing
#' \eqn{O(nF\sum_j K_j)} in both time and memory: neither the design matrix
#' nor the coefficient array is formed anywhere. On the same problem at
#' \eqn{F = 3} it is 0.03 s against 0.83.
#'
#' That second shape is what brings a model with high-order interactions
#' within reach, and what the factorized tensor product spline models of
#' Ruegamer (2024) estimate. Choosing the factors is a modeling decision and
#' belongs to the layer that owns the parameters; evaluating them is basis
#' arithmetic and belongs here.
#'
#' # An array and a vector are read differently
#'
#' A tensor design's columns run with the **last** margin fastest, following
#' [base::kronecker()], and an R array is stored with its **first** index
#' fastest. An array is therefore transposed with `aperm()` before it is
#' flattened, so that `coef[m1, ..., mD]` is the coefficient of the column
#' named for those marginal functions.
#'
#' A plain **vector** has no dimensions to reverse and is taken in the
#' design's own column order. So a vector and an array holding the same
#' numbers in the same storage order describe different functions, and the
#' identity to check against is
#' `basis_eval(b, x) %*% as.numeric(aperm(coef))`. Flattening without the
#' `aperm()` pairs every coefficient with the wrong function and returns a
#' perfectly finite number: 1.16 on the example below.
#'
#' @param basis A basis object, of any class inheriting from [basis].
#' @param x Evaluation points: a numeric vector for a basis of one variable,
#'   or a matrix with one column per variable.
#' @param coef The coefficients. For an ordinary basis, a numeric vector of
#'   length `basis@dimension`, or a matrix of several such columns. For a
#'   [TensorBasis], an array of dimension `(K_1, ..., K_D)`, a vector of
#'   `basis@dimension` values in the design's own column order, or a list of
#'   `D` factor matrices each with `F` columns. A wrong length or a wrong set
#'   of dimensions throws, naming what was wanted.
#' @param ... Passed to methods. The [TensorBasis] method reads `block` from
#'   it, the number of rows processed at once, default `1024`.
#'
#' @return A numeric vector with one value per evaluation point, or a matrix
#'   with one column per column of `coef` when several sets are given to an
#'   ordinary basis.
#'
#' @references
#' Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
#' *Proceedings of AISTATS*.
#'
#' @seealso [tensor_basis()], the case this exists for; [basis_eval()], the
#'   design matrix it avoids forming.
#'
#' @examples
#' b <- tensor_basis(bspline_basis(dimension = 5), bspline_basis(dimension = 4))
#' set.seed(1)
#' x <- cbind(runif(10), runif(10))
#'
#' # A full array of coefficients. The identity holds against the flattening
#' # that matches the column order, which reverses the array's dimensions.
#' cf <- array(rnorm(20), dim = c(5, 4))
#' max(abs(basis_contract(b, x, cf) -
#'         basis_eval(b, x) %*% as.numeric(aperm(cf))))
#'
#' # Flattening without that gives a plausible wrong answer, not an error.
#' max(abs(basis_contract(b, x, cf) - basis_eval(b, x) %*% as.numeric(cf)))
#'
#' # A rank-two factorization, which forms neither the design nor the array.
#' g <- list(matrix(rnorm(10), 5, 2), matrix(rnorm(8), 4, 2))
#' head(basis_contract(b, x, g))
#'
#' # And agrees with the array it stands for.
#' full <- outer(g[[1]][, 1], g[[2]][, 1]) + outer(g[[1]][, 2], g[[2]][, 2])
#' max(abs(basis_contract(b, x, g) - basis_contract(b, x, full)))
#'
#' # On an ordinary basis it is the design matrix times the coefficients.
#' p <- poly_basis(dimension = 4)
#' max(abs(basis_contract(p, c(0.2, 0.8), 1:4) -
#'         basis_eval(p, c(0.2, 0.8)) %*% 1:4))
#'
#' @export
basis_contract <- S7::new_generic(
  "basis_contract", "basis",
  function(basis, x, coef, ...) {
    x <- check_eval_points(basis, x)
    S7::S7_dispatch()
  }
)


#' Contract an Ordinary Basis Against Coefficients
#'
#' @name basis_contract.basis
#'
#' @description
#' Returns `basis_eval(basis, x) %*% coef`. For a basis of one variable the
#' design matrix has `basis@dimension` columns and there is nothing worth
#' avoiding, so the definition is the computation. The method every class
#' inherits except [TensorBasis].
#'
#' @details
#' `coef` is taken in the design's own column order, there being one obvious
#' one. Several sets of coefficients may be given as the columns of a matrix,
#' and the result then has one column each.
#'
#' @param basis A basis object, of any class inheriting from [basis] other
#'   than [TensorBasis].
#' @param x A numeric vector of evaluation points inside the basis interval.
#' @param coef A numeric vector of length `basis@dimension`, or a matrix with
#'   that many rows and one column per set of coefficients.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector with one value per point, or a matrix with one
#'   column per column of `coef`.
#'
#' @seealso [basis_contract()] for the generic;
#'   [basis_contract.TensorBasis()] for the case where the design matrix is
#'   worth avoiding.
#'
#' @keywords internal
S7::method(basis_contract, basis) <- function(basis, x, coef, ...) {
  cf <- as.matrix(coef)
  if (nrow(cf) != basis@dimension) {
    stop(sprintf(
      "'coef' must have %d entries, one per basis function.", basis@dimension
    ), call. = FALSE)
  }
  out <- basis_eval(basis, x) %*% cf
  if (ncol(cf) == 1L) drop(out) else out
}


#' Contract a Tensor Product Basis Against Coefficients
#'
#' @name basis_contract.TensorBasis
#'
#' @description
#' Returns the value of the function the coefficients describe, computed from
#' the marginal evaluations. A list of factor matrices goes to
#' [contract_cp()], which forms nothing; an array or a vector is contracted in
#' blocks of `block` rows, so the peak memory is one block's design matrix and
#' not the whole one.
#'
#' @details
#' # The blocked route
#'
#' Each block of rows is evaluated through [tensor_design()] and multiplied by
#' the flattened coefficients, and the block is then discarded. Measured at
#' four margins of eight functions and 20000 points, where the full design is
#' 625 MB: 32 MB at the default block and 0.83 s, against 625 MB and 1.00 s
#' for the design route, agreeing exactly. Raising `block` past a few thousand
#' buys no speed and costs memory linearly.
#'
#' # The flattening
#'
#' An array is passed through `aperm()` before flattening, because the design's
#' columns run with the last margin fastest while an R array is stored with
#' its first index fastest. A plain vector has no dimensions to reverse and is
#' taken in the design's column order as it stands, so the two shapes describe
#' different functions from the same numbers. See [basis_contract()].
#'
#' A `coef` of the wrong length, or an array whose dimensions are not the
#' margins', throws with the expected values named.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable.
#' @param coef An array of dimension `(K_1, ..., K_D)`, a numeric vector of
#'   `basis@dimension` values in the design's own column order, or a list of
#'   `D` factor matrices in canonical polyadic form.
#' @param block The number of rows processed at once when `coef` is an array
#'   or a vector, default `1024`. It bounds the peak memory, which is
#'   otherwise what forming the design matrix would cost. Ignored for a list.
#' @param ... Unused, and accepted so that the signature matches the generic's.
#'
#' @return A numeric vector with one value per row of `x`.
#'
#' @seealso [contract_cp()] for the factorized route; [basis_contract()] for
#'   the generic and the two shapes of coefficient.
#'
#' @keywords internal
S7::method(basis_contract, TensorBasis) <- function(basis, x, coef,
                                                    block = 1024L, ...) {
  dims <- vapply(basis@marginals, function(m) m@dimension, integer(1))

  if (is.list(coef)) return(contract_cp(basis, x, coef))

  if (length(coef) != basis@dimension) {
    stop(sprintf(
      "'coef' must hold %d values, one per basis function.", basis@dimension
    ), call. = FALSE)
  }
  if (!is.null(dim(coef)) && !identical(as.integer(dim(coef)), dims)) {
    stop(sprintf(
      "'coef' must have dimensions %s.", paste(dims, collapse = " by ")
    ), call. = FALSE)
  }

  # The columns run with the LAST marginal fastest, the convention of
  # kronecker(); an R array is stored with the first index fastest. Reversing
  # the dimensions lines the two up, so that coef[m1, ..., mD] is the
  # coefficient of the column named for those marginal functions. Flattening
  # without it would pair every coefficient with the wrong function, and the
  # answer would be a perfectly finite number.
  cf <- if (is.null(dim(coef)) || length(dims) == 1L) {
    as.numeric(coef)
  } else {
    as.numeric(aperm(coef))
  }

  n <- nrow(x)
  out <- numeric(n)
  starts <- seq.int(1L, max(n, 1L), by = block)
  for (s in starts) {
    idx <- seq.int(s, min(s + block - 1L, n))
    out[idx] <- drop(tensor_design(basis, x[idx, , drop = FALSE]) %*% cf)
  }
  out
}


#' Contract a Tensor Product Basis Against Factor Matrices
#'
#' @description
#' The canonical polyadic contraction. With \eqn{\Gamma_j} of size
#' \eqn{K_j \times F}, the coefficient array is a sum of \eqn{F} outer
#' products and the value is
#' \eqn{\sum_f \prod_j B_j(x_j)^\top \gamma_{j,f}}. Neither the design matrix
#' nor the coefficient array appears anywhere.
#'
#' @details
#' Each margin is evaluated once and multiplied by its own factor matrix,
#' giving \eqn{D} matrices of size \eqn{n \times F}; their elementwise
#' product, summed across columns, is the answer. The cost is
#' \eqn{O(nF\sum_j K_j)} in time and memory, linear in the number of variables
#' where the array is exponential in it: measured at four margins of eight
#' functions, 20000 points and \eqn{F = 3}, 0.03 s against the blocked array
#' route's 0.83 s.
#'
#' Against the array it stands for, the two routes agree to 1.7e-16.
#'
#' @param basis A [TensorBasis] object.
#' @param x A numeric matrix with one column per variable.
#' @param coef A list of `basis_nvar(basis)` numeric matrices, the \eqn{j}th
#'   with `basis@marginals[[j]]@dimension` rows, all with the same number of
#'   columns \eqn{F}.
#'
#' @return A numeric vector with one value per row of `x`.
#'
#' @references
#' Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
#' *Proceedings of AISTATS*.
#'
#' @seealso [basis_contract.TensorBasis()], its only caller.
#'
#' @keywords internal
contract_cp <- function(basis, x, coef) {
  dims <- vapply(basis@marginals, function(m) m@dimension, integer(1))
  if (length(coef) != length(dims)) {
    stop(sprintf(
      "'coef' must hold one factor matrix per marginal (%d).", length(dims)
    ), call. = FALSE)
  }
  coef <- lapply(coef, as.matrix)
  ranks <- vapply(coef, ncol, integer(1))
  if (length(unique(ranks)) != 1L) {
    stop("every factor matrix must have the same number of columns.",
      call. = FALSE
    )
  }
  rows <- vapply(coef, nrow, integer(1))
  if (!identical(rows, dims)) {
    stop(sprintf(
      "factor matrix %d has %d rows, but its marginal has %d functions.",
      which(rows != dims)[1L], rows[which(rows != dims)[1L]],
      dims[which(rows != dims)[1L]]
    ), call. = FALSE)
  }

  phi <- lapply(seq_along(coef), function(j) {
    basis_eval(basis@marginals[[j]], x[, j]) %*% coef[[j]]
  })
  rowSums(Reduce(`*`, phi))
}
