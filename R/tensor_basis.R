#' @include numerical_fallbacks.R
NULL


#' Tensor Product Basis
#'
#' @description
#' The basis of products of the functions of several bases, one per variable:
#' \deqn{B(x_1, \ldots, x_D) = B_1(x_1) \otimes \cdots \otimes B_D(x_D).}
#' Constructed by \code{\link{tensor_basis}}.
#'
#' @details
#' Everything a tensor product needs follows from the marginals, because the
#' product separates. A partial derivative differentiates one marginal and
#' leaves the others alone; the integral over the box from its lower corner is
#' the product of the marginal integrals; and the Gram matrix is the Kronecker
#' product of the marginal Gram matrices, so a tensor of exactly integrated
#' marginals is exactly integrated too, however many variables it has.
#'
#' Columns follow the convention of \code{\link[base]{kronecker}}: the last
#' marginal varies fastest.
#'
#' @inheritParams basis
#' @param marginals The list of bases being multiplied, one per variable.
#'
#' @return An object of class \code{TensorBasis}.
#'
#' @seealso \code{\link{tensor_basis}}, \code{\link{basis_contract}}
#'
#' @examples
#' t2 <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
#' c(basis_nvar(t2), t2@dimension)
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
#' functions.
#'
#' @details
#' The result has \eqn{\prod_j K_j} functions and takes \eqn{D} variables, so
#' the evaluation points become a matrix with one column per variable. That
#' growth is the reason \code{\link{basis_contract}} exists: it computes what a
#' fit needs from the marginal evaluations alone, without ever forming the
#' product.
#'
#' Each marginal must take one variable. A product of products is flattened
#' rather than nested, so the marginals of the result are always the original
#' bases.
#'
#' @param ... The bases to multiply, or a single list of them.
#'
#' @return An object of class \code{\link{TensorBasis}}.
#'
#' @references
#' Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. \emph{Biometrics} 62, 1025-1036.
#'
#' @seealso \code{\link{basis_contract}}, \code{\link{basis_gram}}
#'
#' @examples
#' b <- tensor_basis(bspline_basis(dimension = 4), bspline_basis(dimension = 3, degree = 2))
#' b
#'
#' x <- cbind(c(0.2, 0.5), c(0.7, 0.1))
#' round(basis_eval(b, x), 4)
#'
#' # the Gram matrix is the Kronecker product of the marginal ones
#' g <- basis_gram(b)
#' max(abs(g - kronecker(
#'   basis_gram(bspline_basis(dimension = 4)),
#'   basis_gram(bspline_basis(dimension = 3, degree = 2))
#' )))
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
#' @description
#' The marginal names joined by dots, in the order the columns come out: the
#' last marginal varies fastest.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param ... Unused.
#' @return A character vector of length \code{basis@dimension}.
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
#' @description
#' The row-wise Kronecker product of the marginal evaluations.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param ... Unused.
#' @return A numeric matrix with \code{nrow(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_eval, TensorBasis) <- function(basis, x, ...) {
  name_columns(tensor_design(basis, x, order = NULL, integral = FALSE), basis)
}


#' Partial Derivatives of a Tensor Product Basis
#'
#' @name basis_deriv.TensorBasis
#' @description
#' The mixed partial derivative given by a multi-index, one order per variable.
#' @details
#' The product separates, so the derivative differentiates each marginal to its
#' own order and multiplies the results. An order beyond what a marginal
#' supports makes that factor, and so the whole product, zero.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param order An integer vector with one entry per variable.
#' @param ... Unused.
#' @return A numeric matrix with \code{nrow(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_deriv, TensorBasis) <- function(basis, x, order = 1L, ...) {
  name_columns(tensor_design(basis, x, order = order, integral = FALSE), basis)
}


#' Integral of a Tensor Product Basis
#'
#' @name basis_int.TensorBasis
#' @description
#' The integral over the box from the lower corner to each point.
#' @details
#' The integrand separates, so the multiple integral is the product of the
#' marginal integrals, and the anchor survives: a product in which every factor
#' is zero at the corner is zero at the corner.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param ... Unused.
#' @return A numeric matrix with \code{nrow(x)} rows and
#'   \code{basis@dimension} columns.
#' @keywords internal
S7::method(basis_int, TensorBasis) <- function(basis, x, ...) {
  name_columns(tensor_design(basis, x, order = NULL, integral = TRUE), basis)
}


#' Gram Matrix of a Tensor Product Basis
#'
#' @name basis_gram.TensorBasis
#' @description
#' The Kronecker product of the marginal Gram matrices.
#' @details
#' The integral over the box of a product of separable functions factorizes
#' into a product of one-dimensional integrals, so the matrix is separable and
#' costs one marginal Gram matrix per variable rather than one integration over
#' the box. A tensor of exactly integrated marginals is therefore exact at any
#' number of variables, where a quadrature over the box would not be.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param order An integer vector with one entry per variable.
#' @param at,weight Handled by the generic before dispatch; unused here.
#' @param ... Passed to the marginals.
#' @return A symmetric numeric matrix with \code{basis@dimension} rows and
#'   columns.
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
#' Evaluates each marginal, differentiated or integrated as asked, and
#' multiplies the results row by row.
#'
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param order An integer vector of derivative orders, or \code{NULL}.
#' @param integral Whether to integrate instead.
#'
#' @return A numeric matrix with \code{nrow(x)} rows.
#'
#' @keywords internal
tensor_design <- function(basis, x, order = NULL, integral = FALSE) {
  parts <- marginal_designs(basis, x, order, integral)
  Reduce(khatri_rao, parts)
}


#' Evaluate Every Marginal at Its Own Column
#'
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param order An integer vector of derivative orders, or \code{NULL}.
#' @param integral Whether to integrate instead.
#'
#' @return A list of numeric matrices, one per marginal.
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
#' For matrices with the same number of rows, the matrix whose \eqn{i}th row is
#' the Kronecker product of their \eqn{i}th rows, with the columns of the
#' second varying fastest.
#'
#' @param a,b Numeric matrices with the same number of rows.
#'
#' @return A numeric matrix with \code{ncol(a) * ncol(b)} columns.
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
#' Returns the values of the function the coefficients describe, without
#' necessarily forming the design matrix.
#'
#' @details
#' For an ordinary basis this is \code{basis_eval(basis, x) \%*\% coef} and
#' there is nothing to save. For a tensor product there is: the design matrix
#' has \eqn{\prod_j K_j} columns, so forming it is what makes a model with
#' several variables expensive, while the value it is used to compute needs
#' only the marginal evaluations.
#'
#' Coefficients come in two shapes.
#' \itemize{
#'   \item An \strong{array} of dimension \eqn{(K_1, \ldots, K_D)}, which is
#'     the general case. The rows are processed in blocks, so the peak memory
#'     is bounded by the block size rather than by the number of observations,
#'     however large the product.
#'   \item A \strong{list of factor matrices} \eqn{\Gamma_j} of size
#'     \eqn{K_j \times F}, the canonical polyadic form, in which the
#'     coefficient array is a sum of \eqn{F} outer products. Here the value is
#'     \eqn{\sum_f \prod_j B_j(x_j)^\top \gamma_{j,f}}, which costs
#'     \eqn{O(nF\sum_j K_j)} in both time and memory: neither the design matrix
#'     nor the coefficient array is ever formed.
#' }
#'
#' The second shape is what makes a model with high-order interactions
#' affordable, and what the factorized tensor product spline models of
#' Ruegamer (2024) estimate. Choosing the factors is a modeling decision and
#' belongs to the layer that owns the parameters; evaluating them is basis
#' arithmetic and belongs here.
#'
#' @param basis An object inheriting from class \code{basis}.
#' @param x Evaluation points: a numeric vector, or a matrix with one column
#'   per variable.
#' @param coef The coefficients: a vector, an array with one dimension per
#'   marginal, or a list of factor matrices.
#' @param ... Passed to methods.
#'
#' @return A numeric vector with one value per evaluation point, or a matrix
#'   with one column per column of \code{coef} when several sets are given.
#'
#' @references
#' Ruegamer, D. (2024). Scalable higher-order tensor product spline models.
#' \emph{Proceedings of AISTATS}.
#'
#' @seealso \code{\link{tensor_basis}}, \code{\link{basis_eval}}
#'
#' @examples
#' b <- tensor_basis(bspline_basis(dimension = 5), bspline_basis(dimension = 4))
#' x <- cbind(runif(10), runif(10))
#'
#' # a full array of coefficients
#' set.seed(1)
#' cf <- array(rnorm(20), dim = c(5, 4))
#' max(abs(basis_contract(b, x, cf) - basis_eval(b, x) %*% as.numeric(cf)))
#'
#' # or a rank-two factorization, which never forms either matrix
#' g <- list(matrix(rnorm(10), 5, 2), matrix(rnorm(8), 4, 2))
#' head(basis_contract(b, x, g))
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
#' @description
#' The design matrix times the coefficients. For a basis of one variable there
#' is nothing to avoid forming, so the definition is the computation.
#' @param basis An object inheriting from class \code{basis}.
#' @param x A numeric vector of evaluation points.
#' @param coef A numeric vector, or a matrix of several sets of coefficients.
#' @param ... Unused.
#' @return A numeric vector, or a matrix when \code{coef} is one.
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
#' @description
#' The value of the function the coefficients describe, computed from the
#' marginal evaluations without forming the tensor design matrix.
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param coef An array with one dimension per marginal, or a list of factor
#'   matrices in canonical polyadic form.
#' @param block The number of rows processed at once when \code{coef} is an
#'   array. It bounds the peak memory, which is otherwise what forming the
#'   design matrix would cost.
#' @param ... Unused.
#' @return A numeric vector with one value per row of \code{x}.
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
#' The canonical polyadic contraction: with \eqn{\Gamma_j} of size
#' \eqn{K_j \times F}, the value is
#' \eqn{\sum_f \prod_j B_j(x_j)^\top \gamma_{j,f}}.
#'
#' @details
#' Each marginal is evaluated once and multiplied by its own factor matrix,
#' giving \eqn{D} matrices of size \eqn{n \times F}; their elementwise product,
#' summed across columns, is the answer. Neither the design matrix nor the
#' coefficient array appears anywhere, so the cost is linear in the number of
#' variables where the array is exponential in it.
#'
#' @param basis A \code{\link{TensorBasis}} object.
#' @param x A numeric matrix with one column per variable.
#' @param coef A list of factor matrices, one per marginal, with a common
#'   number of columns.
#'
#' @return A numeric vector with one value per row of \code{x}.
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
