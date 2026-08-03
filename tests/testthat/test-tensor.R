# Tensor products: everything follows from the marginals, so everything is
# checked against them, computed independently.

tb2 <- function() {
  tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
}

kron_rows <- function(mats) {
  n <- nrow(mats[[1L]])
  t(vapply(seq_len(n), function(i) {
    Reduce(kronecker, lapply(mats, function(m) m[i, ]))
  }, numeric(prod(vapply(mats, ncol, integer(1))))))
}


test_that("a product declares its variables, dimension and domain", {
  b <- tb2()
  expect_identical(basis_nvar(b), 2L)
  expect_identical(b@dimension, 12L)
  expect_identical(b@lower, c(0, 0))
  expect_identical(basis_nvar(bspline_basis(dimension = 5)), 1L)

  # a shifted box keeps each marginal's own range
  s <- tensor_basis(
    bspline_basis(-2, 3, dimension = 4),
    fourier_basis(1, 7, dimension = 3)
  )
  expect_identical(s@lower, c(-2, 1))
  expect_identical(s@upper, c(3, 7))
})


test_that("the columns are ordered and named with the LAST marginal fastest", {
  # The convention of kronecker(). A name in the wrong order labels every
  # column but the first and last incorrectly, and nothing else would notice.
  b <- tb2()
  expect_identical(
    basis_colnames(b),
    paste(rep(basis_colnames(bspline_basis(dimension = 4)), each = 3),
      rep(basis_colnames(fourier_basis(dimension = 3)), times = 4),
      sep = "."
    )
  )
  expect_identical(basis_colnames(b)[2L], "bs1.sin1")
})


test_that("evaluation is the row-wise Kronecker product of the marginals", {
  set.seed(41)
  b <- tb2()
  x <- cbind(runif(9), runif(9))
  want <- kron_rows(list(
    basis_eval(bspline_basis(dimension = 4), x[, 1]),
    basis_eval(fourier_basis(dimension = 3), x[, 2])
  ))
  expect_equal(basis_eval(b, x), want, ignore_attr = TRUE)
})


test_that("a partial derivative differentiates one marginal and leaves the rest", {
  set.seed(42)
  b <- tb2()
  m1 <- bspline_basis(dimension = 4)
  m2 <- fourier_basis(dimension = 3)
  x <- cbind(runif(7), runif(7))

  for (o in list(c(1, 0), c(0, 1), c(1, 1), c(2, 1), c(3, 4))) {
    want <- kron_rows(list(
      basis_deriv(m1, x[, 1], order = o[1]),
      basis_deriv(m2, x[, 2], order = o[2])
    ))
    expect_equal(basis_deriv(b, x, order = o), want,
      ignore_attr = TRUE, label = paste(o, collapse = ",")
    )
  }
  # order 0 in every coordinate is the evaluation
  expect_equal(basis_deriv(b, x, order = c(0, 0)), basis_eval(b, x))
  expect_equal(basis_deriv(b, x, order = 0), basis_eval(b, x))
})


test_that("a single non-zero order is refused as ambiguous", {
  # It could mean that order in every coordinate, or that total order, and
  # guessing between them would be a silent choice.
  b <- tb2()
  x <- cbind(0.3, 0.4)
  expect_error(basis_deriv(b, x, order = 1), "one entry per variable")
  expect_error(basis_deriv(b, x, order = c(1, 2, 3)), "one entry per variable")
  expect_silent(basis_deriv(b, x, order = c(1, 2)))
})


test_that("the integral is anchored at the lower CORNER and separates", {
  set.seed(43)
  b <- tensor_basis(
    bspline_basis(-1, 2, dimension = 4),
    bspline_basis(3, 5, dimension = 3, degree = 2)
  )
  expect_true(all(basis_int(b, matrix(b@lower, 1L)) == 0))

  x <- cbind(runif(6, -1, 2), runif(6, 3, 5))
  want <- kron_rows(list(
    basis_int(bspline_basis(-1, 2, dimension = 4), x[, 1]),
    basis_int(bspline_basis(3, 5, dimension = 3, degree = 2), x[, 2])
  ))
  expect_equal(basis_int(b, x), want, ignore_attr = TRUE)
})


test_that("the Gram matrix is the Kronecker product, and matches a box rule", {
  b <- tb2()
  m1 <- bspline_basis(dimension = 4)
  m2 <- fourier_basis(dimension = 3)

  for (o in list(c(0, 0), c(1, 0), c(0, 1), c(2, 1))) {
    expect_equal(basis_gram(b, order = o),
      kronecker(basis_gram(m1, o[1]), basis_gram(m2, o[2])),
      ignore_attr = TRUE, label = paste(o, collapse = ",")
    )
  }
  # and against a product quadrature over the box, which knows nothing about
  # separability
  expect_equal(basis_gram(b), basis7:::numerical_gram(b, 0L, panels = 400L),
    tolerance = 1e-8
  )
})


test_that("a product of products is flattened, and a product of one is itself", {
  b1 <- bspline_basis(dimension = 4)
  b2 <- fourier_basis(dimension = 3)
  b3 <- poly_basis(dimension = 2)

  nested <- tensor_basis(tensor_basis(b1, b2), b3)
  flat <- tensor_basis(b1, b2, b3)
  expect_length(nested@marginals, 3L)
  expect_identical(basis_nvar(nested), 3L)
  expect_identical(basis_colnames(nested), basis_colnames(flat))

  # a list is accepted in place of the arguments
  expect_identical(basis_colnames(tensor_basis(list(b1, b2, b3))),
    basis_colnames(flat)
  )
  # and a single basis comes back unchanged rather than wrapped
  expect_true(S7::S7_inherits(tensor_basis(b1), BsplineBasis))
})


test_that("malformed products and points are refused", {
  b1 <- bspline_basis(dimension = 4)
  expect_error(tensor_basis(), "at least one")
  expect_error(tensor_basis(b1, 1), "must be a basis")
  expect_error(tensor_basis(tb2(), b1) |> tensor_basis(tb2()), NA)

  b <- tb2()
  expect_error(basis_eval(b, cbind(0.5, 0.5, 0.5)), "2 columns")
  expect_error(basis_eval(b, cbind(1.5, 0.5)), "variable 1")
  expect_error(basis_eval(b, cbind(0.5, -1)), "variable 2")
  # a marginal must take one variable
  expect_error(tensor_basis(b, b1), NA) # a product is flattened, so this is fine
})


test_that("a product of partitions of unity is one", {
  b <- tensor_basis(
    bspline_basis(dimension = 5), bspline_basis(dimension = 4, degree = 2)
  )
  set.seed(44)
  x <- cbind(runif(30), runif(30))
  expect_equal(rowSums(basis_eval(b, x)), rep(1, 30), tolerance = 1e-12)
  expect_true(basis7:::basis_partitions_unity(b))
  # and a product with a family that does not claim it does not claim it
  expect_false(basis7:::basis_partitions_unity(tb2()))
})


test_that("a product is exact exactly when its marginals are", {
  expect_false(any(basis_is_numerical(tb2())))

  Bare <- S7::new_class("BareMarginal", parent = basis)
  S7::method(basis_eval, Bare) <- function(basis, x, ...) {
    out <- cbind(1, x)
    colnames(out) <- basis_colnames(basis)
    out
  }
  bare <- Bare(
    basis_name = "bare", dimension = 2L, lower = 0, upper = 1,
    basis_params = list()
  )
  mixed <- tensor_basis(bspline_basis(dimension = 3, degree = 2), bare)
  expect_true(all(basis_is_numerical(mixed)))
})


test_that("check_basis passes on products of every combination", {
  for (b in list(
    tb2(),
    tensor_basis(bspline_basis(dimension = 5), bspline_basis(dimension = 4)),
    tensor_basis(poly_basis(dimension = 4), fourier_basis(dimension = 3)),
    tensor_basis(
      bspline_basis(-2, 3, dimension = 4),
      poly_basis(1, 4, dimension = 3)
    ),
    tensor_basis(
      bspline_basis(dimension = 3, degree = 2),
      fourier_basis(dimension = 3), poly_basis(dimension = 3)
    )
  )) {
    r <- check_basis(b, n = 21L, verbose = FALSE)
    expect_true(all(r, na.rm = TRUE), label = b@basis_name)
    expect_false(isFALSE(r[["shape"]]), label = b@basis_name)
  }
})


test_that("check_basis accepts a point count given as a double", {
  # dim() is always integer, so a count given as 21 rather than 21L would
  # otherwise fail a check about the basis for a reason about the argument.
  expect_true(check_basis(tb2(), n = 21, verbose = FALSE)[["shape"]])
  expect_true(check_basis(bspline_basis(dimension = 5), n = 25,
    verbose = FALSE
  )[["shape"]])
})


test_that("plot refuses a product, and says why", {
  expect_error(plot(tb2()), "one variable")
})
