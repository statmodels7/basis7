# The B-spline family: what it promises beyond the common contract.

test_that("dimension, degree and knot count are consistent", {
  for (degree in 0:4) {
    for (k in (degree + 1L):(degree + 6L)) {
      b <- bspline_basis(dimension = k, degree = degree)
      expect_identical(b@dimension, as.integer(k))
      expect_length(b@basis_params$knots, k - degree - 1L)
      expect_identical(ncol(basis_eval(b, c(0, 0.5, 1))), as.integer(k))
    }
  }
})


test_that("a dimension too small for the degree is refused, at the boundary", {
  # degree + 1 is the smallest basis that exists: the polynomials of that
  # degree, with no interior knot.
  expect_silent(bspline_basis(dimension = 4, degree = 3))
  expect_error(bspline_basis(dimension = 3, degree = 3), "too small")
  expect_error(bspline_basis(dimension = 1, degree = 1), "too small")
  expect_error(bspline_basis(degree = -1), "non-negative integer")
  expect_error(bspline_basis(degree = 2.5), "non-negative integer")
})


test_that("the basis is a partition of unity everywhere", {
  for (degree in 0:3) {
    b <- bspline_basis(dimension = degree + 5L, degree = degree)
    x <- seq(0, 1, length.out = 101)
    expect_equal(rowSums(basis_eval(b, x)), rep(1, 101), tolerance = 1e-12)
  }
})


test_that("derivatives above the degree are zero, and it is the zero matrix", {
  b <- bspline_basis(dimension = 7, degree = 3)
  x <- seq(0, 1, length.out = 11)
  expect_true(all(basis_deriv(b, x, order = 4L) == 0))
  expect_true(all(basis_deriv(b, x, order = 9L) == 0))
  expect_identical(dim(basis_deriv(b, x, order = 4L)), c(11L, 7L))
  # and the Gram matrix of a zero derivative is the zero matrix
  expect_true(all(basis_gram(b, order = 4L) == 0))
})


test_that("the exact Gram matrix agrees with an independent quadrature", {
  # The reference uses a different rule from the one the method uses, so the
  # comparison is not the same arithmetic twice.
  for (degree in 1:3) {
    b <- bspline_basis(dimension = degree + 5L, degree = degree)
    for (d in 0:degree) {
      exact <- basis_gram(b, order = d)
      ref <- basis7:::numerical_gram(b, d, panels = 600L, nodes = 9L)
      expect_equal(exact, ref,
        tolerance = 1e-6,
        label = sprintf("degree %d order %d", degree, d)
      )
    }
  }
})


test_that("the Gram matrix is the quadratic form of the integrated square", {
  # The point of the matrix: beta' G_d beta is the integral of the squared
  # d-th derivative of the function the coefficients describe.
  set.seed(1)
  b <- bspline_basis(dimension = 8, degree = 3)
  beta <- rnorm(8)
  for (d in 0:2) {
    f <- function(t) (basis_deriv(b, t, order = d) %*% beta)^2
    ref <- stats::integrate(function(t) apply(cbind(t), 1, f),
      lower = 0, upper = 1, rel.tol = 1e-10
    )$value
    expect_equal(drop(t(beta) %*% basis_gram(b, order = d) %*% beta), ref,
      tolerance = 1e-6, label = paste("order", d)
    )
  }
})


test_that("the integral is the antiderivative anchored at the lower endpoint", {
  b <- bspline_basis(lower = -1, upper = 2, dimension = 7)
  q <- 0.83
  ref <- vapply(seq_len(b@dimension), function(j) {
    stats::integrate(function(t) basis_eval(b, t)[, j], -1, q,
      rel.tol = 1e-11
    )$value
  }, numeric(1))
  expect_equal(drop(basis_int(b, q)), ref,
    tolerance = 1e-8, ignore_attr = TRUE
  )
})


test_that("a shifted interval behaves like the unit one, rescaled", {
  b0 <- bspline_basis(0, 1, dimension = 6)
  b1 <- bspline_basis(2, 6, dimension = 6)
  u <- seq(0, 1, length.out = 9)
  expect_equal(basis_eval(b0, u), basis_eval(b1, 2 + 4 * u))
  # the first derivative scales by the inverse width
  expect_equal(basis_deriv(b0, u, 1L), 4 * basis_deriv(b1, 2 + 4 * u, 1L))
})
