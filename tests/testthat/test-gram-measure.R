# basis_gram against measures other than Lebesgue.

test_that("a uniform weight reproduces the Lebesgue matrix", {
  for (b in list(
    bspline_basis(dimension = 6), fourier_basis(dimension = 5),
    poly_basis(dimension = 5)
  )) {
    expect_equal(
      basis_gram(b, weight = function(t) rep(1, length(t))),
      basis_gram(b),
      tolerance = 1e-9, label = b@basis_name
    )
  }
})


test_that("a non-uniform weight matches an independent quadrature", {
  b <- bspline_basis(dimension = 5)
  w <- function(t) 2 * t
  got <- basis_gram(b, weight = w)
  want <- outer(seq_len(5), seq_len(5), Vectorize(function(i, j) {
    stats::integrate(
      function(t) basis_eval(b, t)[, i] * basis_eval(b, t)[, j] * w(t),
      0, 1,
      rel.tol = 1e-11
    )$value
  }))
  expect_equal(got, want, tolerance = 1e-9, ignore_attr = TRUE)
})


test_that("the empirical measure is B'B / n", {
  set.seed(11)
  b <- bspline_basis(dimension = 6)
  x <- runif(500)
  expect_equal(basis_gram(b, at = x), crossprod(basis_eval(b, x)) / 500,
    ignore_attr = TRUE
  )
  # it is a probability measure, so a large uniform sample approaches the
  # Lebesgue matrix divided by the width of the interval
  big <- runif(2e5)
  expect_equal(basis_gram(b, at = big), basis_gram(b), tolerance = 0.02)
})


test_that("the empirical measure works at any derivative order", {
  set.seed(12)
  b <- fourier_basis(dimension = 5)
  x <- runif(300)
  for (d in 0:2) {
    expect_equal(basis_gram(b, order = d, at = x),
      crossprod(basis_deriv(b, x, order = d)) / 300,
      ignore_attr = TRUE, label = paste("order", d)
    )
  }
})


test_that("missing points are dropped from the empirical measure", {
  b <- bspline_basis(dimension = 5)
  x <- c(0.1, 0.4, NA, 0.9)
  expect_equal(basis_gram(b, at = x), basis_gram(b, at = c(0.1, 0.4, 0.9)))
  expect_error(basis_gram(b, at = c(NA_real_, NA_real_)), "no usable points")
})


test_that("the two measures cannot be asked for at once, and are validated", {
  b <- bspline_basis(dimension = 5)
  expect_error(
    basis_gram(b, at = c(0.2, 0.5), weight = function(t) rep(1, length(t))),
    "at most one"
  )
  expect_error(basis_gram(b, weight = 3), "must be a function")
  expect_error(basis_gram(b, weight = function(t) -1), "non-negative")
  expect_error(basis_gram(b, weight = function(t) 1), "one .* value per point")
  expect_error(basis_gram(b, at = 2), "outside the basis interval")
})


test_that("a transformed basis carries the measures through", {
  set.seed(13)
  b <- bspline_basis(dimension = 6)
  o <- orthonorm_basis(b)
  x <- runif(400)
  expect_equal(basis_gram(o, at = x), crossprod(basis_eval(o, x)) / 400,
    ignore_attr = TRUE
  )
})
