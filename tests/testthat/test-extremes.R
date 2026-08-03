# The interior is where the tests are comfortable and the defects are not.
# Every check here sits at an endpoint, at a knot, or at a degenerate size.

test_that("evaluation at the exact endpoints works for every family", {
  for (b in list(
    bspline_basis(dimension = 6, degree = 3),
    bspline_basis(dimension = 4, degree = 0),
    bspline_basis(lower = -3, upper = -1, dimension = 5),
    fourier_basis(dimension = 5),
    fourier_basis(lower = -2, upper = 7, dimension = 3)
  )) {
    for (x in c(b@lower, b@upper)) {
      m <- basis_eval(b, x)
      expect_true(all(is.finite(m)), label = b@basis_name)
      expect_identical(dim(m), c(1L, b@dimension), label = b@basis_name)
    }
  }
})


test_that("the numerical derivative works AT the endpoints", {
  # A symmetric stencil centred on an endpoint asks for points outside the
  # interval, where the basis is not defined. Those points get a one-sided
  # stencil instead; without it the step would be clamped to zero and the
  # result would be a division by zero.
  Quad <- S7::new_class("QuadEdge", parent = basis)
  S7::method(basis_eval, Quad) <- function(basis, x, ...) {
    out <- cbind(1, x, x^2)
    colnames(out) <- basis_colnames(basis)
    out
  }
  q <- Quad(
    basis_name = "quad", dimension = 3L, lower = 0, upper = 1,
    basis_params = list()
  )

  for (x in c(0, 1)) {
    d1 <- basis_deriv(q, x, order = 1L)
    expect_true(all(is.finite(d1)), label = paste("x =", x))
    expect_equal(drop(d1), c(0, 1, 2 * x), tolerance = 1e-6, ignore_attr = TRUE)
    d2 <- basis_deriv(q, x, order = 2L)
    expect_true(all(is.finite(d2)), label = paste("x =", x))
    expect_equal(drop(d2), c(0, 0, 2), tolerance = 1e-5, ignore_attr = TRUE)
  }

  # and a whole vector of endpoints and interior points at once
  x <- c(0, 0.5, 1)
  expect_true(all(is.finite(basis_deriv(q, x, order = 1L))))
})


test_that("the integral over the whole interval is finite at the upper end", {
  for (b in list(
    bspline_basis(dimension = 7),
    fourier_basis(lower = -1, upper = 4, dimension = 5)
  )) {
    total <- basis_int(b, b@upper)
    expect_true(all(is.finite(total)), label = b@basis_name)
  }
  # for a B-spline partition of unity the total is the interval width
  b <- bspline_basis(lower = 2, upper = 6, dimension = 7)
  expect_equal(sum(basis_int(b, b@upper)), 4, tolerance = 1e-10)
})


test_that("a degree-zero basis, whose functions are discontinuous, still works", {
  b <- bspline_basis(dimension = 4, degree = 0)
  x <- seq(0, 1, length.out = 21)
  expect_equal(rowSums(basis_eval(b, x)), rep(1, 21), tolerance = 1e-12)
  expect_true(all(basis_deriv(b, x, order = 1L) == 0))
  # the Gram matrix of indicator functions is diagonal, each entry a cell width
  g <- basis_gram(b)
  expect_equal(g[upper.tri(g)], rep(0, 6), ignore_attr = TRUE)
  expect_equal(diag(g), rep(0.25, 4), ignore_attr = TRUE)
})


test_that("the smallest basis of each family works", {
  b <- bspline_basis(dimension = 1, degree = 0)
  expect_equal(drop(basis_eval(b, c(0, 0.5, 1))), c(1, 1, 1),
    ignore_attr = TRUE
  )
  expect_equal(drop(basis_int(b, 1)), 1, ignore_attr = TRUE)

  f <- fourier_basis(dimension = 1)
  expect_identical(basis_colnames(f), "const")
  expect_equal(drop(basis_eval(f, c(0, 1))), c(1, 1), ignore_attr = TRUE)
  expect_equal(drop(basis_gram(f)), 1, ignore_attr = TRUE)
  expect_equal(drop(basis_gram(f, order = 1L)), 0, ignore_attr = TRUE)
  expect_equal(drop(basis_int(f, 1)), 1, ignore_attr = TRUE)
})


test_that("evaluation exactly at an interior knot is well defined", {
  b <- bspline_basis(dimension = 6, degree = 3)
  k <- b@basis_params$knots
  m <- basis_eval(b, k)
  expect_true(all(is.finite(m)))
  expect_equal(rowSums(m), rep(1, length(k)), tolerance = 1e-12)
})


test_that("a very wide and a very narrow interval both behave", {
  wide <- bspline_basis(lower = -1e6, upper = 1e6, dimension = 6)
  x <- seq(-1e6, 1e6, length.out = 9)
  expect_equal(rowSums(basis_eval(wide, x)), rep(1, 9), tolerance = 1e-10)
  expect_true(all(is.finite(basis_deriv(wide, x, order = 1L))))
  expect_true(all(is.finite(basis_gram(wide, order = 2L))))

  narrow <- bspline_basis(lower = 0, upper = 1e-6, dimension = 6)
  y <- seq(0, 1e-6, length.out = 9)
  expect_equal(rowSums(basis_eval(narrow, y)), rep(1, 9), tolerance = 1e-10)
  expect_true(all(is.finite(basis_deriv(narrow, y, order = 1L))))
})


test_that("a high-frequency Fourier basis stays exact", {
  b <- fourier_basis(dimension = 41)
  x <- seq(0, 1, length.out = 17)
  m <- basis_eval(b, x)
  expect_true(all(is.finite(m)))
  expect_equal(m[, "sin20"], sin(40 * pi * x))
  # the Gram matrix stays diagonal at high order despite the large scaling
  g <- basis_gram(b, order = 3L)
  expect_equal(g[upper.tri(g)], rep(0, sum(upper.tri(g))), ignore_attr = TRUE)
  expect_true(all(is.finite(g)))
})
