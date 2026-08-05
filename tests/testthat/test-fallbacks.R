# The bargain: a basis that implements only its evaluation is complete.

# Gaussian bumps, with everything but basis_eval left to the base class.
Bumps <- S7::new_class("Bumps", parent = basis)
S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
  ctr <- basis@basis_params$centres
  s <- basis@basis_params$width
  out <- exp(-0.5 * outer(x, ctr, "-")^2 / s^2)
  colnames(out) <- basis_colnames(basis)
  out
}
bumps <- function(k = 5, lower = 0, upper = 1) {
  Bumps(
    basis_name = "bumps", dimension = as.integer(k),
    lower = lower, upper = upper,
    basis_params = list(
      centres = seq(lower + 0.1, upper - 0.1, length.out = k),
      width = 0.15
    )
  )
}

# The closed forms, for comparison: d/dx exp(-(x-c)^2/2s^2) = -(x-c)/s^2 * f.
bump_d1 <- function(b, x) {
  ctr <- b@basis_params$centres
  s <- b@basis_params$width
  -outer(x, ctr, "-") / s^2 * basis_eval(b, x)
}
bump_d2 <- function(b, x) {
  ctr <- b@basis_params$centres
  s <- b@basis_params$width
  z <- outer(x, ctr, "-")
  (z^2 / s^4 - 1 / s^2) * basis_eval(b, x)
}


test_that("all four generics work from the evaluation alone", {
  b <- bumps()
  x <- seq(0, 1, length.out = 9)
  expect_identical(dim(basis_eval(b, x)), c(9L, 5L))
  expect_identical(dim(basis_deriv(b, x, order = 1L)), c(9L, 5L))
  expect_identical(dim(basis_int(b, x)), c(9L, 5L))
  expect_identical(dim(basis_gram(b)), c(5L, 5L))
})


test_that("the object reports which of its methods are numerical", {
  n <- basis_is_numerical(bumps())
  expect_true(all(n))
  expect_named(n, c("basis_deriv", "basis_int", "basis_gram"))

  # a family that supplies all three reports none
  expect_false(any(basis_is_numerical(bspline_basis(dimension = 5))))
  # and one that supplies some reports the rest
  expect_false(basis_is_numerical(fourier_basis(dimension = 5))[["basis_gram"]])
})


test_that("the numerical derivatives match the closed forms", {
  b <- bumps()
  x <- seq(0.05, 0.95, length.out = 15)
  expect_equal(basis_deriv(b, x, order = 1L), bump_d1(b, x),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  expect_equal(basis_deriv(b, x, order = 2L), bump_d2(b, x),
    tolerance = 1e-4, ignore_attr = TRUE
  )
})


test_that("the numerical integral matches integrate() and starts at zero", {
  b <- bumps()
  expect_true(all(basis_int(b, b@lower) == 0))
  q <- 0.73
  ref <- vapply(seq_len(b@dimension), function(j) {
    stats::integrate(function(t) basis_eval(b, t)[, j], 0, q)$value
  }, numeric(1))
  expect_equal(drop(basis_int(b, q)), ref,
    tolerance = 1e-7, ignore_attr = TRUE
  )
})


test_that("the numerical integral is right for an unsorted, repeated x", {
  # The method integrates over the sorted unique points and maps back, so the
  # order it is given must not matter and a repeat must not shift anything.
  b <- bumps()
  x <- c(0.6, 0.2, 0.6, 0, 0.9)
  got <- basis_int(b, x)
  want <- basis_int(b, sort(unique(x)))[c(3L, 2L, 3L, 1L, 4L), , drop = FALSE]
  expect_equal(got, want)
})


test_that("one stencil is used, never a chain of first differences", {
  # A basis whose third derivative is exactly zero: nested differencing returns
  # noise of order one, a single stencil returns zero. This is the property the
  # whole fallback design turns on.
  Quad <- S7::new_class("Quad", parent = basis)
  S7::method(basis_eval, Quad) <- function(basis, x, ...) {
    out <- cbind(1, x, x^2)
    colnames(out) <- basis_colnames(basis)
    out
  }
  q <- Quad(
    basis_name = "quad", dimension = 3L, lower = -1, upper = 1,
    basis_params = list()
  )
  x <- seq(-0.9, 0.9, length.out = 11)

  expect_equal(basis_deriv(q, x, order = 1L), cbind(0, 1, 2 * x),
    tolerance = 1e-7, ignore_attr = TRUE
  )
  expect_equal(basis_deriv(q, x, order = 2L),
    cbind(rep(0, length(x)), 0, 2),
    tolerance = 1e-6, ignore_attr = TRUE
  )
  expect_lt(max(abs(basis_deriv(q, x, order = 3L))), 1e-3)
})


test_that("the finite-difference weights are exact on polynomials", {
  # The weights come from a Vandermonde solve, so the classic stencils must
  # fall out of it rather than being coincidences.
  expect_equal(numericals7::fd_weights(c(-1, 0, 1), 1L), c(-0.5, 0, 0.5))
  expect_equal(numericals7::fd_weights(c(-1, 0, 1), 2L), c(1, -2, 1))
  expect_equal(numericals7::fd_weights(c(-2, -1, 0, 1, 2), 3L),
    c(-0.5, 1, 0, -1, 0.5),
    tolerance = 1e-10
  )
  expect_equal(numericals7::fd_weights(c(-2, -1, 0, 1, 2), 4L),
    c(1, -4, 6, -4, 1),
    tolerance = 1e-10
  )
  # a one-sided stencil is exact too
  w <- numericals7::fd_weights(0:2, 1L)
  expect_equal(w, c(-1.5, 2, -0.5))
})


test_that("Gauss-Legendre integrates polynomials of the right degree exactly", {
  for (n in 1:8) {
    gl <- basis7:::gauss_legendre(n)
    expect_length(gl$nodes, n)
    expect_equal(sum(gl$weights), 2)
    # exact up to degree 2n - 1
    p <- 2 * n - 1
    got <- sum(gl$weights * gl$nodes^p)
    want <- if (p %% 2 == 0) 2 / (p + 1) else 0
    expect_equal(got, want, tolerance = 1e-10, label = paste("n =", n))
  }
})


test_that("the numerical derivative does not label its rows after the stencil", {
  # The stencil each point uses is chosen per point and held in a character
  # vector. Naming the offsets after it leaks those names through the evaluation
  # points and out as the row names of the result, for any basis whose method
  # propagates the names of x -- outer() does, and so does the one in the
  # vignette.
  Named <- S7::new_class("Named", parent = basis, package = NULL)
  S7::method(basis_eval, Named) <- function(basis, x, ...) {
    out <- exp(-0.5 * outer(x, c(0.25, 0.5, 0.75), "-")^2 / 0.2^2)
    colnames(out) <- basis_colnames(basis)
    out
  }
  nb <- Named(
    basis_name = "named", dimension = 3L, lower = 0, upper = 1,
    basis_params = list()
  )

  # Endpoints included, so the one-sided stencils are exercised too.
  x <- c(0, 0.3, 0.5, 1)
  for (d in 1:3) {
    expect_null(rownames(basis_deriv(nb, x, order = d)),
      label = paste("order", d)
    )
  }
  expect_null(rownames(basis_int(nb, x)))
})
