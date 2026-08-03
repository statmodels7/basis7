# The Legendre family.

test_that("the recurrence reproduces the closed-form polynomials", {
  p <- poly_basis(lower = 0, upper = 1, dimension = 6)
  x <- seq(0, 1, length.out = 21)
  t <- 2 * x - 1
  expect_equal(basis_eval(p, x),
    cbind(
      1, t, (3 * t^2 - 1) / 2, (5 * t^3 - 3 * t) / 2,
      (35 * t^4 - 30 * t^2 + 3) / 8,
      (63 * t^5 - 70 * t^3 + 15 * t) / 8
    ),
    ignore_attr = TRUE
  )
  expect_identical(basis_colnames(p), paste0("P", 0:5))
})


test_that("the Gram matrix is the diagonal the orthogonality gives", {
  for (lim in list(c(0, 1), c(-2, 3))) {
    p <- poly_basis(lim[1], lim[2], dimension = 6)
    g <- basis_gram(p)
    expect_equal(g[upper.tri(g)], rep(0, sum(upper.tri(g))), ignore_attr = TRUE)
    expect_equal(diag(g), (lim[2] - lim[1]) / (2 * (0:5) + 1),
      ignore_attr = TRUE
    )
    # and against an independent quadrature
    expect_equal(g, basis7:::numerical_gram(p, 0L, panels = 200L, nodes = 10L),
      tolerance = 1e-8
    )
  }
})


test_that("derivatives of every order are exact", {
  p <- poly_basis(dimension = 6)
  x <- seq(0.05, 0.95, length.out = 13)
  for (d in 1:5) {
    ana <- basis_deriv(p, x, order = d)
    ref <- basis7:::numerical_deriv_matrix(
      function(z) basis_deriv(p, z, order = d - 1L), x, 1L, 0, 1
    )
    s <- pmax(abs(ana), abs(ref))
    big <- s > 1e-8 * max(s)
    expect_lt(max(abs(ana - ref)[big] / s[big]), 1e-6)
  }
  # above the highest degree everything is zero
  expect_true(all(basis_deriv(p, x, order = 6L) == 0))
  expect_true(all(basis_gram(p, order = 6L) == 0))
})


test_that("the integral is exact and anchored", {
  p <- poly_basis(lower = -1, upper = 4, dimension = 5)
  expect_true(all(basis_int(p, -1) == 0))
  q <- 2.4
  ref <- vapply(seq_len(5), function(j) {
    stats::integrate(function(t) basis_eval(p, t)[, j], -1, q,
      rel.tol = 1e-11
    )$value
  }, numeric(1))
  expect_equal(drop(basis_int(p, q)), ref, tolerance = 1e-9, ignore_attr = TRUE)
})


test_that("the basis spans the polynomials it claims to", {
  # The point of choosing Legendre over raw powers is conditioning, not span.
  p <- poly_basis(dimension = 5)
  x <- seq(0, 1, length.out = 60)
  for (deg in 0:4) {
    fit <- stats::lm.fit(basis_eval(p, x), x^deg)$fitted.values
    expect_equal(fit, x^deg, tolerance = 1e-10, ignore_attr = TRUE)
  }
  # and the conditioning is the reason: raw powers give a Hilbert matrix
  raw <- outer(x, 0:9, "^")
  leg <- basis_eval(poly_basis(dimension = 10), x)
  expect_lt(kappa(crossprod(leg)), kappa(crossprod(raw)) / 1e4)
})


test_that("check_basis passes on the family", {
  r <- check_basis(poly_basis(dimension = 6), verbose = FALSE)
  expect_true(all(r, na.rm = TRUE))
  expect_true(is.na(r[["partition"]]))
})
