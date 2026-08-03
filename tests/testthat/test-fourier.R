# The Fourier family: the shift identity, and the integral convention that the
# identity alone does not deliver.

test_that("an even dimension is refused rather than adjusted", {
  # Growing the dimension silently would hand back a basis of a size the caller
  # did not ask for, and the constructor is the only place the inconsistency
  # can be caught.
  expect_error(fourier_basis(dimension = 4), "must be odd")
  expect_error(fourier_basis(dimension = 2), "must be odd")
  expect_silent(fourier_basis(dimension = 1))
  expect_silent(fourier_basis(dimension = 7))
  expect_error(fourier_basis(omega = 0), "positive")
  expect_error(fourier_basis(omega = -1), "positive")
})


test_that("the columns are the constant and the sine-cosine pairs", {
  b <- fourier_basis(dimension = 5)
  expect_identical(basis_colnames(b), c("const", "sin1", "cos1", "sin2", "cos2"))

  x <- seq(0, 1, length.out = 21)
  m <- basis_eval(b, x)
  expect_equal(m[, "const"], rep(1, 21))
  expect_equal(m[, "sin1"], sin(2 * pi * x))
  expect_equal(m[, "cos2"], cos(4 * pi * x))

  # a dimension of one is the constant alone
  expect_identical(basis_colnames(fourier_basis(dimension = 1)), "const")
  expect_equal(drop(basis_eval(fourier_basis(dimension = 1), c(0, 1))), c(1, 1))
})


test_that("derivatives follow the shift identity at every order", {
  b <- fourier_basis(lower = -1, upper = 3, dimension = 7)
  omega <- b@basis_params$omega
  x <- seq(-1, 3, length.out = 13)
  z <- 2 * pi * (x - b@lower) / omega

  for (d in 1:6) {
    m <- basis_deriv(b, x, order = d)
    expect_equal(m[, "const"], rep(0, 13), label = paste("order", d))
    for (j in 1:3) {
      s <- (2 * pi * j / omega)^d
      expect_equal(m[, paste0("sin", j)], sin(j * z + d * pi / 2) * s,
        label = sprintf("order %d sin%d", d, j)
      )
      expect_equal(m[, paste0("cos", j)], cos(j * z + d * pi / 2) * s,
        label = sprintf("order %d cos%d", d, j)
      )
    }
  }
})


test_that("the integral is anchored at the lower endpoint, on every column", {
  # The shift identity at order -1 gives AN antiderivative, and it is not this
  # one: at the lower endpoint its sine columns are -omega/(2 pi j) while its
  # cosine columns are already zero. Half the columns following one convention
  # and half another is invisible when lower = 0, which is why the check runs
  # on a shifted interval.
  b <- fourier_basis(lower = -1, upper = 3, dimension = 7)
  expect_true(all(basis_int(b, b@lower) == 0))

  q <- 1.37
  ref <- vapply(seq_len(b@dimension), function(j) {
    stats::integrate(function(t) basis_eval(b, t)[, j], b@lower, q,
      rel.tol = 1e-11
    )$value
  }, numeric(1))
  expect_equal(drop(basis_int(b, q)), ref,
    tolerance = 1e-9, ignore_attr = TRUE
  )

  # and the constant column integrates to the elapsed length, not to x
  expect_equal(unname(basis_int(b, q)[, "const"]), q - b@lower)
})


test_that("the Gram matrix is diagonal over a whole period", {
  b <- fourier_basis(lower = 2, upper = 9, dimension = 7)
  omega <- b@basis_params$omega
  for (d in 0:3) {
    g <- basis_gram(b, order = d)
    expect_equal(g[upper.tri(g)], rep(0, sum(upper.tri(g))),
      label = paste("order", d)
    )
    j <- 1:3
    expected <- c(
      if (d == 0L) omega else 0,
      rep(omega / 2 * (2 * pi * j / omega)^(2 * d), each = 2L)
    )
    expect_equal(diag(g), expected, ignore_attr = TRUE, label = paste("order", d))
  }
})


test_that("a period other than the interval width falls back to quadrature", {
  # Orthogonality is a property of a whole period; on a fraction of one the
  # closed form would be wrong, so the basis says so rather than using it.
  b <- fourier_basis(lower = 0, upper = 1, dimension = 5, omega = 0.7)
  expect_false(b@basis_params$full_period)
  g <- basis_gram(b)
  ref <- basis7:::numerical_gram(b, 0L, panels = 400L, nodes = 10L)
  expect_equal(g, ref, tolerance = 1e-8)
  # and it is genuinely not diagonal
  expect_gt(max(abs(g[upper.tri(g)])), 1e-3)
})


test_that("check_basis passes on the family and reports partition as not claimed", {
  r <- check_basis(fourier_basis(dimension = 7), verbose = FALSE)
  expect_true(r[["shape"]])
  expect_true(r[["deriv"]])
  expect_true(r[["integral"]])
  expect_true(r[["gram"]])
  expect_true(is.na(r[["partition"]]))
})
