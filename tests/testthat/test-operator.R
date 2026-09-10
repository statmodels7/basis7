# Linear differential operators: the object, its null space, and the
# roughness matrix it induces on a basis.
#
# The null space is asserted against the DEFINITION -- L applied to it is
# zero -- rather than against the table the package builds, so a mistake in
# the root clustering cannot confirm itself.

set.seed(4)

test_that("the three constructors carry the order and the weights", {
  expect_identical(operator_order(deriv_operator(2)), 2L)
  expect_identical(operator_weights(deriv_operator(3)), c(0, 0, 0))

  expect_identical(operator_order(harmonic_operator(365)), 3L)
  expect_identical(operator_order(harmonic_operator(365, harmonics = 2)), 5L)
  expect_identical(operator_order(oscillator_operator(365)), 2L)
  expect_identical(operator_order(linear_operator(c(1, 2))), 2L)

  # the order is always the dimension of the null space
  for (op in list(deriv_operator(4), harmonic_operator(12),
                  harmonic_operator(12, harmonics = 3),
                  oscillator_operator(7), linear_operator(c(0.5, -1)))) {
    expect_identical(nrow(operator_null(op)), operator_order(op))
  }
})

test_that("the harmonic weights are the ones fda's vec2Lfd is given", {
  # fda writes c(0, (2*pi/T)^2, 0) for the harmonic acceleration operator
  nu2 <- (2 * pi / 365)^2
  expect_equal(operator_weights(harmonic_operator(365)), c(0, nu2, 0))
  expect_equal(operator_weights(linear_operator(c(0, nu2, 0))),
               operator_weights(harmonic_operator(365)))
})

test_that("the harmonic operator is the oscillator with a leading derivative", {
  for (h in 1:3) {
    expect_equal(operator_weights(deriv_operator(1) * oscillator_operator(365, h)),
                 operator_weights(harmonic_operator(365, harmonics = h)))
  }
})

test_that("L applied to its own null space is zero", {
  # THE DEFINITION, checked through a basis that represents the null space
  # exactly: the first three columns of a Fourier basis ARE 1, sin and cos of
  # the fundamental, so L of them must vanish pointwise.
  b <- fourier_basis(lower = 0, upper = 365, dimension = 9)
  xs <- seq(0, 365, length.out = 41)
  lb <- basis7:::operator_eval(b, xs, harmonic_operator(365))
  # RELATIVE, because |L b| scales as the cube of the frequency and an
  # absolute bound would be a statement about the period
  expect_lt(max(abs(lb[, 1:3])) / max(abs(lb[, 4:9])), 1e-12)
  # and the columns above the null space are not themselves zero, or
  # the check would be vacuous
  expect_gt(max(abs(lb[, 4:9])), 0)

  # two harmonics kill five columns
  lb2 <- basis7:::operator_eval(b, xs, harmonic_operator(365, harmonics = 2))
  expect_lt(max(abs(lb2[, 1:5])) / max(abs(lb2[, 6:9])), 1e-12)
  expect_gt(max(abs(lb2[, 6:9])), 0)
})

test_that("composition adds the orders and unites the null spaces", {
  a <- deriv_operator(2)
  o <- oscillator_operator(365)
  p <- a * o
  expect_identical(operator_order(p), 4L)
  expect_identical(sort(operator_null(p)$label),
                   sort(c(operator_null(a)$label, operator_null(o)$label)))
  # and it commutes
  expect_equal(operator_weights(o * a), operator_weights(p))
})

test_that("composing with harmonic_operator raises the order by one", {
  # The documented trap: harmonic_operator() already carries a leading D, so
  # this is D^3(D^2 + nu^2) and its null space holds t^2.
  p <- deriv_operator(2) * harmonic_operator(365)
  expect_identical(operator_order(p), 5L)
  expect_true("t^2" %in% operator_null(p)$label)
  expect_false("t^2" %in% operator_null(deriv_operator(2) *
                                          oscillator_operator(365))$label)
})

test_that("a null-space column is the function its label names", {
  x <- seq(0, 365, length.out = 30)
  d <- operator_null_design(harmonic_operator(365), x)
  expect_equal(d[, 1], rep(1, length(x)))
  expect_equal(d[, 2], sin(2 * pi * x / 365))
  expect_equal(d[, 3], cos(2 * pi * x / 365))

  # a pair must not come out as the same function twice, which is what a
  # sine/cosine decided by the row's parity would give when a real root
  # comes first
  p <- deriv_operator(1) * oscillator_operator(365)
  dp <- operator_null_design(p, x)
  expect_gt(max(abs(dp[, 2] - dp[, 3])), 0.5)
})

test_that("the period may be left to the interval and resolved later", {
  op <- harmonic_operator()
  expect_true(all(is.na(operator_weights(op))))
  expect_error(operator_null(op), "period is NULL")
  expect_error(op * deriv_operator(1), "cannot be composed")

  r <- operator_resolve(op, 0, 365)
  expect_equal(operator_weights(r), operator_weights(harmonic_operator(365)))
  # an operator given a period keeps it whatever the interval is
  expect_equal(operator_weights(operator_resolve(harmonic_operator(12), 0, 365)),
               operator_weights(harmonic_operator(12)))
  # and one already resolved is returned untouched
  expect_identical(operator_resolve(deriv_operator(2), 0, 1), deriv_operator(2))
})

test_that("the constructors refuse what they cannot build", {
  expect_error(harmonic_operator(-1), "positive number")
  expect_error(harmonic_operator(c(1, 2)), "positive number")
  expect_error(harmonic_operator(365, harmonics = 0), "at least 1")
  expect_error(linear_operator(numeric(0)), "at least one finite entry")
  expect_error(linear_operator(c(1, NA)), "at least one finite entry")
  expect_error(deriv_operator(0), "at least 1")
  expect_error(operator_order("D2"), "linear differential operator")
})

test_that("an operator prints its formula and its null space", {
  out <- capture.output(print(harmonic_operator(365)))
  expect_match(out[[1]], "order 3")
  expect_match(out[[2]], "D\\^3 x")
  expect_match(out[[3]], "null space: 1, sin")
  # an unresolved one says so rather than printing weights it does not have
  expect_match(capture.output(print(harmonic_operator()))[[3]], "the smoother's interval")
})


# ---------------------------------------------------------------------------
# The roughness matrix
# ---------------------------------------------------------------------------

test_that("basis_gram routes an operator and reproduces the derivative case", {
  b <- fourier_basis(lower = 0, upper = 1, dimension = 9)
  # deriv_operator(m) through the operator route must equal order = m
  for (m in 1:3) {
    expect_equal(basis_gram(b, order = deriv_operator(m)),
                 basis_gram(b, order = m))
  }
  # on a spline both routes integrate knot interval by knot interval and
  # both are exact, so they agree to rounding
  bs <- bspline_basis(lower = 0, upper = 1, dimension = 10, degree = 3)
  for (m in 1:3) {
    expect_equal(basis_gram(bs, order = deriv_operator(m)),
                 basis_gram(bs, order = m), tolerance = 1e-12)
  }
  # and the base method's evenly spaced panels, which do not line up with
  # the knots, are the looser route the family overrides
  expect_gt(max(abs(basis7:::numerical_operator_gram(bs, deriv_operator(2)) -
                    basis_gram(bs, order = 2L))) /
              max(abs(basis_gram(bs, order = 2L))), 1e-9)
})

test_that("the Fourier roughness matrix of an operator is exactly diagonal", {
  # ⚠️ The book states that the harmonic penalty makes R "structurally
  # different and more complex" than the diagonal matrix of the derivative
  # penalty. It is structurally different -- the null space grows from one
  # to three -- and it is not more complex: it is exactly as diagonal.
  b <- fourier_basis(lower = 0, upper = 1, dimension = 13)
  for (op in list(harmonic_operator(1), harmonic_operator(1, harmonics = 2),
                  oscillator_operator(1), deriv_operator(2) * oscillator_operator(1))) {
    g <- basis_gram(b, order = op)
    expect_lt(max(abs(g - diag(diag(g)))) / max(abs(g)), 1e-12)
  }
})

test_that("the Fourier closed form is the modulus of the characteristic polynomial", {
  # |P(i nu_j)|^2 * T/2, transcribed here by hand rather than read from the
  # package, so the two are independent statements of the same formula
  period <- 365
  k <- 11L
  b <- fourier_basis(lower = 0, upper = period, dimension = k)
  op <- harmonic_operator(period)
  nu <- 2 * pi / period
  j <- seq_len((k - 1L) / 2L)
  want <- c(0, rep(j^2 * (j^2 - 1)^2 * nu^6 * period / 2, each = 2L))
  expect_equal(diag(basis_gram(b, order = op)), want, ignore_attr = TRUE)

  # and it agrees with integrating (Lb)(Lb)' numerically, which shares no
  # arithmetic with it
  gn <- basis7:::numerical_operator_gram(b, op, panels = 300L)
  expect_lt(max(abs(basis_gram(b, order = op) - gn)) /
              max(abs(basis_gram(b, order = op))), 1e-10)
})

test_that("the closed form is not taken where it does not hold", {
  # a period that is not the interval width gives up the orthogonality
  f <- fourier_basis(lower = 0, upper = 1, dimension = 7, omega = 0.7)
  g <- basis_gram(f, order = harmonic_operator(1))
  expect_gt(max(abs(g - diag(diag(g)))) / max(abs(g)), 1e-6)
  # and the empirical and weighted measures go through the quadrature route
  x <- sort(runif(200))
  b <- fourier_basis(lower = 0, upper = 1, dimension = 7)
  ge <- basis_gram(b, order = harmonic_operator(1), at = x)
  expect_equal(dim(ge), c(7L, 7L))
  expect_lt(max(abs(ge[, 1:3])), 1e-9)
  gw <- basis_gram(b, order = harmonic_operator(1),
                   weight = function(t) dnorm(t, 0.5, 0.25))
  expect_equal(dim(gw), c(7L, 7L))
})

test_that("an operator above a spline's degree is refused, not truncated", {
  b <- bspline_basis(lower = 0, upper = 1, dimension = 12, degree = 3)
  expect_error(basis_gram(b, order = harmonic_operator(1, harmonics = 2)),
               "leading\n  term deleted")
  # at the degree it is admitted
  expect_equal(dim(basis_gram(b, order = harmonic_operator(1))), c(12L, 12L))
})

test_that("an operator is refused on a basis of several variables", {
  tb <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
  expect_error(basis_gram(tb, order = deriv_operator(2)), "several")
})

test_that("the null space of the operator is not read off the matrix", {
  # ⚠️ The two questions give different answers, which is why operator_null()
  # is computed from the characteristic polynomial. On a Fourier basis the
  # rank deficiency is the operator's own; on a cubic spline it is not,
  # because a spline represents a sine only approximately.
  op <- harmonic_operator(1)
  rank_def <- function(g) {
    ev <- eigen(g, symmetric = TRUE)$values
    ncol(g) - sum(ev > 1e-10 * max(ev))
  }
  gf <- basis_gram(fourier_basis(0, 1, 13), order = op)
  gb <- basis_gram(bspline_basis(0, 1, 20, 3), order = op)
  expect_identical(rank_def(gf), 3L)
  expect_lt(rank_def(gb), 3L)
  # the operator says three either way
  expect_identical(nrow(operator_null(op)), 3L)
})
