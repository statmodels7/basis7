test_that("a cyclic smoother has k functions and a block one narrower", {
  set.seed(11)
  x <- sort(runif(200))
  sm <- cyclic_smooth(k = 10, lower = 0, upper = 1)
  expect_true(S7::S7_inherits(sm, smoother))
  expect_true(S7::S7_inherits(sm, CyclicSmoother))
  expect_identical(sm@dimension, 10L)
  expect_identical(smoother_basis(sm, x)@dimension, 10L)

  out <- smoother_build(sm, x)
  # the constant is removed and nothing is restored, so k - 1 columns
  expect_identical(ncol(out$X), 9L)
  expect_identical(dim(out$S), c(9L, 9L))
  expect_identical(out$unpenalized, 0L)

  # CONTROL: the non-periodic family of the same k restores the linear
  # column its order-2 penalty leaves free, so it does have one
  ob <- smoother_build(bspline_smooth(k = 10, lower = 0, upper = 1), x)
  expect_identical(ob$unpenalized, 1L)
})

test_that("the block and a fit on it are periodic in every derivative", {
  set.seed(12)
  x <- sort(runif(300, 0, 365))
  sm <- cyclic_smooth(k = 12, degree = 3, lower = 0, upper = 365)
  b <- smoother_basis(sm, x)

  # the BASIS matches at the two ends in value and in the first degree - 1
  # derivatives, which is the definition the construction imposes
  for (j in 0:2) {
    m <- basis_deriv(b, c(0, 365), order = j)
    expect_lt(max(abs(m[1L, ] - m[2L, ])), 1e-8 * max(1, max(abs(m))))
  }

  # and so does the BLOCK the reparametrization produces
  out <- smoother_build(sm, x)
  ends <- smoother_apply(sm, out$blueprint, c(0, 365))
  expect_lt(max(abs(ends[1L, ] - ends[2L, ])), 1e-10)

  # CONTROL: an ordinary B-spline of the same k does not, by a margin of
  # the size of the block itself
  nb <- bspline_smooth(k = 12, lower = 0, upper = 365)
  ob <- smoother_build(nb, x)
  oe <- smoother_apply(nb, ob$blueprint, c(0, 365))
  expect_gt(max(abs(oe[1L, ] - oe[2L, ])), 0.1)
})

test_that("the roughness matrix integrates over ONE period", {
  # The reference shares no code with basis_gram(): a trapezoid of the
  # basis's own second derivatives, on a grid aligned to the knots so that
  # the rule is O(h^2) rather than erratic -- the second derivative of a
  # cubic spline is piecewise linear with kinks at the knots.
  k <- 9L
  sm <- cyclic_smooth(k = k, degree = 3, order = 2, lower = 0, upper = 1)
  b <- smoother_basis(sm, sort(runif(50)))
  g <- basis_gram(b, order = 2L)
  scale <- max(abs(g))
  expect_gt(scale, 0)

  trap <- function(lo, hi, n) {
    xs <- seq(lo, hi, length.out = n)
    d2 <- basis_deriv(b, xs, order = 2L)
    w <- rep(1, n); w[1L] <- w[n] <- 0.5
    crossprod(d2 * (w * (hi - lo) / (n - 1L)), d2)
  }
  n <- as.integer(k) * 2000L + 1L
  expect_lt(max(abs(g - trap(0, 1, n))), 1e-5 * scale)

  # CONTROL: the same integral over 99 per cent of the period is a
  # different matrix, by four orders more than the tolerance above
  expect_gt(max(abs(g - trap(0, 0.99, n))), 1e-2 * scale)
})

test_that("the null space of the pair is the constant at every order", {
  # A non-constant periodic function is never a polynomial, so unlike a
  # B-spline's the null space does not grow with the order.
  x <- sort(runif(80))
  for (ord in 1:3) {
    sm <- cyclic_smooth(k = 9, degree = 3, order = ord, lower = 0, upper = 1)
    g <- basis_gram(smoother_basis(sm, x), order = ord)
    ev <- eigen(g, symmetric = TRUE)
    expect_identical(sum(ev$values <= 1e-8 * max(ev$values)), 1L)
    # and the direction it leaves free really is a constant function
    v <- ev$vectors[, which.min(ev$values)]
    f <- basis_eval(smoother_basis(sm, x), seq(0, 1, length.out = 25L)) %*% v
    expect_lt(max(abs(f - mean(f))), 1e-8)
  }

  # CONTROL: the non-periodic family of the same degree grows its null
  # space with the order, which is what makes the reading above a claim
  nb <- bspline_smooth(k = 9, degree = 3, order = 3, lower = 0, upper = 1)
  gb <- basis_gram(smoother_basis(nb, x), order = 3L)
  eb <- eigen(gb, symmetric = TRUE, only.values = TRUE)$values
  expect_identical(sum(eb <= 1e-8 * max(eb)), 3L)
})

test_that("a cyclic block is reapplied and not rebuilt", {
  set.seed(13)
  x <- sort(runif(250, -2, 5))
  sm <- cyclic_smooth(k = 11, lower = -2, upper = 5)
  out <- smoother_build(sm, x)
  rows <- c(1L, 7L, 60L, 199L, 250L)
  expect_reapplied(smoother_apply(sm, out$blueprint, x[rows]),
                   out$X[rows, , drop = FALSE])
})

test_that("the periodicity constraint has one independent row per degree", {
  for (d in 1:4) {
    p <- bspline_basis(0, 1, dimension = 8L + d, degree = d)
    cm <- periodic_constraint(p, d)
    expect_identical(dim(cm), c(as.integer(d), 8L + as.integer(d)))
    s <- svd(cm, nu = 0L, nv = 0L)
    expect_identical(sum(s$d > 1e-10 * max(s$d)), as.integer(d))
  }
})

test_that("a cyclic smoother declares only what means something to it", {
  # 'order' cannot exceed the degree: the derivative of that order of a
  # spline of that degree is zero and the penalty would be the zero matrix
  expect_error(cyclic_smooth(k = 10, degree = 3, order = 4), "exceeds 'degree'")
  # the constant is removed, so k = 2 leaves one column and k = 1 none
  expect_error(cyclic_smooth(k = 2), "at least 3")
  # 'constrain' has no reading on a periodic basis, so it is not an argument
  expect_error(cyclic_smooth(k = 10, constrain = 1), "unused argument")
  # and the contradictory pair the smoother class refuses
  expect_error(
    cyclic_smooth(k = 10, penalty = penalties7::lasso_penalty,
                  null_space = "shrink"),
    "cannot both be given"
  )
})

test_that("a penalty factory is stored on a cyclic smoother and never called", {
  called <- 0L
  fac <- function(n_coef) {
    called <<- called + 1L
    penalties7::lasso_penalty(n_coef = n_coef)
  }
  sm <- cyclic_smooth(k = 10, penalty = fac, lower = 0, upper = 1)
  expect_identical(called, 0L)
  # building the block does not call it either: basis7 stores the factory
  # and modelterms7 is what calls it
  out <- smoother_build(sm, sort(runif(60)))
  expect_identical(called, 0L)
  expect_identical(ncol(out$X), 9L)
})
