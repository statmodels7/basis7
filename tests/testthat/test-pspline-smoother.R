test_that("a P-spline is a B-spline basis with a difference penalty", {
  set.seed(31)
  x <- sort(runif(300))
  sm <- pspline_smooth(k = 20, degree = 3, diff = 2, lower = 0, upper = 1)
  expect_true(S7::S7_inherits(sm, smoother))
  expect_true(S7::S7_inherits(sm, PsplineSmoother))
  expect_identical(sm@dimension, 20L)
  expect_identical(sm@diff, 2L)

  # the basis is an ordinary B-spline: same object bspline_smooth() builds
  b <- smoother_basis(sm, x)
  bb <- smoother_basis(bspline_smooth(k = 20, degree = 3, lower = 0, upper = 1), x)
  expect_identical(basis_eval(b, x), basis_eval(bb, x))

  # and the roughness matrix is NOT its Gram matrix
  expect_false(isTRUE(all.equal(smoother_gram(sm, b, x),
                                basis_gram(b, order = 2L),
                                check.attributes = FALSE)))
})

test_that("the roughness matrix is the difference operator on the Eilers-Marx coefficients", {
  # The Eilers-Marx basis is built here with splines::splineDesign on knots
  # equally spaced beyond the interval, a route that shares nothing with the
  # package's clamped basis. The clamped coefficients c of a function are
  # carried to its Eilers-Marx coefficients a = M c, and the penalty on c is
  # the difference penalty on a.
  x <- sort(runif(120))
  xs <- seq(0, 1, length.out = 401)
  for (dd in 1:3) {
    k <- 20L
    sm <- pspline_smooth(k = k, degree = 3, diff = dd, lower = 0, upper = 1)
    b <- smoother_basis(sm, x)
    g <- smoother_gram(sm, b, x)
    expect_identical(dim(g), c(k, k))
    h <- 1 / (k - 3)
    be <- splines::splineDesign(h * seq(-3, k), xs, ord = 4, outer.ok = TRUE)
    M <- qr.solve(be, as.matrix(basis_eval(b, xs)))
    D <- base::diff(diag(k), differences = dd)
    expect_lt(max(abs(g - crossprod(D %*% M))), 1e-10 * max(abs(g)))
    ev <- eigen(g, symmetric = TRUE, only.values = TRUE)$values
    expect_identical(sum(ev <= 1e-8 * max(ev)), as.integer(dd))
  }
})

test_that("the null space of a P-spline is exactly the polynomials", {
  # On the Eilers-Marx knots the Greville abscissae are equally spaced, so
  # coefficients polynomial in the index give a polynomial in x. On clamped
  # coefficients (basis7 <= 0.13.1) the R^2 below was 0.9994 at diff = 2.
  x <- seq(0, 1, length.out = 300)
  for (dd in 2:3) {
    sm <- pspline_smooth(k = 20, diff = dd, lower = 0, upper = 1)
    b <- smoother_basis(sm, x)
    e <- eigen(smoother_gram(sm, b, x), symmetric = TRUE)
    N <- e$vectors[, e$values <= 1e-8 * max(e$values), drop = FALSE]
    fn <- as.matrix(basis_eval(b, x)) %*% N
    Q <- outer(x, seq.int(0L, dd - 1L), "^")
    res <- fn - Q %*% qr.solve(Q, fn)
    expect_lt(max(abs(res)), 1e-10 * max(abs(fn)))
  }
})

test_that("the built penalty is the identity and the fit contracts to a line", {
  set.seed(32)
  x <- sort(runif(300))
  sm <- pspline_smooth(k = 25, degree = 3, diff = 2, lower = 0, upper = 1)
  out <- smoother_build(sm, x)
  # smoother_span() constrains against the exact polynomials, so what the
  # Demmler-Reinsch rotation sees is positive definite
  ev <- eigen(out$S, symmetric = TRUE, only.values = TRUE)$values
  pen <- ev[ev > 1e-10]
  expect_length(pen, ncol(out$S) - out$unpenalized)
  expect_lt(max(abs(pen - 1)), 1e-8)

  # and a strongly penalized fit really is a straight line
  y <- sin(3 * pi * x) + 0.6 * x + rnorm(length(x), sd = 0.3)
  X <- cbind(1, out$X)
  S <- rbind(0, cbind(0, out$S))
  f <- X %*% solve(crossprod(X) + 1e10 * S, crossprod(X, y))
  expect_gt(summary(stats::lm(f ~ x))$r.squared, 1 - 1e-8)

  # CONTROL: at a light penalty it is not a straight line, so the check
  # above cannot be satisfied by a construction that fits nothing
  fl <- X %*% solve(crossprod(X) + 1e-6 * S, crossprod(X, y))
  expect_lt(summary(stats::lm(fl ~ x))$r.squared, 0.9)
})

test_that("a P-spline block is reapplied and not rebuilt", {
  set.seed(33)
  x <- sort(runif(200, -1, 3))
  sm <- pspline_smooth(k = 18, lower = -1, upper = 3)
  out <- smoother_build(sm, x)
  rows <- c(1L, 20L, 95L, 200L)
  expect_reapplied(smoother_apply(sm, out$blueprint, x[rows]),
                   out$X[rows, , drop = FALSE])
})

test_that("a P-spline declares only what means something to it", {
  # nothing is integrated, so there is no measure to integrate against and
  # no derivative whose order to name
  expect_error(pspline_smooth(k = 20, measure = "empirical"), "unused argument")
  expect_error(pspline_smooth(k = 20, order = 2), "unused argument")
  # the d-th difference of k coefficients needs k > d
  expect_error(pspline_smooth(k = 4, degree = 3, diff = 4), "smaller than 'k'")
  expect_error(pspline_smooth(k = 3, degree = 3), "too small for 'degree'")
  expect_error(pspline_smooth(k = 20, diff = 0), "at least 1")
})

test_that("a penalty factory is stored on a P-spline and never called", {
  called <- 0L
  fac <- function(n_coef) {
    called <<- called + 1L
    stop("the smoother must not call this")
  }
  sm <- pspline_smooth(k = 15, penalty = fac, lower = 0, upper = 1)
  expect_identical(sm@penalty, fac)
  x <- sort(stats::runif(90))
  out <- smoother_build(sm, x)
  expect_identical(called, 0L)
  ref <- smoother_build(pspline_smooth(k = 15, lower = 0, upper = 1), x)
  expect_identical(out$X, ref$X)
  expect_identical(out$S, ref$S)
})
