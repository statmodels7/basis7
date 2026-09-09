# The smoother contract. What is asserted here is what the construction
# PROMISES, not that it agrees with modelterms7::s(): basis7 sits at the
# bottom of the dependency graph and may not name a package above it, even in
# Suggests. That comparison -- the block and the penalty identical bit for
# bit -- is run against modelterms7 when s() is moved onto this contract.

set.seed(11)
n <- 240L
x <- sort(runif(n, -2, 3))

test_that("a smoother is a recipe carrying the four decisions", {
  sm <- bspline_smooth(k = 10)
  expect_true(S7::S7_inherits(sm, smoother))
  expect_true(S7::S7_inherits(sm, BsplineSmoother))
  expect_identical(sm@dimension, 10L)
  expect_identical(sm@degree, 3L)
  expect_identical(sm@order, 2L)
  expect_identical(sm@null_space, "keep")
  expect_identical(sm@reparam, "dr")
  expect_null(sm@lower)
  expect_null(sm@upper)
  # abstract: there is no direct constructor
  expect_error(smoother(), "abstract")
})

test_that("the constructor rejects what it cannot build", {
  # how many directions the constraint removes depends on 'order' and
  # 'constrain', so the floor on 'k' is not a constant: it is reported
  # against what the constraint actually takes
  expect_error(bspline_smooth(k = 1), "at least 2")
  # a basis too small for its own degree is reported as that, first
  expect_error(bspline_smooth(k = 2), "too small for 'degree'")
  # THE CONSTRAINT LEAVING NOTHING is reachable through 'constrain' alone:
  # with the default the arithmetic forbids it, k being at least degree + 1
  # and order at most degree, so k always exceeds the order.
  expect_error(bspline_smooth(k = 3, degree = 2, order = 2, constrain = 3),
               "leaves nothing to smooth")
  expect_error(bspline_smooth(k = 3, degree = 3), "too small for 'degree'")
  expect_error(bspline_smooth(k = 10, degree = 0), "at least 1")
  expect_error(bspline_smooth(k = 10.5), "whole number")
  expect_error(bspline_smooth(lower = 1, upper = 0), "strictly less")
  expect_error(bspline_smooth(lower = c(0, 1)), "single finite number")
  # ORDER MAY NOT EXCEED DEGREE: above it the roughness matrix is identically
  # zero, so the penalty would penalize nothing and every direction would be
  # free. Reported where the two numbers were written.
  expect_error(bspline_smooth(k = 10, degree = 3, order = 4),
               "exceeds 'degree'")
})

test_that("the block has one column per coordinate and the sizes follow k", {
  for (k in c(5L, 10L, 20L)) {
    keep <- smoother_build(bspline_smooth(k = k), x)
    drop <- smoother_build(bspline_smooth(k = k, null_space = "drop"), x)
    # the constraint against the constant and the linear function removes
    # two directions; keeping the null space adds one back
    expect_identical(ncol(keep$X), k - 1L)
    expect_identical(ncol(drop$X), k - 2L)
    expect_identical(length(keep$names), ncol(keep$X))
    expect_identical(nrow(keep$X), n)
    expect_null(dimnames(keep$X))
  }
})

test_that("the penalty is the identity off the free columns", {
  out <- smoother_build(bspline_smooth(k = 10), x)
  expect_identical(out$S, diag(c(0, rep(1, 8)), 9L))
  expect_identical(out$unpenalized, 1L)
  expect_identical(out$names[[1L]], "lin")

  dr <- smoother_build(bspline_smooth(k = 10, null_space = "drop"), x)
  expect_identical(dr$S, diag(rep(1, 8), 8L))
  expect_identical(dr$unpenalized, 0L)
  expect_identical(dr$names[[1L]], "z1")
})

test_that("the free column is the covariate, and the rest is orthogonal to it", {
  out <- smoother_build(bspline_smooth(k = 10), x)
  # it IS the standardized covariate, so it carries the linear effect
  expect_equal(cor(out$X[, 1L], x), 1, tolerance = 1e-12)
  expect_equal(mean(out$X[, 1L]), 0, tolerance = 1e-12)
  expect_equal(sd(out$X[, 1L]), 1, tolerance = 1e-12)
  # and the constraint of dr_basis makes the smooth part orthogonal to it
  # AND to the constant, over the observed covariate
  expect_lt(max(abs(crossprod(out$X[, 1L], out$X[, -1L]))), 1e-9)
  expect_lt(max(abs(colSums(out$X[, -1L, drop = FALSE]))), 1e-9)
})

test_that("the Demmler-Reinsch columns are ordered from smooth to rough", {
  out <- smoother_build(bspline_smooth(k = 12), x)
  z <- out$X[, -1L, drop = FALSE]
  # orthogonal over the data, which is what the pencil buys
  g <- crossprod(z)
  expect_lt(max(abs(g - diag(diag(g)))), 1e-9)
  # and increasingly oscillatory: the number of sign changes does not fall
  crossings <- apply(z, 2L, function(v) sum(diff(sign(v)) != 0))
  expect_false(is.unsorted(crossings))
})

test_that("smoother_apply reapplies the transform and does not rebuild it", {
  sm <- bspline_smooth(k = 10)
  out <- smoother_build(sm, x)
  i <- c(3L, 17L, 40L, 91L, 150L, 233L)

  # REAPPLIED, the rows are the rows of the block, exactly
  expect_identical(smoother_apply(sm, out$blueprint, x[i]), out$X[i, ])

  # REBUILT on the same rows, they are NOT: the interval is the range of
  # the subset and the rotation diagonalizes a different Gram matrix, so a
  # rebuild is a different function of the covariate. This is the reason
  # the blueprint exists.
  rebuilt <- smoother_build(sm, x[i])$X
  expect_gt(max(abs(rebuilt - out$X[i, ])), 1e-3)
})

test_that("a fixed interval is used, and an unfixed one is read from the data", {
  b <- smoother_basis(bspline_smooth(k = 8, lower = -5, upper = 5), x)
  expect_identical(c(b@lower, b@upper), c(-5, 5))

  bd <- smoother_basis(bspline_smooth(k = 8), x)
  # the observed values lie strictly inside, which is what the padding is for
  expect_lt(bd@lower, min(x))
  expect_gt(bd@upper, max(x))
  expect_equal(bd@lower, min(x) - diff(range(x)) * 0.001, tolerance = 1e-12)

  # a fixed interval is what a prediction outside the observed range needs
  sm <- bspline_smooth(k = 8, lower = -5, upper = 5)
  out <- smoother_build(sm, x)
  expect_silent(smoother_apply(sm, out$blueprint, c(-4.5, 4.5)))
})

test_that("the covariate is checked", {
  sm <- bspline_smooth(k = 8)
  expect_error(smoother_build(sm, numeric(0)), "at least one value")
  expect_error(smoother_build(sm, c(1, NA, 2)), "missing values")
  expect_error(smoother_build(sm, c(1, Inf, 2)), "finite")
  out <- smoother_build(sm, x)
  expect_error(smoother_apply(sm, list(), x), "blueprint element")
})

test_that("a setting this version does not build is rejected, not ignored", {
  # An argument accepted and ignored would report a fit of a model the
  # caller did not ask for. What remains unbuilt is the penalty factory,
  # which is a lot of its own; `order`, `measure`, `constrain`,
  # `null_space` and `reparam` are built and are tested in their own files.
  expect_error(bspline_smooth(k = 10, penalty = function(n_coef) NULL),
               "not built in this version")
  expect_silent(bspline_smooth(k = 10, null_space = "shrink"))
  expect_silent(bspline_smooth(k = 10, reparam = "none"))
  expect_silent(bspline_smooth(k = 10, order = 1))
  # and an unknown value is rejected by name before that
  expect_error(bspline_smooth(k = 10, null_space = "nope"), "arg")
  expect_error(bspline_smooth(k = 10, reparam = "nope"), "arg")
})

test_that("(X, S) is what a penalized least-squares fit needs", {
  # the smoother's whole contract, exercised: solve(X'X + lambda S, X'y)
  #
  # k IS PART OF THE CLAIM, and it was measured rather than assumed: a basis
  # that cannot overfit leaves the penalty nothing to remove at the rough
  # end. At k = 15, fourteen columns over 240 points, the unpenalized fit
  # already reaches 0.1907 against the best 0.1877, so it is not beaten. At
  # k = 60 it reaches 0.2064 and is.
  out <- smoother_build(bspline_smooth(k = 60), x)
  y <- sin(2 * x) + rnorm(n, sd = 0.2)
  rmse <- vapply(c(1e-6, 0.1, 1e6), function(l) {
    b <- solve(crossprod(out$X) + l * out$S, crossprod(out$X, y))
    sqrt(mean((out$X %*% b - sin(2 * x))^2))
  }, numeric(1))
  # the optimum is interior: the penalty helps against a rough fit and hurts
  # when it is strong enough to leave only a line
  expect_lt(rmse[2L], rmse[1L])
  expect_lt(rmse[2L], rmse[3L])

  # AT A LARGE SMOOTHING PARAMETER THE FIT IS A STRAIGHT LINE, which is what
  # order = 2 means: the penalty leaves the linear column alone.
  b <- solve(crossprod(out$X) + 1e10 * out$S, crossprod(out$X, y))
  expect_gt(summary(lm(out$X %*% b ~ x))$r.squared, 1 - 1e-9)
})

test_that("a smoother prints its four decisions", {
  out <- capture.output(print(bspline_smooth(k = 12)))
  expect_match(out[1L], "BsplineSmoother")
  expect_match(paste(out, collapse = " "), "12 functions")
  expect_match(paste(out, collapse = " "), "order 2")
  expect_match(paste(out, collapse = " "), "keep")
  expect_match(paste(out, collapse = " "), "from the data")
  expect_match(paste(capture.output(print(bspline_smooth(k = 5, lower = 0,
                                                         upper = 1))),
                     collapse = " "), "\\[0, 1\\]")
})
