# The last two axes: which coordinates the coefficients live in, and what
# becomes of the directions the penalty does not see.

set.seed(6)
n <- 300L
x <- sort(runif(n, -2, 2))
y <- sin(3 * x) + rnorm(n, sd = 0.3)

test_that("the three coordinate systems describe the same space", {
  o <- lapply(c("dr", "none", "orthonorm"), function(r) {
    smoother_build(bspline_smooth(k = 12, reparam = r), x)
  })
  names(o) <- c("dr", "none", "orthonorm")

  # same number of coordinates, and the same span: an UNPENALIZED fit
  # cannot tell them apart, which is what says a reparametrization changes
  # what a coefficient means and not what the term can express
  expect_identical(vapply(o, function(e) ncol(e$X), integer(1)),
                   c(dr = 11L, none = 11L, orthonorm = 11L))
  f <- lapply(o, function(e) as.vector(e$X %*% qr.solve(e$X, y)))
  expect_lt(max(abs(f$dr - f$none)), 1e-10)
  expect_lt(max(abs(f$dr - f$orthonorm)), 1e-10)
})

test_that("each coordinate system has the property it is named for", {
  dr <- smoother_build(bspline_smooth(k = 12, reparam = "dr"), x)
  # DR: the penalty is the identity on what it covers
  pen <- dr$S[-1L, -1L]
  expect_identical(pen, diag(1, nrow(pen)))

  # NONE: the penalty is the congruence of the roughness matrix, which is
  # not diagonal -- that is the whole difference
  nn <- smoother_build(bspline_smooth(k = 12, reparam = "none"), x)
  pn <- nn$S[-1L, -1L]
  expect_gt(max(abs(pn - diag(diag(pn)))), 1e-6)
  expect_true(isSymmetric(nn$S))

  # ORTHONORM: X'X = I for the block, with the free column dropped. With it
  # kept, a column is prepended afterwards and the block is not orthonormal.
  on <- smoother_build(
    bspline_smooth(k = 12, reparam = "orthonorm", null_space = "drop"), x)
  expect_lt(max(abs(crossprod(on$X) - diag(ncol(on$X)))), 1e-10)
  onk <- smoother_build(bspline_smooth(k = 12, reparam = "orthonorm"), x)
  expect_gt(max(abs(crossprod(onk$X) - diag(ncol(onk$X)))), 1e-6)
})

test_that("every coordinate system is reapplied, not rebuilt", {
  # AN IDENTITY, and it is structural: the block is evaluated from the same
  # basis object smoother_apply() evaluates. Built as a local product
  # instead, the orthonormal route differed in the last bit, new_transformed()
  # flattening a nested transform so that B (T1 T2) is computed where
  # (B T1) T2 was.
  i <- c(3L, 60L, 250L)
  for (r in c("dr", "none", "orthonorm")) {
    sm <- bspline_smooth(k = 12, reparam = r)
    o <- smoother_build(sm, x)
    expect_reapplied(smoother_apply(sm, o$blueprint, x[i]), o$X[i, ])
  }
})

test_that("shrink penalizes the null space so the term can leave the model", {
  o <- lapply(c("keep", "drop", "shrink"), function(ns) {
    smoother_build(bspline_smooth(k = 12, null_space = ns), x)
  })
  names(o) <- c("keep", "drop", "shrink")

  # the weight is mgcv's rule translated: a tenth of what a penalized
  # direction carries, which in these coordinates is exactly 1
  expect_identical(diag(o$keep$S)[[1L]], 0)
  expect_identical(diag(o$shrink$S)[[1L]], 0.1)
  expect_identical(diag(o$drop$S)[[1L]], 1)

  # shrink keeps the column but leaves nothing unpenalized
  expect_identical(ncol(o$shrink$X), ncol(o$keep$X))
  expect_identical(o$shrink$unpenalized, 0L)
  expect_identical(o$keep$unpenalized, 1L)

  # AND THE CONSEQUENCE: at a large smoothing parameter a kept null space
  # survives as a straight line and a shrunk one does not
  sdfit <- vapply(o, function(e) {
    sd(e$X %*% solve(crossprod(e$X) + 1e10 * e$S, crossprod(e$X, y)))
  }, numeric(1))
  expect_gt(sdfit[["keep"]], 1e-3)
  expect_lt(sdfit[["shrink"]], 1e-6)
  expect_lt(sdfit[["drop"]], 1e-6)
})

test_that("shrink is a no-op where there is no null space to shrink", {
  # a Fourier smoother removes the constant and restores nothing, so its
  # block is entirely penalized whatever null_space says
  xs <- sort(runif(200))
  a <- smoother_build(fourier_smooth(k = 9, lower = 0, upper = 1), xs)
  b <- smoother_build(
    fourier_smooth(k = 9, null_space = "shrink", lower = 0, upper = 1), xs)
  expect_identical(a$X, b$X)
  expect_identical(a$S, b$S)
})

test_that("the axes compose", {
  # order, constrain, measure, null_space and reparam on one object
  sm <- bspline_smooth(k = 16, degree = 5, order = 3, constrain = 3,
                       measure = "empirical", null_space = "shrink",
                       reparam = "none")
  o <- smoother_build(sm, x)
  # 16 functions less 4 constrained, plus the two restored
  expect_identical(ncol(o$X), 14L)
  expect_identical(o$unpenalized, 0L)
  expect_identical(diag(o$S)[1:2], c(0.1, 0.1))
  expect_reapplied(smoother_apply(sm, o$blueprint, x[1:5]), o$X[1:5, ])
})
