# The axes of a smoother -- order, measure, constrain -- and the two
# families whose bases the package already ships. Every assertion here is a
# statement about the construction that can be read off the result; the
# comparison with modelterms7 is the identity gate, which lives outside the
# package because basis7 may not name a package above it.

pfit <- function(X, S, y, lambda) {
  solve(crossprod(X) + lambda * S, crossprod(X, y))
}

test_that("order says what a strongly penalized fit contracts toward", {
  set.seed(7)
  x <- sort(runif(400, -2, 2))
  y <- sin(3 * x) + rnorm(400, sd = 0.3)

  # ORDER 1: a constant. The residual variation vanishes as 1/lambda, so
  # what is asserted is that it is gone, not an R-squared -- the ratio of
  # two null quantities is what the plan records reading wrongly once.
  o1 <- smoother_build(bspline_smooth(k = 20, degree = 5, order = 1), x)
  f1 <- as.vector(o1$X %*% pfit(o1$X, o1$S, y, 1e12))
  expect_lt(sd(f1), 1e-6)
  expect_identical(o1$unpenalized, 0L)

  # ORDER 2: a straight line, and one free column
  o2 <- smoother_build(bspline_smooth(k = 20, degree = 5, order = 2), x)
  f2 <- as.vector(o2$X %*% pfit(o2$X, o2$S, y, 1e12))
  expect_gt(sd(f2), 1e-3)
  expect_gt(summary(stats::lm(f2 ~ x))$r.squared, 1 - 1e-6)
  expect_identical(o2$unpenalized, 1L)

  # ORDER 3: a parabola, which a line does NOT explain, and two free columns
  o3 <- smoother_build(bspline_smooth(k = 20, degree = 5, order = 3), x)
  f3 <- as.vector(o3$X %*% pfit(o3$X, o3$S, y, 1e12))
  expect_gt(summary(stats::lm(f3 ~ poly(x, 2)))$r.squared, 1 - 1e-6)
  expect_lt(summary(stats::lm(f3 ~ x))$r.squared, 0.99)
  expect_identical(o3$unpenalized, 2L)
  expect_identical(o3$names[1:2], c("lin", "poly2"))
})

test_that("the free columns are orthogonal to each other and to the block", {
  set.seed(4)
  x <- sort(runif(300, -1, 4))
  o <- smoother_build(bspline_smooth(k = 15, degree = 5, order = 3), x)
  fr <- o$X[, 1:2, drop = FALSE]
  # standardized, mutually orthogonal, and orthogonal to the penalized part
  expect_equal(colMeans(fr), c(0, 0), tolerance = 1e-10)
  expect_equal(apply(fr, 2L, sd), c(1, 1), tolerance = 1e-10)
  expect_lt(abs(crossprod(fr[, 1L], fr[, 2L])), 1e-8)
  expect_lt(max(abs(crossprod(fr, o$X[, -(1:2)]))), 1e-8)
  # the first is the linear one
  expect_equal(abs(cor(fr[, 1L], x)), 1, tolerance = 1e-12)
  # AND THE PAIR SPANS the null space with the constant, which is the
  # property "keep restores the null space minus the constant" means. The
  # second column on its own correlates only 0.38 with x^2 -- it is the
  # part of x^2 the constant and x do not already explain -- so a
  # correlation is the wrong reading and the span is the right one.
  # read as a residual rather than an R-squared: x is exactly the first
  # column up to scaling, so lm()'s summary warns that the fit is perfect
  expect_lt(max(abs(residuals(stats::lm(x ~ fr)))), 1e-10)
  expect_lt(max(abs(residuals(stats::lm(x^2 ~ fr)))), 1e-10)
})

test_that("the free columns are reapplied and not rebuilt", {
  set.seed(4)
  x <- sort(runif(300, -1, 4))
  sm <- bspline_smooth(k = 15, degree = 5, order = 3)
  o <- smoother_build(sm, x)
  i <- c(5L, 44L, 200L)
  expect_reapplied(smoother_apply(sm, o$blueprint, x[i]), o$X[i, ])
})

test_that("the reapplied check refuses a rebuilt basis", {
  # THE POSITIVE CONTROL for expect_reapplied(). A tolerance nothing can
  # fail is not a check, so the defect the helper exists for is put to it:
  # a basis REBUILT from the new points carries different knots, a
  # different empirical Gram and a different rotation, and must be
  # rejected. The other edge is asserted beside it -- a perturbation of the
  # size the platforms were measured to differ by must still pass, so the
  # tolerance cannot be tightened back into redness without a test saying
  # so.
  set.seed(4)
  x <- sort(runif(300, -1, 4))
  sm <- bspline_smooth(k = 15, degree = 5, order = 3)
  o <- smoother_build(sm, x)
  i <- c(5L, 44L, 200L)
  ref <- o$X[i, ]

  rebuilt <- smoother_build(sm, x[i])$X
  expect_identical(dim(rebuilt), dim(ref))
  # of the size of the quantity itself, not of its last bits
  expect_gt(max(abs(rebuilt - ref)) / max(abs(ref)), 0.1)
  # and above the helper own threshold, READ FROM THE DEFAULT rather than
  # copied, so the control cannot drift from the rule it controls
  tol <- eval(formals(expect_reapplied)[["tolerance"]])
  expect_gt(max(abs(rebuilt - ref)), tol * max(abs(ref)))

  expect_reapplied(ref * (1 + 2 * .Machine$double.eps), ref)
})

test_that("smoother_apply reaches no part of the construction", {
  # THE ROUTING HALF, WHICH A TOLERANCE CANNOT CARRY. The test above says
  # the numbers agree; this says the route is the stored blueprint, and a
  # rebuild that happened to land close would pass the first and must not
  # pass this. Both halves together are what the identity used to assert
  # before a BLAS's shape-dependence made the identity unavailable.
  set.seed(4)
  x <- sort(runif(300, -1, 4))
  sm <- bspline_smooth(k = 15, degree = 5, order = 3)
  o <- smoother_build(sm, x)

  # the two steps the CONSTRUCTION goes through and the reapplication must
  # not: the empirical Gram, which depends on every observation, and the
  # reparametrization built from it
  local_mocked_bindings(
    smoother_gram = function(...) stop("reached the Gram", call. = FALSE),
    smoother_reparam = function(...) stop("reached the reparam", call. = FALSE)
  )
  expect_no_error(smoother_apply(sm, o$blueprint, x[c(5L, 44L, 200L)]))

  # AND THE MOCK IS PROVEN LIVE, or the assertion above passes vacuously
  # whether or not the bindings were ever replaced
  expect_error(smoother_build(sm, x), "reached the")
})

test_that("constrain must contain the null space and may exceed it", {
  set.seed(9)
  x <- sort(runif(400, -2, 2))

  # BELOW the null space is rejected where the two numbers were written
  expect_error(bspline_smooth(k = 20, degree = 5, order = 3, constrain = 1),
               "does not contain the null space")
  # at or above it is accepted, and each degree costs one dimension
  d <- vapply(list(NULL, 2L, 3L), function(cc) {
    ncol(smoother_build(bspline_smooth(k = 20, constrain = cc,
                                       null_space = "drop"), x)$X)
  }, integer(1))
  expect_identical(d, c(18L, 17L, 16L))

  # and constraining against x^2 makes the block orthogonal to it
  o0 <- smoother_build(bspline_smooth(k = 20, null_space = "drop"), x)
  o2 <- smoother_build(bspline_smooth(k = 20, constrain = 2,
                                      null_space = "drop"), x)
  expect_gt(max(abs(cor(x^2, o0$X))), 0.9)
  expect_lt(max(abs(cor(x^2, o2$X))), 1e-9)
})

test_that("the measure changes the penalty and so the block", {
  set.seed(5)
  x <- sort(runif(300, -2, 2))
  a <- bspline_smooth(k = 12, lower = -2, upper = 2)
  b <- bspline_smooth(k = 12, lower = -2, upper = 2, measure = "empirical")
  cf <- bspline_smooth(k = 12, lower = -2, upper = 2,
                       measure = function(u) dnorm(u, 0, 0.25))
  xa <- smoother_build(a, x)$X
  expect_gt(max(abs(xa - smoother_build(b, x)$X)), 1e-6)
  expect_gt(max(abs(xa - smoother_build(cf, x)$X)), 1e-6)

  # a measure that is not one of the two names and not a function is refused
  expect_error(bspline_smooth(k = 12, measure = "gaussian"), "must be")
  expect_error(bspline_smooth(k = 12, measure = 3), "must be")
})

test_that("a Fourier smoother is periodic, which is what it is for", {
  set.seed(3)
  x <- sort(runif(300))
  f <- function(u) sin(2 * pi * u) + 0.4 * cos(4 * pi * u)
  y <- f(x) + rnorm(300, sd = 0.2)

  # THE DEFAULT PENALTY IS THE HARMONIC OPERATOR since operators
  # arrived, so the constant AND the fundamental leave the penalized
  # part; 'shrink' then gives the fundamental a weight of its own rather
  # than leaving it free, and no column is unpenalized.
  sm <- fourier_smooth(k = 9, lower = 0, upper = 1)
  o <- smoother_build(sm, x)

  expect_identical(ncol(o$X), 8L)
  expect_identical(o$unpenalized, 0L)
  expect_identical(head(o$names, 2L), c("sin1", "cos1"))

  # and the derivative penalty is the construction this family had
  # before: one constant removed, nothing restored
  od <- smoother_build(fourier_smooth(k = 9, order = 2, lower = 0,
                                      upper = 1), x)
  expect_identical(ncol(od$X), 8L)
  expect_identical(od$unpenalized, 0L)
  expect_identical(od$names, paste0("z", 1:8))

  # EVERY COLUMN is periodic, hence so is any fit built on them
  ends <- smoother_apply(sm, o$blueprint, c(0, 1))
  expect_lt(max(abs(ends[1L, ] - ends[2L, ])), 1e-12)
  b <- pfit(o$X, o$S, y, 0.01)
  expect_lt(abs(diff(as.vector(ends %*% b))), 1e-12)

  # and it recovers a periodic truth
  expect_lt(sqrt(mean((o$X %*% b - f(x))^2)), 0.05)
})

test_that("the Fourier null space is the constant at every order", {
  set.seed(3)
  x <- sort(runif(200))
  for (m in 1:3) {
    sp <- smoother_span(fourier_smooth(k = 9, order = m, lower = 0, upper = 1), x)
    expect_identical(ncol(sp$constraint), 1L)
    expect_identical(ncol(sp$free), 0L)
    o <- smoother_build(fourier_smooth(k = 9, order = m, lower = 0, upper = 1), x)
    # the dimension does NOT fall with the order, as it does for a B-spline
    expect_identical(ncol(o$X), 8L)
  }
})

test_that("a Fourier smoother declines the arguments it has no reading for", {
  # not by refusing them: they are not on the signature, so R reports them
  expect_error(fourier_smooth(k = 9, degree = 3), "unused argument")
  expect_error(fourier_smooth(k = 9, constrain = 2), "unused argument")
  # and an even k would leave half a sine-cosine pair
  expect_error(fourier_smooth(k = 8), "must be odd")
})

test_that("a Legendre smoother has a polynomial null space like a B-spline", {
  set.seed(3)
  x <- sort(runif(200, -1, 1))
  o2 <- smoother_build(legendre_smooth(k = 8), x)
  expect_identical(ncol(o2$X), 7L)
  expect_identical(o2$unpenalized, 1L)
  expect_identical(o2$names[[1L]], "lin")

  o3 <- smoother_build(legendre_smooth(k = 8, order = 3), x)
  expect_identical(o3$unpenalized, 2L)

  # order above the highest degree would give the zero matrix
  expect_error(legendre_smooth(k = 3, order = 4), "exceeds the highest degree")

  # the two families fit the same smooth truth comparably
  y <- sin(3 * x) + rnorm(200, sd = 0.2)
  rb <- sqrt(mean((smoother_build(bspline_smooth(k = 12), x)$X %*%
    pfit(smoother_build(bspline_smooth(k = 12), x)$X,
         smoother_build(bspline_smooth(k = 12), x)$S, y, 0.1) - sin(3 * x))^2))
  rl <- sqrt(mean((o2$X %*% pfit(o2$X, o2$S, y, 0.1) - sin(3 * x))^2))
  expect_lt(abs(rb - rl), 0.1)
})

test_that("a smoother of another family prints its own decisions", {
  out <- paste(capture.output(print(fourier_smooth(k = 9, lower = 0,
                                                   upper = 12))),
               collapse = " ")
  expect_match(out, "FourierSmoother")
  expect_match(out, "9 functions")
  expect_match(out, "\\[0, 12\\]")
  # a Fourier has no degree to print
  expect_false(grepl("degree", out))
})
