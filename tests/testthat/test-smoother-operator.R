# A smoother penalized by a linear differential operator, and what the new
# default of fourier_smooth() changed.
#
# The claim that nothing ELSE moved is pinned in test-smoother-axes.R and by
# the identity below: `order = m` and `order = deriv_operator(m)` build the
# same object bit for bit, so every construction written before operators
# existed takes the route it always took.

set.seed(9)
n <- 300L
x <- sort(runif(n))
f <- function(u) sin(2 * pi * u) + 0.4 * cos(4 * pi * u)
y <- f(x) + rnorm(n, sd = 0.2)

test_that("a whole order and the operator that spells it build the same thing", {
  for (sm in list(
    list(bspline_smooth(k = 10, order = 2), bspline_smooth(k = 10, order = deriv_operator(2))),
    list(bspline_smooth(k = 12, order = 3), bspline_smooth(k = 12, order = deriv_operator(3))),
    list(legendre_smooth(k = 8, order = 2), legendre_smooth(k = 8, order = deriv_operator(2))),
    list(fourier_smooth(k = 9, order = 2, lower = 0, upper = 1),
         fourier_smooth(k = 9, order = deriv_operator(2), lower = 0, upper = 1)),
    list(cyclic_smooth(k = 10, order = 2, lower = 0, upper = 1),
         cyclic_smooth(k = 10, order = deriv_operator(2), lower = 0, upper = 1))
  )) {
    a <- smoother_build(sm[[1L]], x)
    b <- smoother_build(sm[[2L]], x)
    expect_identical(a$X, b$X)
    expect_identical(a$S, b$S)
    expect_identical(a$unpenalized, b$unpenalized)
    expect_identical(a$names, b$names)
  }
})

test_that("the escape hatch is the derivative Gram, at the numbers", {
  # ⚠️ A PROPERTY IS NOT AN IDENTITY. The block's column count, its
  # unpenalized count and its coefficient names can all be right while the
  # numbers have shifted, and `order = 2` promises the numbers: it is the
  # construction the family had before operators existed. What can be
  # asserted here is the numerical claim that promise rests on -- the
  # smoother's roughness matrix under a whole order IS the plain
  # derivative Gram of its basis, with the operator machinery contributing
  # nothing to it.
  #
  # The identity against the version that shipped cannot live here, that
  # version not being installed at test time. It was taken once, at the
  # release: nine smoother shapes and a fitted model, block, penalty,
  # reapplied block, roughness, logLik, edf, coef, fitted and vcov, all
  # identical between 0.12.0 and 0.13.0.
  for (m in 1:3) {
    for (ms in list("lebesgue", "empirical")) {
      sm <- fourier_smooth(k = 11, order = m, lower = 0, upper = 1,
                           measure = ms)
      b <- smoother_basis(sm, x)
      g <- smoother_gram(sm, b, x)
      want <- if (identical(ms, "empirical")) {
        basis_gram(b, order = m, at = x)
      } else {
        basis_gram(b, order = m)
      }
      expect_identical(g, want)
    }
  }
  # and the same for the two other families that carry an order
  for (sm in list(bspline_smooth(k = 12, order = 2),
                  legendre_smooth(k = 9, order = 2),
                  cyclic_smooth(k = 12, order = 2, lower = 0, upper = 1))) {
    b <- smoother_basis(sm, x)
    expect_identical(smoother_gram(sm, b, x), basis_gram(b, order = 2L))
  }
})

test_that("fourier_smooth defaults to the harmonic operator and shrinks", {
  sm <- fourier_smooth(k = 9, lower = 0, upper = 365)
  expect_false(basis7:::is_deriv_operator(sm@order))
  expect_identical(sm@null_space, "shrink")
  # the period is the interval's, resolved at build
  expect_equal(operator_weights(basis7:::smoother_operator(sm, c(0, 365))),
               operator_weights(harmonic_operator(365)))
})

test_that("the null space is the constant and the fundamental", {
  sk <- fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep")
  o <- smoother_build(sk, x)
  # the constant goes to the model's intercept and is not restored; the
  # fundamental pair is, and the penalty leaves exactly those two alone
  expect_identical(o$unpenalized, 2L)
  expect_identical(head(o$names, 2L), c("sin1", "cos1"))
  expect_identical(ncol(o$X), 8L)

  # with two harmonics it is four
  s2 <- fourier_smooth(k = 13, lower = 0, upper = 1, null_space = "keep",
                       order = harmonic_operator(harmonics = 2))
  o2 <- smoother_build(s2, x)
  expect_identical(o2$unpenalized, 4L)
  expect_identical(head(o2$names, 4L), c("sin1", "cos1", "sin2", "cos2"))

  # and a derivative penalty restores nothing, as it always did
  od <- smoother_build(fourier_smooth(k = 9, order = 2, lower = 0, upper = 1), x)
  expect_identical(od$unpenalized, 0L)
  expect_identical(od$names[[1L]], "z1")
})

test_that("a strongly penalized fit contracts to the fundamental", {
  # THE PROPERTY THE OPERATOR EXISTS FOR. Under the harmonic penalty with the
  # null space kept, a large smoothing parameter leaves a pure sinusoid;
  # under the derivative penalty it leaves a constant.
  sinu <- function(v) sqrt(mean(resid(lm(v ~ sin(2 * pi * x) + cos(2 * pi * x)))^2))
  fit <- function(sm, lam) {
    o <- smoother_build(sm, x)
    b <- solve(crossprod(o$X) + lam * o$S, crossprod(o$X, y - mean(y)))
    as.vector(o$X %*% b)
  }
  fh <- fit(fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep"), 1e8)
  fd <- fit(fourier_smooth(k = 9, order = 2, lower = 0, upper = 1), 1e8)

  expect_lt(sinu(fh) / stats::sd(fh), 1e-5)   # a sinusoid, and not a flat line
  expect_gt(stats::sd(fh), 0.3)
  expect_lt(stats::sd(fd), 1e-3)              # the derivative penalty flattens
})

test_that("every column stays periodic under either penalty", {
  # the one property a Fourier basis is chosen for, and the reason a
  # non-periodic null function is refused
  for (sm in list(
    fourier_smooth(k = 9, lower = 0, upper = 1),
    fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep"),
    fourier_smooth(k = 13, lower = 0, upper = 1, null_space = "keep",
                   order = harmonic_operator(harmonics = 2)),
    fourier_smooth(k = 9, order = 2, lower = 0, upper = 1)
  )) {
    o <- smoother_build(sm, x)
    ends <- smoother_apply(sm, o$blueprint, c(0, 1))
    expect_lt(max(abs(ends[1L, ] - ends[2L, ])), 1e-12)
  }
})

test_that("the block is reapplied at new values and not rebuilt", {
  for (sm in list(
    fourier_smooth(k = 9, lower = 0, upper = 1, null_space = "keep"),
    fourier_smooth(k = 13, lower = 0, upper = 1, null_space = "keep",
                   order = harmonic_operator(harmonics = 2)),
    bspline_smooth(k = 20, order = harmonic_operator(1))
  )) {
    o <- smoother_build(sm, x)
    i <- c(3L, 44L, 180L, 299L)
    expect_equal(smoother_apply(sm, o$blueprint, x[i]), o$X[i, , drop = FALSE],
                 tolerance = 1e-12, ignore_attr = TRUE)
  }
})

test_that("a periodic basis refuses a null function it cannot carry", {
  # t is not periodic, so restoring it would cost the fit the property the
  # basis was chosen for
  expect_error(
    smoother_build(fourier_smooth(k = 9, lower = 0, upper = 1,
                                  order = deriv_operator(2) * oscillator_operator(1)), x),
    "cannot carry"
  )
  # and neither is a frequency that is not a multiple of the basis's
  expect_error(
    smoother_build(fourier_smooth(k = 9, lower = 0, upper = 1,
                                  order = harmonic_operator(0.37)), x),
    "cannot carry"
  )
  # while a multiple of it is carried
  expect_silent(smoother_build(fourier_smooth(k = 13, lower = 0, upper = 1,
                                              order = harmonic_operator(0.5)), x))
})

test_that("every family that takes constrain still builds with it", {
  # ⚠️ REGRESSION GUARD. check_constrain() takes an operator since
  # operators arrived, and pspline_smooth() and adaptive_smooth() went on
  # handing it the raw integer `diff`, so is_deriv_operator() read
  # @operator_params off an integer and both families raised the moment a
  # `constrain` was given. Nothing covered it: `constrain` was exercised
  # on bspline alone, and the two difference-penalty families reach the
  # same helper by another path.
  #
  # The four blocks were measured identical to what 0.12.0 built, so what
  # this asserts is that they build at all and that the constraint removes
  # what it says it removes.
  for (sm in list(
    bspline_smooth(k = 20, constrain = 2L, lower = 0, upper = 1),
    legendre_smooth(k = 12, constrain = 2L),
    pspline_smooth(k = 20, constrain = 2L, lower = 0, upper = 1),
    adaptive_smooth(k = 20, m = 4, constrain = 2L, lower = 0, upper = 1)
  )) {
    o <- smoother_build(sm, x)
    expect_true(is.matrix(o$X))
    expect_identical(o$unpenalized, 1L)
    expect_identical(o$names[[1L]], "lin")
  }
})

test_that("constrain and a non-polynomial null space are refused together", {
  expect_error(bspline_smooth(k = 12, order = harmonic_operator(1), constrain = 2),
               "not the polynomials")
  # and the polynomial case is untouched
  expect_silent(bspline_smooth(k = 12, order = 2, constrain = 2))
})

test_that("a spline carries a periodic operator, with the null space declared", {
  # ADMISSIBLE THOUGH THE BASIS DOES NOT CONTAIN THE NULL SPACE EXACTLY: the
  # free columns come from the operator, so the pencil is well posed even
  # where the roughness matrix is only nearly singular.
  sm <- bspline_smooth(k = 20, order = harmonic_operator(1))
  o <- smoother_build(sm, x)
  expect_identical(o$unpenalized, 2L)
  expect_identical(head(o$names, 2L), c("sin1", "cos1"))
  b <- solve(crossprod(o$X) + 1e7 * o$S, crossprod(o$X, y - mean(y)))
  fv <- as.vector(o$X %*% b)
  expect_lt(sqrt(mean(resid(lm(fv ~ sin(2 * pi * x) + cos(2 * pi * x)))^2)) /
              stats::sd(fv), 1e-4)
})

test_that("a penalty factory still works on the family whose default shrinks", {
  # ⚠️ check_available() refuses a factory together with null_space =
  # "shrink", so a family defaulting to "shrink" would break a spelling that
  # has always worked. The default resolves to "keep" where a factory is
  # given, and the two are still refused when both are named.
  f <- function(n_coef) n_coef
  sm <- fourier_smooth(k = 9, lower = 0, upper = 1, penalty = f)
  expect_identical(sm@null_space, "keep")
  expect_error(fourier_smooth(k = 9, lower = 0, upper = 1, penalty = f,
                              null_space = "shrink"), "cannot both be given")
})

test_that("the smoother prints which operator it penalizes with", {
  out <- capture.output(print(fourier_smooth(k = 9, lower = 0, upper = 365)))
  expect_true(any(grepl("harmonic acceleration operator of order 3", out)))
  out2 <- capture.output(print(fourier_smooth(k = 9, order = 2, lower = 0, upper = 1)))
  expect_true(any(grepl("derivative of order 2", out2)))
})
