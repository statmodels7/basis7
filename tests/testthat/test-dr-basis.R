# The Demmler-Reinsch construction: three properties, and the design it
# survives that the alternatives do not.

test_that("the three defining properties hold exactly", {
  set.seed(31)
  x <- sort(stats::runif(300))
  b <- bspline_basis(dimension = 12)
  d <- dr_basis(b, x)

  z <- basis_eval(d, x)
  zz <- crossprod(z)
  p <- basis_gram(b, order = 2L)
  tpt <- crossprod(d@transform, p %*% d@transform)

  # 1. the empirical inner product is diagonal
  expect_lt(max(abs(zz[upper.tri(zz)])), 1e-9)
  # 2. the penalty is a scalar multiple of the identity
  expect_lt(max(abs(tpt[upper.tri(tpt)])), 1e-9)
  expect_lt(diff(range(diag(tpt))) / mean(diag(tpt)), 1e-9)
  # 3. empirically orthogonal to a constant and to x
  expect_lt(max(abs(crossprod(cbind(1, x), z))), 1e-9)

  # the dimension falls by the rank of the constraint
  expect_identical(d@dimension, b@dimension - 2L)
  # and the scaling puts the trace at one
  expect_equal(sum(diag(zz)) / length(x), 1, tolerance = 1e-10)
})


test_that("the properties survive a RANK-DEFICIENT design", {
  # Equally spaced knots with a gap in the data: the Schoenberg-Whitney
  # interlacing fails and B loses rank. Any construction that factorizes the
  # design cannot start here, which is the reason this one factorizes the
  # penalty instead. Checked rather than believed.
  # The gap is placed deterministically rather than drawn, so the test asserts
  # a property of the construction and not of one random seed.
  b <- bspline_basis(dimension = 14, degree = 3)
  x <- sort(c(
    seq(0, 0.25, length.out = 120),
    seq(0.78, 1, length.out = 120)
  ))

  bm <- basis_eval(b, x)
  expect_lt(qr(bm)$rank, ncol(bm)) # genuinely rank-deficient
  expect_s3_class(try(chol(crossprod(bm)), silent = TRUE), "try-error")

  d <- dr_basis(b, x)
  z <- basis_eval(d, x)
  zz <- crossprod(z)
  tpt <- crossprod(d@transform, basis_gram(b, order = 2L) %*% d@transform)

  expect_true(all(is.finite(z)))
  expect_lt(max(abs(zz[upper.tri(zz)])), 1e-8)
  expect_lt(max(abs(tpt[upper.tri(tpt)])), 1e-8)
  expect_lt(diff(range(diag(tpt))) / mean(diag(tpt)), 1e-8)
  expect_lt(max(abs(crossprod(cbind(1, x), z))), 1e-8)
})


test_that("prediction at new points goes through the kept transform", {
  # Nothing of the size of the sample is stored, so a new point costs the
  # parent's evaluation and one matrix multiplication.
  set.seed(33)
  x <- sort(stats::runif(200))
  b <- bspline_basis(dimension = 10)
  d <- dr_basis(b, x)

  xnew <- seq(0, 1, length.out = 17)
  expect_equal(basis_eval(d, xnew), basis_eval(b, xnew) %*% d@transform,
    ignore_attr = TRUE
  )
  expect_true(all(is.finite(basis_eval(d, xnew))))
})


test_that("a supplied penalty is used, and a discrete one works", {
  set.seed(34)
  x <- sort(stats::runif(250))
  b <- bspline_basis(dimension = 12)
  dd <- crossprod(diff(diag(12), differences = 2L))

  d <- dr_basis(b, x, penalty = dd)
  tpt <- crossprod(d@transform, dd %*% d@transform)
  expect_lt(max(abs(tpt[upper.tri(tpt)])), 1e-9)
  expect_lt(diff(range(diag(tpt))) / mean(diag(tpt)), 1e-9)

  # a different penalty gives a different basis
  expect_gt(max(abs(d@transform - dr_basis(b, x)@transform)), 1e-6)
})


test_that("the constraints are configurable, and default to 1 and x", {
  set.seed(35)
  x <- sort(stats::runif(200))
  b <- bspline_basis(dimension = 10)

  # A constraint has to remove the penalty's null space, or the directions left
  # over are neither penalized nor identified. Removing the constant alone goes
  # with the first-derivative penalty, whose null space is the constants.
  d1 <- dr_basis(b, x,
    penalty = basis_gram(b, order = 1L),
    constraints = matrix(1, 1L, length(x))
  )
  expect_identical(d1@dimension, 9L)
  expect_lt(max(abs(colSums(basis_eval(d1, x)))), 1e-9)

  # and the default removes two, which is what the second-derivative penalty
  # needs: its null space is the straight lines.
  expect_identical(dr_basis(b, x)@dimension, 8L)

  # The mismatch is refused rather than regularized, and deterministically:
  # the linear direction is exactly unpenalized here, so whether chol() happens
  # to survive a zero pivot must not be what decides.
  expect_error(
    dr_basis(b, x, constraints = matrix(1, 1L, length(x))),
    "neither penalized nor identified"
  )
})


test_that("scaling can be turned off and changes only the scale", {
  set.seed(36)
  x <- sort(stats::runif(200))
  b <- bspline_basis(dimension = 9)
  a <- dr_basis(b, x, scale = TRUE)
  n <- dr_basis(b, x, scale = FALSE)

  ratio <- n@transform / a@transform
  expect_lt(diff(range(ratio)) / abs(mean(ratio)), 1e-9)
  expect_equal(sum(diag(crossprod(basis_eval(a, x)))) / length(x), 1,
    tolerance = 1e-10
  )
})


test_that("a singular penalty on the constrained space is refused by name", {
  # A direction that is neither penalized nor identified is not a basis
  # function, and saying so beats regularizing quietly.
  expect_error(
    dr_basis(bspline_basis(dimension = 5), seq(0, 1, length.out = 40),
      penalty = matrix(0, 5, 5)
    ),
    "neither penalized nor identified"
  )
})


test_that("malformed arguments are refused", {
  b <- bspline_basis(dimension = 6)
  x <- seq(0, 1, length.out = 30)
  expect_error(dr_basis(b, x, penalty = diag(5)), "6 by 6")
  expect_error(dr_basis(b, c(x, NA)), "missing values")
  expect_error(dr_basis(b, x, constraints = matrix(1, 1L, 5L)), "one column per")
  expect_error(dr_basis(b, c(0, 2)), "outside the basis interval")
})


test_that("check_basis passes on a DR basis", {
  set.seed(37)
  d <- dr_basis(bspline_basis(dimension = 10), sort(stats::runif(200)))
  r <- check_basis(d, verbose = FALSE)
  expect_true(all(r, na.rm = TRUE))
})
