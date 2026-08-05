# Contracting against coefficients: the same answer as the design matrix, and
# without building it.

tb <- function(dims = c(5L, 4L)) {
  # The degree is capped at what the dimension allows, so a small marginal is
  # a lower-degree spline rather than an error.
  tensor_basis(lapply(dims, function(k) {
    bspline_basis(dimension = k, degree = min(3L, k - 1L))
  }))
}


test_that("an ordinary basis contracts to the design matrix times the coefficients", {
  set.seed(51)
  b <- bspline_basis(dimension = 6)
  x <- runif(10)
  cf <- rnorm(6)
  expect_equal(basis_contract(b, x, cf), drop(basis_eval(b, x) %*% cf))

  # several sets of coefficients at once
  cm <- matrix(rnorm(18), 6, 3)
  expect_equal(basis_contract(b, x, cm), basis_eval(b, x) %*% cm)
  expect_error(basis_contract(b, x, rnorm(5)), "6 entries")
})


test_that("a full array contracts to exactly what the design matrix gives", {
  # The array's dimensions must line up with the column order, which runs with
  # the LAST marginal fastest while an R array runs with the first fastest.
  # Flattening without accounting for that pairs every coefficient with the
  # wrong function and returns a perfectly finite wrong number.
  set.seed(52)
  b <- tb(c(5L, 4L))
  x <- cbind(runif(12), runif(12))
  cf <- array(rnorm(20), dim = c(5, 4))

  want <- drop(basis_eval(b, x) %*% as.numeric(aperm(cf)))
  expect_equal(basis_contract(b, x, cf), want)

  # the coefficient of column "bs1.bs2" is cf[1, 2], read off directly
  j <- match("bs1.bs2", basis_colnames(b))
  probe <- array(0, dim = c(5, 4))
  probe[1, 2] <- 1
  expect_equal(basis_contract(b, x, probe), basis_eval(b, x)[, j],
    ignore_attr = TRUE
  )
})


test_that("the block size changes nothing but the memory", {
  set.seed(53)
  b <- tb(c(4L, 3L))
  x <- cbind(runif(37), runif(37))
  cf <- array(rnorm(12), dim = c(4, 3))
  full <- basis_contract(b, x, cf)
  for (blk in c(1L, 2L, 7L, 1000L)) {
    expect_equal(basis_contract(b, x, cf, block = blk), full,
      label = paste("block", blk)
    )
  }
})


test_that("the canonical polyadic form gives the array it encodes", {
  set.seed(54)
  b <- tb(c(5L, 4L))
  x <- cbind(runif(15), runif(15))

  for (rank in 1:3) {
    g <- list(matrix(rnorm(5 * rank), 5, rank), matrix(rnorm(4 * rank), 4, rank))
    arr <- Reduce(`+`, lapply(seq_len(rank), function(f) {
      outer(g[[1L]][, f], g[[2L]][, f])
    }))
    expect_equal(basis_contract(b, x, g), basis_contract(b, x, arr),
      label = paste("rank", rank)
    )
  }
})


test_that("a full-rank factorization reproduces any array, in three variables", {
  # The exit test of the form: with one rank-one term per coefficient the
  # factorization is exact, so nothing is lost by the representation itself.
  set.seed(55)
  dims <- c(3L, 3L, 3L)
  b <- tensor_basis(lapply(dims, function(k) {
    bspline_basis(dimension = k, degree = 2)
  }))
  x <- cbind(runif(8), runif(8), runif(8))
  arr <- array(rnorm(27), dim = dims)

  idx <- expand.grid(seq_len(3), seq_len(3), seq_len(3))
  g <- lapply(seq_len(3), function(j) {
    m <- matrix(0, 3, 27)
    m[cbind(idx[[j]], seq_len(27))] <- 1
    m
  })
  g[[1L]] <- g[[1L]] * rep(as.numeric(arr), each = 3)

  expect_equal(basis_contract(b, x, g), basis_contract(b, x, arr))
  expect_equal(basis_contract(b, x, arr),
    drop(basis_eval(b, x) %*% as.numeric(aperm(arr)))
  )
})


test_that("the factorized form never builds the design matrix", {
  # The claim the form exists for. A product of six bases of ten functions has
  # a million columns, so a design matrix for a thousand points would be eight
  # gigabytes; the factorization touches neither it nor the coefficient array.
  set.seed(56)
  d <- 6L
  k <- 10L
  b <- tensor_basis(lapply(seq_len(d), function(i) {
    bspline_basis(dimension = k)
  }))
  expect_identical(b@dimension, 1000000L)

  n <- 1000L
  x <- matrix(runif(n * d), n, d)
  rank <- 4L
  g <- lapply(seq_len(d), function(i) matrix(rnorm(k * rank), k, rank))

  before <- gc(reset = TRUE)
  value <- basis_contract(b, x, g)
  peak_mb <- sum(gc()[, "max used"] * c(56, 8)) / 1024^2

  expect_length(value, n)
  expect_true(all(is.finite(value)))
  # the design matrix alone would be n * 10^6 * 8 bytes, about 7600 MB
  expect_lt(peak_mb, 500)

  # and the answer is right: checked against the definition, one row at a time
  phi <- lapply(seq_len(d), function(j) {
    basis_eval(b@marginals[[j]], x[, j]) %*% g[[j]]
  })
  expect_equal(value, rowSums(Reduce(`*`, phi)))
})


test_that("malformed coefficients are refused by name", {
  b <- tb(c(5L, 4L))
  x <- cbind(0.3, 0.7)

  expect_error(basis_contract(b, x, rnorm(19)), "20 values")
  expect_error(basis_contract(b, x, array(0, dim = c(4, 5))), "5 by 4")
  expect_error(
    basis_contract(b, x, list(matrix(0, 5, 2))),
    "one factor matrix per marginal"
  )
  expect_error(
    basis_contract(b, x, list(matrix(0, 5, 2), matrix(0, 4, 3))),
    "same number of columns"
  )
  expect_error(
    basis_contract(b, x, list(matrix(0, 5, 2), matrix(0, 3, 2))),
    "3 rows, but its marginal has 4"
  )
})


test_that("a flat vector of coefficients works, in the column order", {
  set.seed(57)
  b <- tb(c(3L, 2L))
  x <- cbind(runif(5), runif(5))
  cf <- rnorm(6)
  expect_equal(basis_contract(b, x, cf), drop(basis_eval(b, x) %*% cf))
})
