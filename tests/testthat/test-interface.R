# The contract every basis keeps, whatever family it belongs to.

bases <- function() {
  list(
    bspline = bspline_basis(dimension = 6, degree = 3),
    bspline_shifted = bspline_basis(lower = -2, upper = 5, dimension = 8),
    bspline_deg1 = bspline_basis(dimension = 4, degree = 1),
    fourier = fourier_basis(dimension = 5),
    fourier_shifted = fourier_basis(lower = -1, upper = 3, dimension = 7)
  )
}


test_that("every generic returns a matrix of the declared shape and names", {
  for (nm in names(bases())) {
    b <- bases()[[nm]]
    x <- seq(b@lower, b@upper, length.out = 9)
    k <- b@dimension

    for (m in list(
      basis_eval(b, x),
      basis_deriv(b, x, order = 1L),
      basis_deriv(b, x, order = 2L),
      basis_int(b, x)
    )) {
      expect_true(is.matrix(m), label = nm)
      expect_identical(dim(m), c(9L, k), label = nm)
      expect_identical(colnames(m), basis_colnames(b), label = nm)
    }

    # a single point still gives a one-row matrix, not a dropped vector
    expect_identical(dim(basis_eval(b, x[1])), c(1L, k), label = nm)

    g <- basis_gram(b)
    expect_identical(dim(g), c(k, k), label = nm)
    expect_equal(g, t(g), label = nm)
  }
})


test_that("order 0 of basis_deriv is the basis itself", {
  b <- bspline_basis(dimension = 6)
  x <- seq(0, 1, length.out = 7)
  expect_identical(basis_deriv(b, x, order = 0L), basis_eval(b, x))
})


test_that("the integral is zero at the lower endpoint, exactly", {
  # Any antiderivative satisfies the differentiation check, so the constant has
  # to be pinned by the contract; a basis that returned a different one would
  # be silently unusable in a sum with one that did not.
  for (nm in names(bases())) {
    b <- bases()[[nm]]
    expect_true(all(basis_int(b, b@lower) == 0), label = nm)
  }
})


test_that("evaluation outside the interval is refused", {
  b <- bspline_basis(lower = 0, upper = 1, dimension = 5)
  expect_error(basis_eval(b, 1.5), "outside the basis interval")
  expect_error(basis_eval(b, -0.1), "outside the basis interval")
  expect_error(basis_deriv(b, 2), "outside the basis interval")
  expect_error(basis_int(b, 2), "outside the basis interval")
  expect_error(basis_eval(b, "a"), "must be numeric")

  # a point that is an endpoint up to rounding is accepted and clamped
  expect_silent(basis_eval(b, 1 + 1e-12))
})


test_that("missing values give missing rows and nothing else", {
  for (nm in names(bases())) {
    b <- bases()[[nm]]
    x <- c(b@lower + 0.1, NA_real_, b@upper - 0.1)
    for (m in list(
      basis_eval(b, x), basis_deriv(b, x, order = 1L), basis_int(b, x)
    )) {
      expect_true(all(is.na(m[2L, ])), label = nm)
      expect_false(anyNA(m[c(1L, 3L), ]), label = nm)
    }
  }
})


test_that("an invalid order is refused", {
  b <- bspline_basis(dimension = 5)
  expect_error(basis_deriv(b, 0.5, order = -1), "non-negative integer")
  expect_error(basis_deriv(b, 0.5, order = 1.5), "non-negative integer")
  expect_error(basis_deriv(b, 0.5, order = c(1, 2)), "non-negative integer")
  expect_error(basis_gram(b, order = -1), "non-negative integer")
})


test_that("the constructors refuse an impossible interval or dimension", {
  expect_error(bspline_basis(lower = 1, upper = 0), "strictly less")
  expect_error(bspline_basis(dimension = 0), "positive integer")
  expect_error(bspline_basis(dimension = 2.5), "positive integer")
  expect_error(bspline_basis(lower = NA), "finite")
  expect_error(fourier_basis(lower = 0, upper = 0), "strictly less")
})


test_that("print reports the family, the interval and what is numerical", {
  out <- paste(capture.output(print(bspline_basis(dimension = 6))),
    collapse = "\n"
  )
  expect_match(out, "bspline")
  expect_match(out, "Functions: 6")
  expect_match(out, "Numerical: none")

  Bare <- S7::new_class("Bare", parent = basis)
  S7::method(basis_eval, Bare) <- function(basis, x, ...) {
    matrix(x, length(x), 1L, dimnames = list(NULL, basis_colnames(basis)))
  }
  bare <- Bare(
    basis_name = "bare", dimension = 1L, lower = 0, upper = 1,
    basis_params = list()
  )
  out2 <- paste(capture.output(print(bare)), collapse = "\n")
  expect_match(out2, "basis_deriv")
})


test_that("plot draws without error for every order it accepts", {
  b <- bspline_basis(dimension = 6)
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({
    grDevices::dev.off()
    unlink(f)
  })
  expect_silent(plot(b))
  expect_silent(plot(b, order = 1))
  expect_silent(plot(b, order = -1))
  expect_error(plot(b, order = -2), "must be -1")
})
