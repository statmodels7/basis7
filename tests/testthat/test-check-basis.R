# The validator, and the two failure modes it has to tell apart.

broken <- function(generic, f, parent = BsplineBasis) {
  cls <- S7::new_class(paste0("Broken", sample.int(1e6, 1L)), parent = parent)
  S7::method(generic, cls) <- f
  b <- bspline_basis(dimension = 6)
  cls(
    basis_name = "broken", dimension = b@dimension, lower = b@lower,
    upper = b@upper, basis_params = b@basis_params
  )
}


test_that("the shipped bases pass every check that applies to them", {
  r <- check_basis(bspline_basis(dimension = 6), verbose = FALSE)
  expect_true(all(r, na.rm = TRUE))
  expect_false(anyNA(r))

  r2 <- check_basis(fourier_basis(dimension = 5), verbose = FALSE)
  expect_true(all(r2, na.rm = TRUE))
  expect_true(is.na(r2[["partition"]])) # the family never claimed it
})


test_that("a WRONG derivative is caught", {
  b <- broken(basis_deriv, function(basis, x, order = 1L, ...) {
    1.05 * basis7:::bspline_design(basis, x, derivs = order)
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["deriv"]])
})


test_that("a MISSING derivative is reported as unchecked, never as passed", {
  # Absent and wrong are different failure modes, and a validator needs a probe
  # for each. Reporting a numerical value as verified would be worse than
  # reporting nothing: it claims a check was made that was not.
  Bare <- S7::new_class("BareBasis", parent = basis)
  S7::method(basis_eval, Bare) <- function(basis, x, ...) {
    out <- cbind(1, x, x^2)
    colnames(out) <- basis_colnames(basis)
    out
  }
  bare <- Bare(
    basis_name = "bare", dimension = 3L, lower = 0, upper = 1,
    basis_params = list()
  )
  r <- check_basis(bare, verbose = FALSE)
  expect_true(is.na(r[["deriv"]]))
  expect_true(is.na(r[["integral"]]))
  expect_false(isTRUE(r[["deriv"]]))
  expect_true(all(attr(r, "numerical")))
})


test_that("an integral with the wrong constant is caught", {
  # It still differentiates back to the basis, so only the anchor check finds
  # it: this is exactly why the convention is part of the contract.
  b <- broken(basis_int, function(basis, x, ...) {
    basis7:::bspline_design(basis, x, integral = TRUE) + 0.1
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["integral"]])
})


test_that("an integral of the wrong function is caught", {
  b <- broken(basis_int, function(basis, x, ...) {
    1.05 * basis7:::bspline_design(basis, x, integral = TRUE)
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["integral"]])
})


test_that("a wrong Gram matrix is caught", {
  b <- broken(basis_gram, function(basis, order = 0L, ...) {
    1.05 * basis7:::numerical_gram(basis, order, panels = 60L, nodes = 10L)
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["gram"]])
})


test_that("a non-symmetric Gram matrix is caught", {
  b <- broken(basis_gram, function(basis, order = 0L, ...) {
    g <- basis7:::numerical_gram(basis, order, panels = 60L, nodes = 10L)
    g[1L, 2L] <- g[1L, 2L] + 1
    g
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["gram"]])
})


test_that("a basis that loses its column names is caught", {
  b <- broken(basis_eval, function(basis, x, ...) {
    unname(basis7:::bspline_design(basis, x))
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["shape"]])
})


test_that("a basis that drops missing values is caught", {
  b <- broken(basis_eval, function(basis, x, ...) {
    m <- basis7:::bspline_design(basis, x)
    m[is.na(x), ] <- 0
    colnames(m) <- basis_colnames(basis)
    m
  })
  r <- check_basis(b, verbose = FALSE)
  expect_false(r[["missing"]])
})


test_that("the printed report names each outcome", {
  out <- paste(capture.output(check_basis(bspline_basis(dimension = 5))),
    collapse = "\n"
  )
  expect_match(out, "check_basis: bspline")
  expect_match(out, "\\[PASSED\\]")

  out2 <- paste(capture.output(check_basis(fourier_basis(dimension = 5))),
    collapse = "\n"
  )
  # a property the family never claimed is not the same as one left unchecked
  expect_match(out2, "\\[not claimed\\]")

  b <- broken(basis_deriv, function(basis, x, order = 1L, ...) {
    1.05 * basis7:::bspline_design(basis, x, derivs = order)
  })
  out3 <- paste(capture.output(check_basis(b)), collapse = "\n")
  expect_match(out3, "\\[FAILED\\]")
})


test_that("the comparison is relative to the values themselves", {
  # A B-spline is zero on most of its interval, so a denominator floored at one
  # would flatten a real disagreement between small numbers into agreement.
  # Two values half of one another disagree whatever their size:
  expect_false(basis7:::rel_close(matrix(1e-9), matrix(1.5e-9), 1e-6))
  expect_false(basis7:::rel_close(
    cbind(c(1, 1e-4)), cbind(c(1, 1.5e-4)), 1e-6
  ))
  # ... while a value at the level of rounding error against the largest in the
  # matrix carries no information and is not compared.
  expect_true(basis7:::rel_close(
    cbind(c(1, 1e-14)), cbind(c(1, 5e-14)), 1e-6
  ))
  expect_true(basis7:::rel_close(matrix(1e-9), matrix(1e-9), 1e-6))
})
