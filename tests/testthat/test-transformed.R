# The one wrapper for every linear reparametrisation.

parents <- function() {
  list(
    bspline = bspline_basis(dimension = 7),
    fourier = fourier_basis(lower = -1, upper = 2, dimension = 5),
    legendre = poly_basis(dimension = 5)
  )
}


test_that("orthonormalisation makes the Gram matrix the identity", {
  # The exit test of the construction, by both routes: the exact congruence
  # the transform performs, and an independent quadrature that knows nothing
  # about it.
  for (nm in names(parents())) {
    o <- orthonorm_basis(parents()[[nm]])
    k <- o@dimension
    expect_equal(basis_gram(o), diag(k),
      tolerance = 1e-10, ignore_attr = TRUE, label = nm
    )
    expect_equal(basis7:::numerical_gram(o, 0L, panels = 300L, nodes = 10L),
      diag(k),
      tolerance = 1e-8, ignore_attr = TRUE, label = nm
    )
  }
})


test_that("orthonormalisation preserves the span", {
  b <- bspline_basis(dimension = 6)
  o <- orthonorm_basis(b)
  x <- seq(0, 1, length.out = 80)
  # any function the parent represents, the transform represents too
  set.seed(21)
  beta <- rnorm(6)
  target <- basis_eval(b, x) %*% beta
  fit <- stats::lm.fit(basis_eval(o, x), target)$fitted.values
  expect_equal(drop(fit), drop(target), tolerance = 1e-9, ignore_attr = TRUE)
})


test_that("orthonormalising at a derivative order is possible where it exists", {
  # At order 1 the Gram matrix of a Fourier basis is singular, because the
  # constant differentiates to zero, so the refusal must name that.
  expect_error(orthonorm_basis(fourier_basis(dimension = 5), order = 1L),
    "singular"
  )
})


test_that("a constraint removes exactly its rank from the dimension", {
  b <- bspline_basis(dimension = 8)
  x <- seq(0, 1, length.out = 150)

  c1 <- colSums(basis_eval(b, x))
  cs <- constrain_basis(b, c1)
  expect_identical(cs@dimension, 7L)
  expect_lt(max(abs(colSums(basis_eval(cs, x)))), 1e-10)

  c2 <- rbind(c1, colSums(basis_eval(b, x) * x))
  cs2 <- constrain_basis(b, c2)
  expect_identical(cs2@dimension, 6L)
  expect_lt(max(abs(crossprod(cbind(1, x), basis_eval(cs2, x)))), 1e-9)

  # a repeated constraint has rank one, and removes one function
  expect_identical(constrain_basis(b, rbind(c1, 2 * c1))@dimension, 7L)
})


test_that("an impossible or malformed constraint is refused", {
  b <- bspline_basis(dimension = 5)
  expect_error(constrain_basis(b, diag(5)), "leaves no functions")
  expect_error(constrain_basis(b, rep(1, 4)), "5 columns")
  expect_error(constrain_basis(b, c(NA, 1, 1, 1, 1)), "missing values")
})


test_that("transforms compose by multiplication rather than by nesting", {
  b <- bspline_basis(dimension = 6)
  x <- seq(0, 1, length.out = 40)

  o1 <- orthonorm_basis(b)
  o2 <- orthonorm_basis(o1)
  expect_true(S7::S7_inherits(o2, TransformedBasis))
  # one level deep, whatever the chain length: the parent of a transform of a
  # transform is the original basis, not the intermediate one
  expect_identical(o2@parent_basis@basis_name, "bspline")
  expect_false(S7::S7_inherits(o2@parent_basis, TransformedBasis))
  expect_identical(dim(o2@transform), c(6L, 6L))

  # orthonormalising something already orthonormal changes nothing
  expect_equal(abs(basis_eval(o2, x)), abs(basis_eval(o1, x)),
    tolerance = 1e-8
  )

  # a constraint on top of an orthonormalisation is still one object
  cs <- constrain_basis(o1, colSums(basis_eval(o1, x)))
  expect_identical(cs@parent_basis@basis_name, "bspline")
  expect_identical(dim(cs@transform), c(6L, 5L))
})


test_that("derivatives and integrals transform by the same matrix", {
  b <- poly_basis(dimension = 5)
  o <- orthonorm_basis(b)
  x <- seq(0.1, 0.9, length.out = 9)
  for (d in 1:3) {
    expect_equal(basis_deriv(o, x, order = d),
      basis_deriv(b, x, order = d) %*% o@transform,
      ignore_attr = TRUE, label = paste("order", d)
    )
  }
  expect_equal(basis_int(o, x), basis_int(b, x) %*% o@transform,
    ignore_attr = TRUE
  )
  # the anchor survives: a combination of columns all zero there is zero there
  expect_true(all(basis_int(o, o@lower) == 0))
})


test_that("a transformed basis reports the PARENT's numerical status", {
  # Its own methods delegate and multiply, so what is numerical about it is
  # whatever was numerical about the parent. Saying otherwise would let
  # check_basis compare a finite difference with a finite difference.
  Bare <- S7::new_class("BareForTransform", parent = basis)
  S7::method(basis_eval, Bare) <- function(basis, x, ...) {
    out <- cbind(1, x, x^2)
    colnames(out) <- basis_colnames(basis)
    out
  }
  bare <- Bare(
    basis_name = "bare", dimension = 3L, lower = 0, upper = 1,
    basis_params = list()
  )
  o <- orthonorm_basis(bare)
  expect_true(all(basis_is_numerical(o)))

  r <- check_basis(o, verbose = FALSE)
  expect_true(is.na(r[["deriv"]]))
  expect_true(is.na(r[["integral"]]))

  # while a transform of an exact parent stays exact
  expect_false(any(basis_is_numerical(orthonorm_basis(poly_basis(dimension = 4)))))
})


test_that("check_basis passes on transformed bases", {
  for (tb in list(
    orthonorm_basis(bspline_basis(dimension = 6)),
    constrain_basis(bspline_basis(dimension = 6), rep(1, 6)),
    orthonorm_basis(poly_basis(dimension = 5))
  )) {
    r <- check_basis(tb, verbose = FALSE)
    expect_true(all(r, na.rm = TRUE), label = tb@basis_name)
    expect_false(isFALSE(r[["deriv"]]), label = tb@basis_name)
  }
})
