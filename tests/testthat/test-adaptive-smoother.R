test_that("an adaptive smoother carries one component per smoothing parameter", {
  set.seed(41)
  x <- sort(runif(300))
  sm <- adaptive_smooth(k = 30, degree = 3, diff = 2, m = 4,
                        lower = 0, upper = 1)
  expect_true(S7::S7_inherits(sm, smoother))
  expect_true(S7::S7_inherits(sm, AdaptiveSmoother))
  expect_identical(sm@dimension, 30L)
  expect_identical(sm@m, 4L)
  expect_identical(sm@diff, 2L)
  # ORDER RECORDS THE DIFFERENCE ORDER, which is what smoother_span() reads
  expect_identical(sm@order, 2L)

  # the basis is an ordinary B-spline, the same object the sibling families
  # build, so nothing about the basis is particular to this construction
  b <- smoother_basis(sm, x)
  bb <- smoother_basis(pspline_smooth(k = 30, degree = 3, lower = 0, upper = 1), x)
  expect_identical(basis_eval(b, x), basis_eval(bb, x))

  g <- smoother_gram(sm, b, x)
  expect_true(is.list(g))
  expect_length(g, 4L)
  expect_identical(dim(g[[1L]]), c(30L, 30L))
})

test_that("the components sum to the difference penalty exactly", {
  # THE IDENTITY THAT MAKES THIS FAMILY A P-SPLINE AT EQUAL SMOOTHING
  # PARAMETERS. The weight functions are a B-spline basis, hence a partition
  # of unity, so sum_i diag(v_i) is the identity and sum_i S_i is D'D. It
  # holds for every m, including the two smallest, where the weight basis
  # falls to a quadratic and then to a linear one -- and where mgcv's own
  # two-component case, cbind(1, index), is NOT a partition of unity.
  x <- sort(runif(150))
  k <- 24L
  for (m in c(2L, 3L, 5L, 9L)) {
    for (dd in 1:3) {
      sm <- adaptive_smooth(k = k, diff = dd, m = m, lower = 0, upper = 1)
      g <- smoother_gram(sm, smoother_basis(sm, x), x)
      expect_length(g, m)
      d <- base::diff(diag(1, k), differences = dd)
      tot <- Reduce(`+`, g)
      expect_lt(max(abs(tot - crossprod(d))), 1e-10 * max(abs(crossprod(d))))
    }
  }
})

test_that("one component alone is not the difference penalty", {
  # the negative control the test above needs: without it, components that
  # were all D'D / m would satisfy the sum and mean nothing
  x <- sort(runif(150))
  k <- 24L
  sm <- adaptive_smooth(k = k, diff = 2, m = 5, lower = 0, upper = 1)
  g <- smoother_gram(sm, smoother_basis(sm, x), x)
  d <- crossprod(base::diff(diag(1, k), differences = 2L))
  for (i in seq_along(g)) {
    expect_gt(max(abs(g[[i]] * length(g) - d)), 1e-3 * max(abs(d)))
  }
  # and they are genuinely LOCALIZED: each is nearly all null space, its
  # weight vanishing off the support of its own weight function, while the
  # sum is not
  nulls <- vapply(g, function(P) {
    ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    sum(ev <= 1e-10 * max(ev))
  }, integer(1))
  expect_true(all(nulls > 2L))
  ev <- eigen(Reduce(`+`, g), symmetric = TRUE, only.values = TRUE)$values
  expect_identical(sum(ev <= 1e-10 * max(ev)), 2L)
})

test_that("every component is positive semidefinite", {
  # what additive_penalty() requires of each of its arguments, asserted here
  # because basis7 cannot name penalties7 to ask it
  x <- sort(runif(150))
  sm <- adaptive_smooth(k = 30, m = 5, lower = 0, upper = 1)
  for (P in smoother_gram(sm, smoother_basis(sm, x), x)) {
    expect_lt(max(abs(P - t(P))), 1e-12 * max(abs(P)))
    ev <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    expect_gt(max(ev), 0)
    expect_gt(min(ev), -1e-10 * max(ev))
  }
})

test_that("smoother_build carries the components through", {
  set.seed(42)
  x <- sort(runif(300))
  sm <- adaptive_smooth(k = 30, m = 4, lower = 0, upper = 1)
  out <- smoother_build(sm, x)
  expect_true(is.list(out$S))
  expect_length(out$S, 4L)
  expect_identical(dim(out$S[[1L]]), c(ncol(out$X), ncol(out$X)))
  expect_null(dimnames(out$S[[1L]]))
  # the free column is read off the SUM, a coordinate being unpenalized only
  # where NO component touches it. ⚠️ For this family the two readings
  # coincide -- the built components are 93 per cent nonzero, localized in
  # rank rather than in sparsity -- so this asserts the answer and not the
  # route, and an injection reading one component instead stays green. The
  # measurement is in smoother_build()'s own comment.
  expect_identical(out$unpenalized, 1L)
  expect_identical(out$names[[1L]], "lin")
  expect_identical(
    vapply(out$S, function(P) as.integer(sum(cumprod(colSums(abs(P)) == 0))),
           integer(1)),
    rep(1L, length(out$S))
  )

  # AND THE BLOCK IS THE P-SPLINE'S, bit for bit: the two families differ in
  # the penalty and in nothing else, so at equal smoothing parameters they
  # are one model
  ps <- smoother_build(pspline_smooth(k = 30, reparam = "none",
                                      lower = 0, upper = 1), x)
  expect_identical(out$X, ps$X)
  expect_lt(max(abs(Reduce(`+`, out$S) - ps$S)), 1e-10 * max(abs(ps$S)))
})

test_that("the block is reapplied at new values rather than rebuilt", {
  set.seed(43)
  x <- sort(runif(200))
  for (rp in c("none", "orthonorm")) {
    sm <- adaptive_smooth(k = 25, m = 4, reparam = rp)
    out <- smoother_build(sm, x)
    i <- c(2L, 37L, 111L, 200L)
    expect_reapplied(smoother_apply(sm, out$blueprint, x[i]), out$X[i, ])
  }
})

test_that("a rotation is not a change of model", {
  # the coordinates differ and the fitted values cannot: the transform is
  # invertible, so the two span the same space
  set.seed(44)
  x <- sort(runif(250))
  y <- sin(2 * pi * x) + rnorm(250, sd = 0.2)
  lam <- c(3, 0.01, 10, 0.05)
  fit <- function(rp) {
    o <- smoother_build(adaptive_smooth(k = 25, m = 4, reparam = rp,
                                        lower = 0, upper = 1), x)
    K <- crossprod(o$X) + Reduce(`+`, Map(function(P, l) l * P, o$S, lam))
    o$X %*% solve(K, crossprod(o$X, y))
  }
  a <- fit("none")
  b <- fit("orthonorm")
  expect_gt(sd(a), 0)
  expect_lt(max(abs(a - b)), 1e-8 * sd(a))
})

test_that("the null space may be kept or dropped", {
  set.seed(45)
  x <- sort(runif(200))
  keep <- smoother_build(adaptive_smooth(k = 25, m = 4), x)
  drop <- smoother_build(adaptive_smooth(k = 25, m = 4, null_space = "drop"), x)
  expect_identical(ncol(keep$X), ncol(drop$X) + 1L)
  expect_identical(keep$unpenalized, 1L)
  expect_identical(drop$unpenalized, 0L)
  # at order 3 the null space is the lines AND the parabolas, so two
  # columns are free rather than one
  three <- smoother_build(adaptive_smooth(k = 25, m = 4, diff = 3), x)
  expect_identical(three$unpenalized, 2L)
})

test_that("reparam = \"dr\" is rejected, and by the generic as well", {
  expect_error(adaptive_smooth(k = 30, m = 4, reparam = "dr"),
               "not available for an adaptive smoother")
  # THE REFUSAL IS ALSO IN smoother_reparam(), which is what makes it hold
  # for a family written later that answers smoother_gram() with a list.
  # Reached with a smoother whose own constructor admits "dr".
  x <- sort(runif(80))
  sm <- bspline_smooth(k = 10, reparam = "dr", lower = 0, upper = 1)
  b <- smoother_basis(sm, x)
  sp <- smoother_span(sm, x)
  g <- smoother_gram(sm, b, x)
  expect_error(smoother_reparam(sm, b, x, list(g, g), sp$constraint),
               "needs one roughness matrix")
  # and the same call with ONE matrix is fine, so the error is about the
  # list and not about the arguments
  expect_silent(smoother_reparam(sm, b, x, g, sp$constraint))
})

test_that("the arguments a sum of components has no reading for are rejected", {
  expect_error(adaptive_smooth(k = 30, m = 4, null_space = "shrink"),
               "not available for an adaptive smoother")
  expect_error(adaptive_smooth(k = 30, m = 4,
                               penalty = function(n_coef) n_coef),
               "'penalty' is not available")
  # a family whose penalty integrates nothing declares no measure
  expect_error(adaptive_smooth(k = 30, m = 4, measure = "empirical"),
               "unused argument")
})

test_that("m is validated against the differences it is carried on", {
  expect_error(adaptive_smooth(k = 30, m = 1), "pspline_smooth")
  expect_error(adaptive_smooth(k = 30, m = 0), "pspline_smooth")
  # the profile lives on the k - diff differences and cannot hold more
  # functions than there are of them: mgcv's own condition
  expect_error(adaptive_smooth(k = 10, diff = 2, m = 8), "too large for 'k'")
  expect_error(adaptive_smooth(k = 10, diff = 4, m = 6), "too large for 'k'")
  expect_silent(adaptive_smooth(k = 10, diff = 2, m = 7))
  # and the arguments the family shares with its sibling keep their checks
  expect_error(adaptive_smooth(k = 3, degree = 3, m = 2), "too small for 'degree'")
  expect_error(adaptive_smooth(k = 30, diff = 30, m = 4), "must be smaller than")
})

test_that("print says what the penalty is", {
  # a difference penalty integrates nothing, so naming a measure would
  # report a construction the smoother does not run
  out <- utils::capture.output(print(adaptive_smooth(k = 30, m = 5)))
  expect_true(any(grepl("weighted by 5 components", out, fixed = TRUE)))
  expect_false(any(grepl("measure", out, fixed = TRUE)))
  ps <- utils::capture.output(print(pspline_smooth(k = 20)))
  expect_true(any(grepl("difference of order 2", ps, fixed = TRUE)))
  expect_false(any(grepl("measure", ps, fixed = TRUE)))
  # and a family that DOES integrate still says so
  bs <- utils::capture.output(print(bspline_smooth(k = 20)))
  expect_true(any(grepl("derivative of order 2", bs, fixed = TRUE)))
  expect_true(any(grepl("lebesgue measure", bs, fixed = TRUE)))
})
