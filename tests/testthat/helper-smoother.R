# Helpers for the smoother tests.

# The rows a blueprint reapplies, against the rows of the block they came
# from.
#
# NOT AN IDENTITY, and the reason is measured rather than assumed. Both
# routes end in the same product on the SAME basis object -- the parent's
# design times the transform, which is what smoother_build() means when its
# own comment says the two are the same arithmetic rather than the same
# formula. What differs is the SHAPE: the build multiplies a matrix with one
# row per observation and the reapplication one with a row per new point,
# and a BLAS is free to block, vectorize and accumulate a product
# differently by shape.
#
# Measured 2026-09-09: bit-identical under the reference BLAS here and under
# R release and oldrel-1 on the CI's Ubuntu image, and two ulps apart under
# R-devel on that same image, where it separated three of nine numbers of
# one comparison and failed the check.
#
# THE TOLERANCE IS CHOSEN SO THE DEFECT THE CHECK EXISTS FOR STILL FAILS IT.
# That defect is a basis REBUILT from the new points rather than reapplied,
# which carries different knots, a different empirical Gram and a different
# rotation. Measured on the axes case a rebuild sits 3.025 from a block
# whose largest entry is 2.1, so the gap is of the size of the quantity
# itself; twelve orders separate it from the tolerance below, and four
# separate that tolerance from the platform's last bits.
expect_reapplied <- function(applied, built, tolerance = 1e-12) {
  testthat::expect_identical(dim(applied), dim(built))
  s <- max(abs(built))
  # a comparison against a block of zeros would pass whatever was handed to
  # it, so the scale is asserted rather than floored
  testthat::expect_gt(s, 0)
  testthat::expect_lt(max(abs(applied - built)), tolerance * s)
}
