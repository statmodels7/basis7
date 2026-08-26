# Reject Points Outside a Range, and Clamp Those On Its Edge

Throws if any entry of `z` lies outside `[lo, hi]` by more than `1e-8`
times the width of the range, and otherwise returns `z` with anything
inside that tolerance moved onto the nearer endpoint. The one place the
range rule of the package is written;
[`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md)
calls it once per variable.

## Usage

``` r
clamp_to_range(z, lo, hi, what)
```

## Arguments

- z:

  A numeric vector. `NA` entries are ignored.

- lo, hi:

  The endpoints, single finite numbers with `lo < hi`. Not validated
  here; the class validator has already required it.

- what:

  A phrase naming the range, spliced into the error message between
  "fall outside" and the interval. Callers pass `"the basis interval"`
  or `"the range of variable 2"`.

## Value

`z`, of the same length, with entries within tolerance of an endpoint
set to that endpoint exactly.

## Details

The tolerance is relative to the width, `1e-8 * (hi - lo)`, so the rule
means the same on \\\[0, 1\]\\ and on \\\[0, 1000\]\\. Entries that are
`NA` are excluded from both the test and the clamp, and travel through.

The error message names the count and the range, as in
`2 of 5 evaluation points fall outside the basis interval [0, 1].`, and
is raised with `call. = FALSE`, so the phrase `what` supplies is what a
user sees in place of this function's name.

## See also

[`check_eval_points()`](https://statmodels7.github.io/basis7/reference/check_eval_points.md),
its only caller.
