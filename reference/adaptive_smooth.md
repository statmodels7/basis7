# An Adaptive Smoother

A difference penalty whose weight varies along the covariate, so that a
function may be smoothed hard where it is quiet and left free where it
is not. Where
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
penalizes every difference alike under one smoothing parameter, this one
carries `m` of them and lets the data say how the roughness is
distributed.

## Usage

``` r
adaptive_smooth(
  k = 40,
  degree = 3,
  diff = 2,
  m = 5,
  constrain = NULL,
  null_space = "keep",
  reparam = "none",
  penalty = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- k:

  The number of basis functions. A rich basis is the point of the
  construction, so the default is larger than
  [`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)'s.

- degree:

  The degree of the B-spline, `3` for a cubic.

- diff:

  The order of difference the penalty takes. The fit contracts to a
  polynomial of degree `diff - 1`.

- m:

  The number of weight components, and so the number of smoothing
  parameters. mgcv calls it `m` as well, and takes the same default. Two
  gives one weight rising and one falling across the index; more give a
  finer profile at the price of a smoothing parameter each.

- constrain:

  The directions the smooth is made orthogonal to, `NULL` for the null
  space of the penalty.

- null_space:

  What becomes of the directions the penalty does not see. `"shrink"` is
  rejected here; see the note below.

- reparam:

  The coordinates the coefficients live in, `"none"` or `"orthonorm"`.
  `"dr"` is rejected.

- penalty:

  Rejected here, and accepted only to say so: the penalty of this family
  is the sum of its components, and a factory replaces it with one
  penalty over the coefficients, which is
  [`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
  with that factory.

- lower, upper:

  The interval, `NULL` to read it from the data.

## Value

An S7 object of class
[AdaptiveSmoother](https://statmodels7.github.io/basis7/reference/AdaptiveSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## The construction

Write \\D\\ for the matrix of `diff`-th differences of the coefficients.
A P-spline penalizes \\\lambda\\ c^\top D^\top D\\ c\\; this one gives
each difference a weight of its own, \$\$c^\top D^\top
\mathrm{diag}(w)\\ D\\ c, \qquad w = \sum\_{i=1}^{m} \lambda_i v_i,\$\$
where \\v_1, \ldots, v_m\\ is a B-spline basis evaluated over the
coefficient INDEX. Because \\\mathrm{diag}\\ is linear the whole penalty
is the sum \\\sum_i \lambda_i S_i\\ with \\S_i = D^\top
\mathrm{diag}(v_i) D\\, so
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md)
answers with those `m` components and whichever layer places the smooth
turns them into one penalty with `m` smoothing parameters. The weight
profile is itself a spline, and its coefficients are those smoothing
parameters.

The index basis is built over the index's own range, so an affine
relabelling of the index – \\1, \ldots, n\\, or \\i/n\\, or mgcv's
\\i/k\\ – gives the identical profile, measured to 1e-16. The
construction carries no arbitrary constant.

## At equal smoothing parameters it IS a P-spline

The weight functions are a B-spline basis, hence a partition of unity,
so \\\sum_i v_i = 1\\ and therefore \\\sum_i S_i = D^\top D\\ exactly –
measured at 8.9e-16 to 2.7e-15 for `m` from 2 to 12. Holding every
\\\lambda_i\\ at one value gives the P-spline penalty at that value, so
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
is the interior point of this family rather than a different
construction, and the extra freedom is spent only where the data pay for
it.

## What it buys, and what it costs

Measured against a single-lambda P-spline through mgcv's REML, which
shares no code with this package, on 400 observations at `k = 40` and
eight seeds: on a truth of variable roughness – a sine with a narrow
bump – the adaptive wins on 8 seeds of 8, at a median root mean square
error of 0.0357 against 0.0443, and does so at FEWER effective degrees
of freedom, 17.4 against 23.6. On a truth of constant roughness it wins
on 0 seeds of 8, 0.0315 against 0.0310: about 1.6 per cent worse, which
is what the extra smoothing parameters cost where there is nothing to
adapt to. That second measurement is the control, without which the
first would show only that more parameters fit better.

Where the gain comes from is not where it is first looked for. Split by
region on one sample, the adaptive is better on the quiet left (0.0242
against 0.0305) and on the smooth right (0.0437 against 0.0476) and
slightly WORSE at the feature itself (0.0553 against 0.0503). What it
buys is not a sharper peak but less noise chasing where the function is
quiet.

## The coordinates

`reparam = "dr"` is rejected, and by construction rather than by choice:
Demmler-Reinsch diagonalizes the pencil of the Gram matrix against a
single penalty, and here there are `m` of them, so a rotation making one
component the identity leaves the others arbitrary. The default is
`"none"`, which the measurement prefers on both axes it can be judged
on. At a spread of smoothing parameters an adaptive really reaches –
1.3e8, measured on mgcv – the condition number of the system solved is
3.4e3 in the raw coordinates against 4.7e4 orthonormalized, and the
components' scales, which decide whether the `m` smoothing parameters
are comparable with one another, spread by a factor of 1.7 against 2.8.

## The rank is the family's

Each \\S_i\\ is nearly all null space, its weight vanishing off the
support of its own weight function: measured at `k = 40`, `m = 5`, the
five components have null dimensions 20, 1, 1, 1 and 19 out of 38. The
null space of the SUM is the intersection of theirs and does not move
with the smoothing parameters, which is why the rank must be read from
the components rather than from the assembled \\S(\lambda)\\. Measured,
a rank counted off the assembled matrix reads 38, 38, 26 and 19 as one
parameter is raised through 1, 1e6, 1e12 and 1e24, where the family's
own is 38 throughout.

## References

Ruppert, D. and Carroll, R. J. (2000). Spatially-adaptive penalties for
spline fitting. *Australian and New Zealand Journal of Statistics*
42(2), 205-223.

Krivobokova, T., Crainiceanu, C. M. and Kauermann, G. (2008). Fast
adaptive penalized splines. *Journal of Computational and Graphical
Statistics* 17(1), 1-20.

## See also

[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md),
which is this family at equal smoothing parameters, and
[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md),
which returns the components.

## Examples

``` r
set.seed(4)
x <- sort(runif(300))
out <- smoother_build(adaptive_smooth(k = 30, m = 4), x)

# the penalty comes back as one component per smoothing parameter
length(out$S)
#> [1] 4
dim(out$S[[1]])
#> [1] 29 29

# and at equal smoothing parameters their sum is the P-spline penalty
ps <- smoother_build(pspline_smooth(k = 30, reparam = "none"), x)
max(abs(Reduce(`+`, out$S) - ps$S))
#> [1] 4.440892e-15

# the coordinates cannot be Demmler-Reinsch: there are several pencils
try(adaptive_smooth(k = 30, m = 4, reparam = "dr"))
#> Error : reparam = "dr" is not available for an adaptive smoother: Demmler-Reinsch
#>   diagonalizes the pencil of the Gram matrix against ONE penalty, and this family
#>   carries several, so a rotation making one of them the identity leaves the
#>   others arbitrary. Use "none" (the default) or "orthonorm".
```
