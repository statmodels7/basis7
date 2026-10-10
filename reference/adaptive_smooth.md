# An Adaptive Smoother

A difference penalty whose weight varies along the covariate, so that a
function may be smoothed hard where it is quiet and left free where it
is not. Where
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
penalizes every difference alike under one smoothing parameter, this
family carries `m` of them, and the data determine how the roughness is
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

  The directions to which the smooth is made orthogonal, `NULL` for the
  null space of the penalty.

- null_space:

  What happens to the directions that the penalty does not see.
  `"shrink"` is rejected, because the shrinkage is a fraction of the
  eigenvalues of one roughness matrix and this penalty is a sum of
  components.

- reparam:

  The coordinates in which the coefficients are expressed, `"none"` or
  `"orthonorm"`. `"dr"` is rejected.

- penalty:

  Any value other than `NULL` signals an error. The penalty of this
  family is the sum of its components, and a factory would replace it
  with one penalty over the coefficients, which is
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
returns those `m` components and the layer that places the smooth
combines them into one penalty with `m` smoothing parameters. The weight
profile is itself a spline, and its coefficients are those smoothing
parameters.

The index basis is built over the range of the index, so an affine
relabeling of the index (\\1, \ldots, n\\, \\i/n\\, or the \\i/k\\ of
mgcv) gives the same profile up to rounding.

## Equal smoothing parameters

The weight functions are a B-spline basis, hence a partition of unity,
so \\\sum_i v_i = 1\\ and therefore \\\sum_i S_i = D^\top D\\, up to
rounding. Holding every \\\lambda_i\\ at one value gives the P-spline
penalty at that value, so
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md)
is a special case of this family, and the additional freedom is used
only where the data require it.

## Comparison with a single smoothing parameter

On a function whose roughness varies along the covariate the adaptive
penalty smooths the quiet stretches more strongly than a P-spline with
one smoothing parameter, and follows the noise less there. On a function
of constant roughness its additional smoothing parameters have nothing
to adapt to, and the fit is slightly worse than the P-spline fit.

The differences are taken in the Eilers-Marx coordinates of
[`pspline_smooth()`](https://statmodels7.github.io/basis7/reference/pspline_smooth.md),
as in the adaptive P-spline of mgcv's `bs = "ad"`. At equal smoothing
parameters the two sums are the same P-spline penalty; the weight
functions differ, mgcv building them from a different basis over the
index, so the two fits are close and not identical.

## The coordinates

`reparam = "dr"` is rejected: the Demmler-Reinsch rotation diagonalizes
the pencil of the Gram matrix against a single penalty, this family has
`m` penalties, and a rotation that makes one component the identity
leaves the others unspecified. The default is `"none"`, and
`"orthonorm"` is also available.

## The rank of the penalty

Each \\S_i\\ has a large null space, made of the directions on which its
weight function vanishes. The null space of the sum is the intersection
of the null spaces of the components and does not depend on the
smoothing parameters, so the rank is read from the components and not
from the assembled \\S(\lambda)\\, whose numerical rank falls as one
smoothing parameter grows and depends on the tolerance used.

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
#> [1] 1.421085e-14

# the coordinates cannot be Demmler-Reinsch: there are several pencils
try(adaptive_smooth(k = 30, m = 4, reparam = "dr"))
#> Error : reparam = "dr" is not available for an adaptive smoother: Demmler-Reinsch
#>   diagonalizes the pencil of the Gram matrix against ONE penalty, and this family
#>   carries several, so a rotation making one of them the identity leaves the
#>   others arbitrary. Use "none" (the default) or "orthonorm".
```
