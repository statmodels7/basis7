# A P-spline Smoother

The Eilers-Marx smoother: a B-spline basis of `k` functions over equally
spaced knots, penalized by the sum of squared `diff`-th differences of
its coefficients instead of an integrated squared derivative. The basis
is rich and the penalty cheap to compute: `k` is chosen large enough not
to limit the fit, and the smoothing parameter controls the roughness.

## Usage

``` r
pspline_smooth(
  k = 20,
  degree = 3,
  diff = 2,
  constrain = NULL,
  null_space = "keep",
  reparam = "dr",
  penalty = NULL,
  lower = NULL,
  upper = NULL
)
```

## Arguments

- k:

  The number of basis functions, a whole number greater than `degree`
  and greater than `diff`.

- degree:

  The degree of the B-spline, `3` for a cubic.

- diff:

  The order of difference the penalty takes, `2` for the usual
  construction. The fit contracts to a polynomial of degree `diff - 1`.

- constrain:

  The directions to which the smooth is made orthogonal, `NULL` for the
  null space of the penalty.

- null_space:

  What happens to the directions that the penalty does not see.

- reparam:

  The coordinates in which the coefficients are expressed.

- penalty:

  `NULL` for the quadratic difference penalty, or a factory building a
  penalty from a coefficient count. See the section on the smoother's
  own page.

- lower, upper:

  The interval, `NULL` to read it from the data.

## Value

An S7 object of class
[PsplineSmoother](https://statmodels7.github.io/basis7/reference/PsplineSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## The penalty is a functional of the coefficients

Where
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
integrates \\f^{(m)}\\ against a measure, this penalizes \\\sum_j
(\Delta^d c_j)^2\\, so the roughness matrix is \\D_d^\top D_d\\ with
\\D_d\\ the difference operator on the coefficient vector. Nothing is
integrated, which is why the family has no `measure` and no `order`:
there is no measure to integrate against and no derivative whose order
to name.

The argument `diff` belongs to this family alone. A difference penalty
reads the coefficients as an ordered sequence in which neighbors are
comparable. This holds for the coefficients of a B-spline and not for
those of a Fourier basis, where adjacent coefficients are a sine, a
cosine and the next sine, and their difference has no meaning.

## The null space

\\D_d c = 0\\ exactly when the coefficients are a polynomial of degree
below \\d\\ in their index, so the null space of the roughness matrix
has dimension `diff`.

The differences are taken on the coefficients of the basis of Eilers and
Marx, whose knots are equally spaced with the same step beyond the
interval as inside it, and the penalty is carried onto the clamped basis
the package evaluates (the two span the same splines on the interval, so
only the coordinates change). On those knots the Greville abscissae are
equally spaced, and Marsden's identity makes the null space exactly the
polynomials of degree below `diff`. It is also the penalty of mgcv's
`bs = "ps"` smooths, up to the normalization that mgcv applies to its
penalty matrices.

## Comparison with the integrated penalty on the same basis

The two are different penalties, and neither contains the other. The
difference operator carries no factor of the knot spacing, so the raw
P-spline roughness matrix is on a scale several orders of magnitude
smaller than the integrated-derivative matrix on the same basis. After
the Demmler-Reinsch rotation both penalties are the identity on the
penalized directions, so equal smoothing parameters give similar amounts
of smoothing.

## References

Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
B-splines and penalties. *Statistical Science* 11(2), 89-121.

## See also

[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the integrated-derivative penalty on the same basis.

## Examples

``` r
set.seed(3)
x <- sort(runif(200))
out <- smoother_build(pspline_smooth(k = 20), x)
dim(out$X)
#> [1] 200  19

# the roughness matrix is a difference operator, so it has no measure
try(pspline_smooth(k = 20, measure = "empirical"))
#> Error in pspline_smooth(k = 20, measure = "empirical") : 
#>   unused argument (measure = "empirical")

# and 'diff' may not reach the number of functions
try(pspline_smooth(k = 4, degree = 3, diff = 4))
#> Error : 'diff' (4) must be smaller than 'k' (4): the 4-th difference of 4
#>   coefficients has no rows, so the penalty would be the zero matrix.
```
