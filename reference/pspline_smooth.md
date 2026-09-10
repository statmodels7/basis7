# A P-spline Smoother

The Eilers-Marx smoother: a B-spline basis of `k` functions over equally
spaced knots, penalized by the sum of squared `diff`-th differences of
its coefficients rather than by an integrated squared derivative. A rich
basis and a cheap penalty, which is the construction's own argument: `k`
is chosen large enough not to matter and the smoothing parameter does
the rest.

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

  The directions the smooth is made orthogonal to, `NULL` for the null
  space of the penalty.

- null_space:

  What becomes of the directions the penalty does not see.

- reparam:

  The coordinates the coefficients live in.

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

That is also the reason the argument belongs to this family and not to
every one. A difference penalty reads the coefficients as an ordered
sequence in which neighbours are comparable, which a B-spline's are and
a Fourier basis's are not – there "adjacent" is a sine, a cosine and the
next sine, and their difference means nothing.

## What it contracts to

\\D_d c = 0\\ exactly when the coefficients are a polynomial of degree
below \\d\\ in their index, and the null space of the roughness matrix
has dimension exactly `diff`: measured at `k = 20`, `degree = 3`, it is
1, 2 and 3 at `diff` of 1, 2 and 3.

That null space is only APPROXIMATELY the polynomials, and the reason is
the boundary knots. Marsden's identity gives \\\sum_j \xi_j B_j(x) = x\\
with \\\xi_j\\ the Greville abscissae, so coefficients affine in the
index give a straight line exactly where \\\xi_j\\ is itself affine in
\\j\\ – which fails at the ends of a clamped sequence, whose boundary
knots are repeated. Measured, the \\R^2\\ of \\\xi_j\\ against \\j\\ is
0.9893, 0.9979 and 0.9997 at `k` of 10, 20 and 40, the departure being a
fixed number of knots out of `k`; and the functions spanning the null
space are the polynomials of degree below `diff` to an \\R^2\\ of
1.0000000000, 0.9994 and 0.9957.

⚠️ It costs the construction nothing, which is the measurement that
matters rather than the one above.
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
constrains the block against the exact polynomials, so the
Demmler-Reinsch rotation runs on their complement, where the difference
penalty is positive definite: the built penalty is the identity to
1.0000000000 on every one of its 23 penalized directions, exactly as
[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)'s
is, and a strongly penalized fit contracts to a straight line with an
\\R^2\\ against \\(1, x)\\ of 1.0000000000 at \\\lambda = 10^{10}\\ for
both.

## Against the integrated penalty on the same basis

The two are different penalties and neither contains the other. Measured
at `k = 25`, `degree = 3` over 300 observations, the raw roughness
matrices correlate at 0.5075 and their scales differ by four orders – 6
against 2.556e+05 – the difference operator carrying no factor of the
knot spacing.

⚠️ The smoothing parameters nevertheless mean the same thing, and that
is the reparametrization doing what it is for. At matched effective
degrees of freedom of 5, 8 and 12 the two smoothing parameters stand in
a ratio of 1.0, 0.9 and 0.8, not the four orders the raw matrices differ
by, because after the Demmler-Reinsch rotation both penalties are the
identity. The fitted functions differ by a root mean square of 0.0135,
0.0243 and 0.0180 against a signal whose own standard deviation is
0.7061.

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
