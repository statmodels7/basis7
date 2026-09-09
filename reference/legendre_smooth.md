# A Legendre Smoother

The orthogonal polynomial smoother: `k` shifted Legendre polynomials
over an interval, penalized by the integrated squared derivative of
order `order`, rotated to the Demmler-Reinsch coordinates. It is a
global basis where a B-spline is local: every function is supported on
the whole interval, so a feature at one end moves the fit at the other.

## Usage

``` r
legendre_smooth(
  k = 8,
  order = 2,
  measure = "lebesgue",
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

  The number of polynomials, a whole number of at least 2. The highest
  degree is `k - 1`.

- order:

  The order of derivative the penalty integrates, at most `k - 1`.

- measure:

  The measure the roughness is integrated against.

- constrain:

  The directions the smooth is made orthogonal to, `NULL` for the null
  space of the penalty.

- null_space:

  What becomes of the directions the penalty does not see.

- reparam:

  The coordinates the coefficients live in.

- penalty:

  `NULL` for the quadratic roughness penalty, or a factory.

- lower, upper:

  The interval, `NULL` to read it from the data.

## Value

An S7 object of class
[LegendreSmoother](https://statmodels7.github.io/basis7/reference/LegendreSmoother.md),
inheriting from
[smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## Details

The null space of the roughness matrix is the polynomials of degree
below `order`, exactly as for a B-spline, so `order` and `constrain`
mean what they mean there and `null_space = "keep"` restores `order - 1`
free columns. Measured, the null space of the order-2 Gram matrix of an
eight-function Legendre basis is spanned by the constant and the linear
function to an R-squared of 1.0000000000.

## See also

[`bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.md)
for the local family,
[`poly_basis()`](https://statmodels7.github.io/basis7/reference/poly_basis.md)
for the basis alone.

## Examples

``` r
set.seed(3)
x <- sort(runif(200, -1, 1))
out <- smoother_build(legendre_smooth(k = 8), x)
dim(out$X)
#> [1] 200   7
out$names
#> [1] "lin" "z1"  "z2"  "z3"  "z4"  "z5"  "z6" 

# 'order' may not exceed what the basis can differentiate away.
try(legendre_smooth(k = 3, order = 4))
#> Error : 'order' (4) exceeds the highest degree the basis carries (2): the
#>   derivative of that order is zero, so the penalty would be the zero matrix.
```
