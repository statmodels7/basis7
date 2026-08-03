# B-Spline Basis

The S7 class of B-spline bases. Constructed by
[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md).

## Usage

``` r
BsplineBasis(
  basis_name = character(0),
  dimension = integer(0),
  lower = integer(0),
  upper = integer(0),
  basis_params = list()
)
```

## Arguments

- basis_name:

  A short name for the family, used when printing.

- dimension:

  The number of functions in the basis, a positive integer.

- lower, upper:

  The endpoints of the interval the basis lives on, or one endpoint per
  variable for a basis of several.

- basis_params:

  A named list of whatever else the subclass needs.

## Value

An object of class `BsplineBasis`. Use
[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)
rather than calling the class directly, so that the knots are placed and
the arguments checked.

## Details

Evaluation, derivatives and integrals come from splines2, which computes
all three from the recurrence rather than by differencing. The Gram
matrix is integrated exactly, interval by interval: on each knot
interval the integrand is a polynomial of known degree, so a
Gauss-Legendre rule sized from that degree leaves no quadrature error.

## See also

[`bspline_basis`](https://statmodels7.github.io/basis7/reference/bspline_basis.md)

## Examples

``` r
b <- bspline_basis(dimension = 5)
S7::S7_inherits(b, BsplineBasis)
#> [1] TRUE
```
