# Plot a Basis

[`autoplot()`](https://statmodels7.github.io/basis7/reference/autoplot.md)
returns a graph of all `@dimension` functions of a basis on one panel,
over an equally spaced grid covering the whole interval, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints that
graph. `order` selects what is drawn: the functions themselves, a
derivative of any order, or the integral anchored at the lower endpoint.
One line is drawn per basis function, each in its own colour. The graph
has no legend, and its title names the basis and the number of its
functions.

## Usage

``` r
# S3 method for class 'basis'
plot(x, order = 0L, n = 200L, ...)

autoplot_basis(object, order = 0L, n = 200L, ...)
```

## Arguments

- x, object:

  A basis object of one variable, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md). A
  basis of several variables signals an error.

- order:

  What to draw. `0`, the default, draws the basis functions; a positive
  whole number draws that derivative; `-1` draws the integral from the
  lower endpoint. Any other value signals an error, including a
  fraction, a value below `-1` and a vector of length other than one.

- n:

  The number of grid points, default `200`. A basis with many knots or a
  high frequency needs more, because 200 points leave a curve visibly
  polygonal there.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/basis7/reference/autoplot.md)
returns a ggplot2 object;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## What is drawn

The grid is `n` equally spaced points from `@lower` to `@upper`,
endpoints included, and the curves are whichever of
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
or
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
is the one that `order` names. The vertical axis is labeled to match:
\\B(x)\\ at order `0`, \\B'(x)\\ at order `1`, \\B^{(k)}(x)\\ above
that, and \\\int B(t)\\\mathrm{d}t\\ for the integral.

The graph is drawn with ggplot2, which is suggested; the methods signal
an error that names it when it is not installed. The title, the labels
of the axes and the colours are changed on the returned object, for
instance with `+ ggplot2::labs(title = "my title")`.

## Only one variable

A basis of several variables signals an error. A product of two bases is
a surface over a rectangle; a margin can be plotted instead, which is
`tb@marginals[[1]]` for a
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md).

## Orders above the smoothness of the family

An order above the smoothness of the family is accepted and draws a flat
line at zero, because the fourth derivative of a piecewise cubic is zero
away from the knots.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
[`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md)
and
[`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
for the numbers behind the three cases of `order`, and
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md)
for the object's summary.

## Examples

``` r
b <- bspline_basis(dimension = 6)

# The three graphs show the six cubic B-splines, their first derivative,
# and their integrals.
plot(b)

plot(b, order = 1)

plot(b, order = -1)


# Every curve in the third graph starts at zero, because the integral is
# anchored at the lower endpoint.
basis_int(b, 0)
#>      bs1 bs2 bs3 bs4 bs5 bs6
#> [1,]   0   0   0   0   0   0

# The returned graph is modified like any other ggplot2 object.
p <- plot(b)

p + ggplot2::labs(title = "six cubic B-splines on [0, 1]")


# A Fourier basis at a high frequency needs a finer grid than the default.
plot(fourier_basis(dimension = 21), n = 1000)
```
