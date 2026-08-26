# How Many Variables a Basis Takes

Returns the number of variables a basis is a function of: `1` for the
three shipped families and for any wrapper over them, and the number of
margins for a
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md).
It is the number of columns
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
expects its `x` to have, so it answers what shape of input the object
takes.

## Usage

``` r
basis_nvar(basis)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

## Value

A single positive integer, of length one and never `NA`.

## Details

The count is `length(basis@lower)`, the endpoints carrying one entry per
variable. A basis therefore declares its input dimension by
construction: there is no separate property that could disagree with the
interval, and the validator already requires `@lower` and `@upper` to be
of the same length.

At `basis_nvar(b) == 1` the evaluation points are a plain numeric vector
and the generics return an `n` by `@dimension` matrix. Above one they
are a matrix of `basis_nvar(b)` columns, one per variable; a plain
vector is reshaped by row, so `c(0.1, 0.2, 0.5, 0.6)` on a two-variable
basis is the two points `(0.1, 0.2)` and `(0.5, 0.6)`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) refuses a basis
of more than one variable, having no single picture to draw.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
whose input shape this describes;
[`tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.md),
the only shipped family that answers more than `1`;
[`print.basis()`](https://statmodels7.github.io/basis7/reference/print.basis.md),
which shows it on the `Variables:` line.

## Examples

``` r
# One variable for every family and every wrapper over one.
basis_nvar(bspline_basis(dimension = 5))
#> [1] 1
basis_nvar(orthonorm_basis(bspline_basis(dimension = 5)))
#> [1] 1

# A product answers with its number of margins.
tb <- tensor_basis(bspline_basis(dimension = 4), fourier_basis(dimension = 3))
basis_nvar(tb)
#> [1] 2

# Which is the number of columns basis_eval() reads: four numbers are two
# points on a two-variable basis, one row each.
dim(basis_eval(tb, c(0.1, 0.2, 0.5, 0.6)))
#> [1]  2 12
dim(basis_eval(tb, cbind(c(0.1, 0.5), c(0.2, 0.6))))
#> [1]  2 12
```
