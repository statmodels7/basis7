# Check a Polynomial Family's Constraint Argument

Validates `constrain` for a family whose null space is the polynomials
of degree below `order`. The constraint must **contain** that null
space, so a degree below `order - 1` is rejected.

## Usage

``` r
check_constrain(constrain, order)
```

## Arguments

- constrain:

  The value given, `NULL` or a whole number.

- order:

  The penalty's derivative order.

## Value

`constrain` as a length-one integer, or `NULL`.

## Details

The requirement is not a convention. A direction the penalty does not
see and the constraint does not remove is neither penalized nor
identified, and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
signals an error there: with `order = 3` and a constraint spanning only
the constant and the linear function, the quadratic direction is left
free and unpenalized. Reported here, where the two numbers were written,
rather than several frames down.

Constraining **beyond** the null space is legitimate. It buys
orthogonality to a parametric term written in the formula. Measured on
`y ~ 1 + x + x^2 + s(x)` at `k = 20` and 400 observations, the largest
correlation between a column of the block and `x^2` falls from 0.995 at
the default to 1.9e-15 at `constrain = 2`, which is exact by
construction, and the standard error of the quadratic coefficient falls
with it by more than an order of magnitude. The cost is one dimension.
