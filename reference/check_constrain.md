# Check a Polynomial Family's Constraint Argument

Validates `constrain` for a family whose null space is the polynomials
of degree below `order`. The constraint must **contain** that null
space, so a degree below `order - 1` is rejected.

## Usage

``` r
check_constrain(constrain, op)
```

## Arguments

- constrain:

  The value given, `NULL` or a whole number.

- op:

  The penalty's operator.

## Value

`constrain` as a length-one integer, or `NULL`.

## Details

The requirement follows from identifiability. A direction that the
penalty does not see and that the constraint does not remove is neither
penalized nor identified, and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
signals an error there: with `order = 3` and a constraint spanning only
the constant and the linear function, the quadratic direction is left
free and unpenalized. The error is reported here, at the constructor.

Constraining **beyond** the null space is admitted. It makes the block
orthogonal, over the observed covariate, to the polynomials up to degree
`constrain`, which avoids collinearity with a parametric term of that
degree written in the formula. The cost is one dimension per degree.
