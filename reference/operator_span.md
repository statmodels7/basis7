# The Constraint and the Free Columns of an Operator's Null Space

Splits the null space of a general operator into what is constrained
away and what is restored as free columns: the constraint is the whole
null space, and the free columns are all of it except the constant.

## Usage

``` r
operator_span(op, x)
```

## Arguments

- op:

  A
  [LinearOperator](https://statmodels7.github.io/basis7/reference/LinearOperator.md),
  with its period resolved.

- x:

  The covariate, a numeric vector.

## Value

A list of `constraint`, `free` and `params`, as
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md)
returns.

## Details

The rule is the one the polynomial families have always followed, read
for an arbitrary operator. The penalized part is made orthogonal to
every direction the penalty does not see, because a direction that is
neither penalized nor constrained is one the pencil cannot separate; and
the constant is not restored, a model carrying an intercept already
spanning it.

Each free column is divided by its root mean square over the covariate
and is **not centered**. Scaling by a positive constant is what a badly
scaled column such as \\t^3\\ over \\\[0, 365\]\\ needs, and it
preserves every property the function has; centering does not. A
periodic family restores sines and cosines here, and subtracting a
constant from a sine gives a column that is no longer periodic, which is
the property the basis was chosen for.

## See also

[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md),
which calls it, and
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
for the functions involved.
