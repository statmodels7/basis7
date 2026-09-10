# Column Names for an Operator's Restored Columns

Turns the labels
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
gives into names a design matrix can carry: `sin(0.0172 t)` becomes
`sin1`, `cos(0.0344 t)` becomes `cos2`, and anything else keeps a
syntactically safe form of its label.

## Usage

``` r
operator_free_names(labels)
```

## Arguments

- labels:

  The `label` column of
  [`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md),
  subset to the functions that were restored.

## Value

A character vector of one name per label.

## Details

The number is the rank of the frequency among those restored, not the
frequency itself, so the columns of a two-harmonic smooth read `sin1`,
`cos1`, `sin2`, `cos2` and match what a reader of a Fourier basis
expects. A column from a real root keeps its label, `t` and `t^2`
passing through [`make.names()`](https://rdrr.io/r/base/make.names.html)
as they stand.
