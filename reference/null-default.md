# The Left Value Unless It Is NULL

Returns `a` unless it is `NULL`, in which case `b`. Written here instead
of imported: of the toolkit this package imports numericals7 only, and
the operator reached base R only in 4.4.0, later than the version
`DESCRIPTION` requires.

## Usage

``` r
a %||% b
```

## Arguments

- a, b:

  Any two objects.

## Value

`a`, or `b` where `a` is `NULL`.
