# The Left Value Unless It Is NULL

Returns `a` unless it is `NULL`, in which case `b`. The function is
written here instead of imported, because this package imports only
numericals7 from the toolkit, and the operator reached base R only in
4.4.0, later than the version that `DESCRIPTION` requires.

## Usage

``` r
a %||% b
```

## Arguments

- a, b:

  Any two objects.

## Value

`a`, or `b` where `a` is `NULL`.
