# The Left Value Unless It Is NULL

Returns `a` unless it is `NULL`, in which case `b`. Written here rather
than imported: this package depends on numericals7 alone, and the
operator reached base R only in 4.4.0, later than the version
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
