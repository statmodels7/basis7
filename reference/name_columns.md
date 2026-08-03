# Name the Columns of a Basis Matrix

Gives a basis matrix the column names the basis declares, so that
evaluation, derivatives and integrals of the same basis always agree on
them.

## Usage

``` r
name_columns(m, basis)
```

## Arguments

- m:

  A numeric matrix with as many columns as the basis has functions.

- basis:

  An object inheriting from class `basis`.

## Value

`m`, with column names set.
