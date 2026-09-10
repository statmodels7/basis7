# The Names of the Free Columns

Names the columns `null_space = "keep"` restores where the null space is
the polynomials: `lin` for the linear one, which is the only one at
`order = 2`, and `poly2`, `poly3` and so on for the higher powers an
order above 2 leaves unpenalized.

## Usage

``` r
free_names(n)
```

## Arguments

- n:

  How many columns were restored.

## Value

A character vector of length `n`.

## Details

A smoother whose null space is not the polynomials names its own columns
instead, through the `free_names` element of
[`smoother_span()`](https://statmodels7.github.io/basis7/reference/smoother_span.md);
see
[`operator_free_names()`](https://statmodels7.github.io/basis7/reference/operator_free_names.md),
which a periodic null space uses.
