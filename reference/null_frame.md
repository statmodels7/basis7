# Assemble the Null-Space Table

Builds the data frame
[`operator_null()`](https://statmodels7.github.io/basis7/reference/operator_null.md)
returns, writing the label of each function from its rate, frequency and
power of \\t\\.

## Usage

``` r
null_frame(rate, freq, degree, part = rep("", length(rate)))
```

## Arguments

- rate, freq, degree:

  The three descriptions, of equal length.

- part:

  Which half of a conjugate pair the row is, `"sin"` or `"cos"`, and
  `""` for a row that comes from a real root. It is carried as a column
  rather than deduced from the row's position: a real root ahead of a
  pair shifts every position after it, and a pair read off the parity
  would then come out as two sines.

## Value

A data frame of five columns.
