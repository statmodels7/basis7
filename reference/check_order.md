# Validate a Derivative Order

A single non-negative integer for a basis of one variable, and one per
variable for a basis of several.

## Usage

``` r
check_order(order, nvar = 1L)
```

## Arguments

- order:

  The value supplied by the caller.

- nvar:

  The number of variables the basis takes.

## Value

`order`, as an integer vector of length `nvar`.

## Details

A multi-index is required rather than recycled from a scalar, because a
scalar has two plausible readings for a product basis – that order in
every coordinate, or that total order – and guessing between them would
be a silent choice. The single exception is zero, which means no
derivative under either reading.
