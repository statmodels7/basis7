# Gram Matrix Against the Empirical Measure

\\B^{(d)\top} B^{(d)} / n\\ at the given points: the inner products a
design matrix produces, rather than those of the functions on their
interval.

## Usage

``` r
empirical_gram(basis, order, at)
```

## Arguments

- basis:

  An object inheriting from class `basis`.

- order:

  The derivative order.

- at:

  A numeric vector of points.

## Value

A symmetric numeric matrix with `basis@dimension` rows and columns.
