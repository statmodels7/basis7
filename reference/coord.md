# One Coordinate of the Evaluation Points

The `j`th variable of the points, and the points with that variable
replaced. A basis of one variable has a vector of points and no
coordinate to pick, so both are the identity there.

## Usage

``` r
coord(x, j)

replace_coord(x, j, z)
```

## Arguments

- x:

  A numeric vector or matrix of evaluation points.

- j:

  The coordinate.

- z:

  The replacement values.

## Value

A numeric vector for `coord`, and points of the same shape as `x` for
`replace_coord`.
