# Build a Transformed Basis

Wraps a basis in a
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md),
collapsing the transform into the parent's when the parent is already
one, so that a chain of transformations is stored as a single matrix.
The one constructor
[`orthonorm_basis()`](https://statmodels7.github.io/basis7/reference/orthonorm_basis.md),
[`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
and
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
all go through.

## Usage

``` r
new_transformed(basis, transform, name, prefix, params = list())
```

## Arguments

- basis:

  The basis to transform, any object inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- transform:

  The matrix \\T\\, with one row per function of `basis`.

- name:

  The `basis_name` of the result. The callers pass `orthonorm(...)`,
  `constrained(...)` or `dr(...)` wrapped around the parent's own name,
  so the name nests even where the matrices collapse.

- prefix:

  The two-letter prefix its column names are numbered under: `"on"`,
  `"cn"` or `"dr"`.

- params:

  Extra entries for `basis_params`, merged with `prefix`.
  [`constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.md)
  adds `constraint_rank` and
  [`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
  adds that and `empirical_variance`.

## Value

An object of class
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md).

## Details

When `basis` is itself a `TransformedBasis`, the result holds
`basis@transform %*% transform` and the original's parent. The dimension
of the result is `ncol(transform)`, and its interval is the parent's, a
linear map of the functions leaving the domain alone.

Nothing is validated. The callers have already checked the shape of
their own matrix, and a `transform` whose row count disagrees with the
parent's dimension gives R's own non-conformable-arguments error at the
first evaluation.

## See also

[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md),
the class it builds.
