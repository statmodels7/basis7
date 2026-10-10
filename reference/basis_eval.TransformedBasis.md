# Evaluate a Transformed Basis

Evaluates the parent and multiplies by the transform, \\\tilde{B}(x) =
B(x)\\T\\, with no further computation. The result spans a subspace of
what the parent spans, and all of it when \\T\\ is square and
invertible. Properties of the individual columns, such as the partition
of unity of a B-spline, are in general not kept.

## Arguments

- basis:

  A
  [TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval, which
  is the parent's.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names `on1`, `cn1` or `dr1` and so on.

## Details

Cost is one parent evaluation plus one `length(x)` by `K` by `m` matrix
product. Because
[`new_transformed()`](https://statmodels7.github.io/basis7/reference/new_transformed.md)
collapses a chain into one matrix, that is the cost however many
transformations were composed.

## See also

[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
for the generic;
[TransformedBasis](https://statmodels7.github.io/basis7/reference/TransformedBasis.md)
for the class.
