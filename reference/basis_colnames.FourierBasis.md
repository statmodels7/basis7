# Column Names of a Fourier Basis

Names the columns `const`, `sin1`, `cos1`, `sin2`, `cos2`, and so on:
the constant first, then the sine and cosine of each frequency in the
order the matrix holds them, so the name of a coefficient gives its
harmonic. The default of
[`basis_colnames.basis()`](https://statmodels7.github.io/basis7/reference/basis_colnames.basis.md)
would be `fo1` to `fo5`.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A character vector of length `basis@dimension`, `"const"` first.

## Details

At `n_pairs == 0`, which is `dimension = 1`, the result is the single
name `"const"` and the general branch is skipped. Falling through it
would give `"sin"` and `"cos"` with no number,
[`paste0()`](https://rdrr.io/r/base/paste.html) recycling a zero-length
argument to the empty string.

## See also

[`basis_colnames()`](https://statmodels7.github.io/basis7/reference/basis_colnames.md)
for the generic.
