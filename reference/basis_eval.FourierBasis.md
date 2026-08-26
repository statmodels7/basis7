# Evaluate a Fourier Basis

Evaluates the constant and the sine-cosine pairs at the given points,
from the phase-shift identity of
[FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
taken at order zero, where it is the sinusoid itself. The first column
is `1` at every point; column \\2j\\ is \\\sin(j z)\\ and column \\2j +
1\\ is \\\cos(j z)\\, with \\z = 2\pi (x - \ell)/\omega\\.

## Arguments

- basis:

  A
  [FourierBasis](https://statmodels7.github.io/basis7/reference/FourierBasis.md)
  object.

- x:

  A numeric vector of evaluation points inside the basis interval.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric matrix with `length(x)` rows and `basis@dimension` columns,
with column names `const`, `sin1`, `cos1`, and so on.

## Details

Cost is two [`sin()`](https://rdrr.io/r/base/Trig.html) and
[`cos()`](https://rdrr.io/r/base/Trig.html) calls per frequency per
point, with no recurrence and no accumulation, so a high frequency is
evaluated as accurately as a low one. Missing points give missing rows.

## See also

[`fourier_trig()`](https://statmodels7.github.io/basis7/reference/fourier_trig.md),
which builds the trigonometric columns;
[`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
for the generic.
