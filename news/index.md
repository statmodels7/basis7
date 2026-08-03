# Changelog

## basis7 0.1.0

First release.

- The `basis` class and four generics:
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md),
  [`basis_deriv()`](https://statmodels7.github.io/basis7/reference/basis_deriv.md),
  [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md).
  Derivative order is an argument, so no order is privileged, and
  [`basis_gram()`](https://statmodels7.github.io/basis7/reference/basis_gram.md)
  returns the inner products a roughness penalty integrates.

- Two families with exact formulas:
  [`bspline_basis()`](https://statmodels7.github.io/basis7/reference/bspline_basis.md),
  whose Gram matrix is integrated exactly knot interval by knot
  interval, and
  [`fourier_basis()`](https://statmodels7.github.io/basis7/reference/fourier_basis.md),
  whose derivatives of every order and whose antiderivative come from
  one phase-shift identity, and whose Gram matrix is diagonal in closed
  form over a whole period.

- Numerical fallbacks registered on the base class, so a basis that
  implements only
  [`basis_eval()`](https://statmodels7.github.io/basis7/reference/basis_eval.md)
  is complete. Derivatives use one finite-difference stencil built from
  a Vandermonde solve, never a chain of first differences, and switch to
  a one-sided stencil at the interval endpoints.

- [`basis_int()`](https://statmodels7.github.io/basis7/reference/basis_int.md)
  is anchored: its value at the lower endpoint is exactly zero, for
  every family. Any antiderivative satisfies a differentiation check, so
  the constant is fixed by the contract instead.

- [`check_basis()`](https://statmodels7.github.io/basis7/reference/check_basis.md)
  runs six numerical checks and reports a quantity that came from a
  fallback as unchecked rather than as passed, and a property the family
  never claimed as not claimed.

- [`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
  answers the same question programmatically.

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods; the
  plot draws the basis, any derivative, or the integral.
