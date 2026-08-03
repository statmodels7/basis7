# basis7 0.1.0

First release.

* The `basis` class and four generics: `basis_eval()`, `basis_deriv()`,
  `basis_int()` and `basis_gram()`. Derivative order is an argument, so no
  order is privileged, and `basis_gram()` returns the inner products a
  roughness penalty integrates.

* Two families with exact formulas: `bspline_basis()`, whose Gram matrix is
  integrated exactly knot interval by knot interval, and `fourier_basis()`,
  whose derivatives of every order and whose antiderivative come from one
  phase-shift identity, and whose Gram matrix is diagonal in closed form over
  a whole period.

* Numerical fallbacks registered on the base class, so a basis that implements
  only `basis_eval()` is complete. Derivatives use one finite-difference
  stencil built from a Vandermonde solve, never a chain of first differences,
  and switch to a one-sided stencil at the interval endpoints.

* `basis_int()` is anchored: its value at the lower endpoint is exactly zero,
  for every family. Any antiderivative satisfies a differentiation check, so
  the constant is fixed by the contract instead.

* `check_basis()` runs six numerical checks and reports a quantity that came
  from a fallback as unchecked rather than as passed, and a property the
  family never claimed as not claimed.

* `basis_is_numerical()` answers the same question programmatically.

* `print()` and `plot()` methods; the plot draws the basis, any derivative, or
  the integral.
