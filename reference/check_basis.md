# Validate a Basis

Runs the numerical checks that a basis must pass, and prints the outcome
of each. It is meant above all for a basis written outside the package,
where a hand-derived derivative or a misplaced constant of integration
is the likeliest mistake.

## Usage

``` r
check_basis(basis, n = 41L, orders = 1:2, tol = 1e-06, verbose = TRUE)
```

## Arguments

- basis:

  A basis object, of any class inheriting from
  [basis](https://statmodels7.github.io/basis7/reference/basis.md).

- n:

  The number of points to test at, default `41`. They are equally spaced
  across the interval with five per cent trimmed off each end, because a
  one-sided stencil at an endpoint carries a larger error that would
  read as a failure of the basis. Passed through
  [`as.integer()`](https://rdrr.io/r/base/integer.html).

- orders:

  The derivative orders to check, default `1:2`. Each is compared
  against one differentiation of the order below it, never against a
  chain of lower-order differences. Raising it tests more and costs one
  pass per order.

- tol:

  The relative tolerance for the derivative and integral checks, default
  `1e-6`. Relative to the values themselves, with the denominator taken
  from them; see
  [`rel_close()`](https://statmodels7.github.io/basis7/reference/rel_close.md).
  The Gram comparison uses `1e-6` regardless.

- verbose:

  Whether to print the table, default `TRUE`. The result is returned
  invisibly either way.

## Value

Invisibly, a named logical vector of length 6, in the order `shape`,
`deriv`, `integral`, `partition`, `gram`, `missing`, with `NA` for a
check that was not run. It carries the attribute `"numerical"`, the
named logical vector
[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
returns.

## Details

The identities verified, for a basis \\\varphi_1, \dots, \varphi_d\\ on
\\\[a, b\]\\, are

\$\$\varphi_j^{(k)}(x) = \frac{\mathrm{d}}{\mathrm{d}x}\\
\varphi_j^{(k-1)}(x), \qquad \frac{\mathrm{d}}{\mathrm{d}x} \int_a^x
\varphi_j(t)\\\mathrm{d}t = \varphi_j(x), \qquad \int_a^a \varphi_j =
0,\$\$

together with \\\sum_j \varphi_j(x) = 1\\ where the family has that
property and \\G\_{jl} = \int \varphi_j \varphi_l \\\mathrm{d}\mu =
G\_{lj}\\ with \\G \succeq 0\\ for the Gram matrix.

The checks are:

1.  **shape**: every generic returns a matrix of the declared size, with
    the declared column names, for a vector and for a single point;

2.  **derivatives**: each analytic order agrees with one numerical
    differentiation of the order below it;

3.  **integral**: it differentiates back to the basis, and is exactly
    zero at the lower endpoint;

4.  **partition of unity**: the rows sum to one, for the families that
    have that property;

5.  **Gram**: symmetric, positive semidefinite, and equal to an
    independent quadrature;

6.  **missing values**: a missing evaluation point gives a missing row
    and nothing else.

## What is skipped, and what is weakened

Where a quantity comes from the numerical fallback, checking it against
a numerical reference would be the same arithmetic twice, agreeing
however wrong the basis is. The `deriv` and `integral` checks are
therefore not run at all in that case: they report `NA` and print
`[numerical]`.

The `gram` check is **weakened, not skipped**. Symmetry and positive
semidefiniteness are properties of the matrix whatever produced it, so
they are still tested and the row still prints `[PASSED]`; what is
dropped is the comparison against an independent quadrature. The `shape`
and `missing` checks read nothing about the fallback and always run.

The independent quadrature is 401 panels of 7 nodes, a rule no basis in
the package uses for its own answer, so a family whose own Gram matrix
is a quadrature is still compared against different arithmetic.

`partition` is `NA` and prints `[not claimed]` for a family that is not
a partition of unity. That is a different thing from `[numerical]`, and
the printout distinguishes them: a property the family never claimed was
not applicable, where a numerical value was simply not verified.

## How the reference's own error is allowed for

A central difference assumes derivatives the function may not have. At a
knot a spline's third derivative jumps, and a stencil straddling it
returns a number of the order of the jump, which compared against an
exact value reads as a failure of the basis. Each reference is therefore
computed twice, at a step and at half of it, and the gap between them
bounds its own uncertainty; the comparison is allowed that much slack
point by point, so each point contributes exactly the accuracy its
reference supports. See
[`fd_reference()`](https://statmodels7.github.io/basis7/reference/fd_reference.md).

A deliberate error of five per cent is still caught by four orders of
magnitude, which is the check that the allowance has not blunted
anything.

## See also

[`basis_is_numerical()`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)
for the attribute and which route each quantity takes;
[`vignette("defining-a-basis")`](https://statmodels7.github.io/basis7/articles/defining-a-basis.md),
which uses this function to develop one.

## Examples

``` r
# Every shipped family passes all six checks.
invisible(check_basis(bspline_basis(dimension = 6)))
#> check_basis: bspline (6 functions)
#>   shape       shapes and column names                        [PASSED]
#>   deriv       derivatives against finite differences         [PASSED]
#>   integral    integral differentiates back, zero at lower    [PASSED]
#>   partition   partition of unity                             [PASSED]
#>   gram        Gram symmetric, PSD, matches quadrature        [PASSED]
#>   missing     missing values give missing rows               [PASSED]
invisible(check_basis(fourier_basis(dimension = 5)))
#> check_basis: fourier (5 functions)
#>   shape       shapes and column names                        [PASSED]
#>   deriv       derivatives against finite differences         [PASSED]
#>   integral    integral differentiates back, zero at lower    [PASSED]
#>   partition   partition of unity                             [not claimed]
#>   gram        Gram symmetric, PSD, matches quadrature        [PASSED]
#>   missing     missing values give missing rows               [PASSED]

# The result is a logical vector, and the attribute says which quantities
# were computed from a formula.
r <- check_basis(poly_basis(dimension = 5), verbose = FALSE)
r
#>     shape     deriv  integral partition      gram   missing 
#>      TRUE      TRUE      TRUE        NA      TRUE      TRUE 
#> attr(,"numerical")
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 
attr(r, "numerical")
#> basis_deriv   basis_int  basis_gram 
#>       FALSE       FALSE       FALSE 

# A basis defined from its evaluation alone: the two checks that would
# compare a difference against a difference are not run, and the Gram check
# keeps its symmetry and definiteness half.
Bumps <- S7::new_class("Bumps", parent = basis)
S7::method(basis_eval, Bumps) <- function(basis, x, ...) {
  out <- exp(-0.5 * outer(x, seq(0, 1, length.out = basis@dimension),
                          "-")^2 / 0.12^2)
  colnames(out) <- basis_colnames(basis)
  out
}
#> Overwriting method basis_eval(<Bumps>)
invisible(check_basis(Bumps(basis_name = "bumps", dimension = 4L,
                            lower = 0, upper = 1)))
#> check_basis: bumps (4 functions)
#>   shape       shapes and column names                        [PASSED]
#>   deriv       derivatives against finite differences         [numerical]
#>   integral    integral differentiates back, zero at lower    [numerical]
#>   partition   partition of unity                             [not claimed]
#>   gram        Gram symmetric, PSD, matches quadrature        [PASSED]
#>   missing     missing values give missing rows               [PASSED]
#>   computed numerically: basis_deriv, basis_int, basis_gram

# A derivative five per cent wrong is caught.
Wrong <- S7::new_class("Wrong", parent = BsplineBasis)
S7::method(basis_deriv, Wrong) <- function(basis, x, order = 1L, ...) {
  1.05 * S7::method(basis_deriv, BsplineBasis)(basis, x, order = order)
}
b <- bspline_basis(dimension = 6)
invisible(check_basis(Wrong(basis_name = "wrong", dimension = b@dimension,
                            lower = 0, upper = 1,
                            basis_params = b@basis_params)))
#> check_basis: wrong (6 functions)
#>   shape       shapes and column names                        [PASSED]
#>   deriv       derivatives against finite differences         [FAILED]
#>   integral    integral differentiates back, zero at lower    [PASSED]
#>   partition   partition of unity                             [PASSED]
#>   gram        Gram symmetric, PSD, matches quadrature        [PASSED]
#>   missing     missing values give missing rows               [PASSED]
```
