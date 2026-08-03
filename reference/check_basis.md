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

  An object inheriting from class `basis`.

- n:

  The number of points at which to test.

- orders:

  The derivative orders to check.

- tol:

  The relative tolerance for the derivative and integral checks.

- verbose:

  Whether to print the outcome.

## Value

A named logical vector, invisibly, with `NA` for a check that was not
run. The attribute `"numerical"` records which generics fell back.

## Details

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

An order whose value comes from the numerical fallback is reported as
`[numerical]` rather than as passed. Checking such a value against a
numerical reference would be the same arithmetic twice, agreeing however
wrong the basis is, and a validator that reports agreement in that case
is worse than one that reports nothing: it says a thing was verified
when it was not.

The Gram check is run against a rule the basis does not itself use, so
that a basis whose own method is quadrature is still compared with
something independent.

The derivative and integral checks allow for the accuracy of their own
reference. A central difference assumes derivatives the function may not
have: at a knot a spline's third derivative jumps, and a stencil
straddling it returns the jump rather than the truncation error, which
would read as a failure of the basis. Each reference is therefore
computed twice, at a step and at half of it, and the gap between them
bounds its uncertainty; the comparison is allowed that much slack, point
by point. A deliberate error of five per cent is still caught by four
orders of magnitude, which is the check that the allowance has not
blunted anything.

## See also

[`basis_is_numerical`](https://statmodels7.github.io/basis7/reference/basis_is_numerical.md)

## Examples

``` r
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
```
