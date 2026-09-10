# One Operation on a Penalty, However Many Components It Has

Applies `f` to a roughness matrix, or to each component of a list of
them, and returns the same shape it was given.

## Usage

``` r
over_penalty(s, f)
```

## Arguments

- s:

  A numeric matrix, or a list of them.

- f:

  A function of one matrix.

## Value

What `f` returns, or a list of what it returns for each component.

## Details

[`smoother_gram()`](https://statmodels7.github.io/basis7/reference/smoother_gram.md)
answers with a matrix for every family whose roughness is one quadratic
form, and with a list for one whose roughness is a sum of localized
components –
[`adaptive_smooth()`](https://statmodels7.github.io/basis7/reference/adaptive_smooth.md),
whose components carry a smoothing parameter each. Every step between
that answer and the built penalty is the same operation on each
component: the congruence of a reparametrization, the border a kept null
space adds, the removal of the dimnames. Writing the branch once here is
what keeps those steps from each growing one of their own.
