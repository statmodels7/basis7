# What a Smoother Removes and What It Gives Back

Returns the directions a smoother constrains its block against, and the
ones `null_space = "keep"` restores as free columns. Both are evaluated
at the covariate, one column per direction.

## Usage

``` r
smoother_span(sm, x, ...)

smoother_span_apply(sm, params, newx, ...)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

- x:

  The covariate, a numeric vector.

- ...:

  Passed to methods.

- params:

  For `smoother_span_apply()`, the `params` element of the
  `smoother_span()` result.

- newx:

  For `smoother_span_apply()`, the new covariate values.

## Value

`smoother_span()` returns a list of three elements: `constraint`, a
numeric matrix of one column per direction removed, or `NULL` for none;
`free`, a numeric matrix of one column per direction restored, with zero
columns where none is; and `params`, what `smoother_span_apply()` needs
to rebuild `free` at new values. `smoother_span_apply()` returns that
matrix at `newx`.

## Why the pair and not the basis

The null space is a property of the basis and the penalty **together**,
not of the basis alone. The second-derivative Gram matrix of a cubic
B-spline has the polynomials of degree below 2 in its null space and the
order-3 matrix the polynomials of degree below 3, while a Fourier basis
has the constant alone at every order, the basis containing no linear
function. Measured, the null function of a Fourier Gram matrix has a
standard deviation of exactly zero at orders 1, 2 and 3. So the answer
belongs to the smoother, which carries both, and a family declares it by
registering a method here.

## The default, which serves every polynomial family

The base method removes the polynomials of degree below `order`, or up
to `constrain` when that is given and larger, and restores all but the
constant. The constant is not restored because a model carrying an
intercept already spans it, which is the convention
[`dr_basis()`](https://statmodels7.github.io/basis7/reference/dr_basis.md)
has always followed.

The restored columns are the raw powers made orthogonal to one another
and to the constant over the observed covariate, then standardized. The
first is `(x - mean(x)) / sd(x)`, so at `order = 2` there is exactly one
and it is the standardized covariate.

Directions of the constraint **beyond** the null space are removed and
not restored. Returning them as free columns would put a parametric term
inside the smooth instead of in the formula, where a reader sees it.

## See also

[`smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.md),
which calls both.

## Examples

``` r
set.seed(1)
x <- sort(runif(60))
sp <- smoother_span(bspline_smooth(k = 8), x)
c(removed = ncol(sp$constraint), restored = ncol(sp$free))
#>  removed restored 
#>        2        1 

# the restored column is the standardized covariate
round(c(mean = mean(sp$free), sd = sd(sp$free), cor = cor(sp$free, x)), 12)
#> mean   sd  cor 
#>    0    1    1 

# reapplied at new values rather than rebuilt
i <- c(2L, 20L, 55L)
max(abs(smoother_span_apply(bspline_smooth(k = 8), sp$params, x[i]) -
        sp$free[i, , drop = FALSE]))
#> [1] 0
```
