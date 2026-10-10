# The Weight of a Shrunk Null Space

The penalty that `null_space = "shrink"` puts on the directions that the
roughness matrix does not see: one tenth of the weight of a penalized
direction.

## Usage

``` r
shrink_weight()
```

## Value

A single number.

## Details

The weight follows the shrinkage construction of the `bs = "ts"` smooths
of mgcv, which keeps the positive eigenvalues of the penalty and
replaces each zero eigenvalue with a tenth of the smallest positive one.
In Demmler-Reinsch coordinates every penalized direction has eigenvalue
1, so the rule is the single number returned here.

With this weight the fitted values of the term tend to a constant as the
smoothing parameter grows, so the term can leave the model, while a
linear component is shrunk more slowly than the oscillating ones.
