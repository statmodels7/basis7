# The Weight a Shrunk Null Space Carries

The penalty `null_space = "shrink"` puts on the directions the roughness
matrix does not see: one tenth of what a penalized direction carries.

## Usage

``` r
shrink_weight()
```

## Value

A single number.

## Details

It is mgcv's rule translated rather than a number chosen here. Measured
on `mgcv::s(bs = "ts")`, the shrinkage construction leaves the positive
eigenvalues of the penalty exactly as they were – 7708.76 down to 20.82,
identical to the unshrunk `"tp"` – and replaces each zero with 2.082,
which is a tenth of the smallest positive one. In Demmler-Reinsch
coordinates every penalized direction has eigenvalue exactly 1, so the
rule is the single number below.

What the tenth buys is measured. Both a weight of 0.1 and a weight of 1
let the term leave the model, the fitted values reaching a standard
deviation under 1e-6 at a large smoothing parameter. They differ in
rate: on a genuinely linear truth at a smoothing parameter of 100, the
root mean square error against that truth is 0.0723 at 0.1 and 0.4049 at
1, so the heavier weight destroys a real linear trend at a smoothing
parameter chosen to smooth the wiggles.
