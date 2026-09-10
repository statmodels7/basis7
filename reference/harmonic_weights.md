# The Weights of a Harmonic Acceleration Operator

Multiplies out \\r \prod\_{i \le h} (r^2 + (i\nu)^2)\\ and drops the
leading coefficient, which is 1.

## Usage

``` r
harmonic_weights(period, harmonics)
```

## Arguments

- period:

  The length of one cycle.

- harmonics:

  How many harmonics to leave in the null space.

## Value

A numeric vector of `2 * harmonics` entries, the weights \\w_0, \ldots,
w\_{m-1}\\.
