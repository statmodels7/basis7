# The Weights of an Oscillator Operator

Multiplies out \\\prod\_{i \le h} (r^2 + (i\nu)^2)\\ and drops the
leading coefficient, which is 1.

## Usage

``` r
oscillator_weights(period, harmonics)
```

## Arguments

- period:

  The length of one cycle.

- harmonics:

  How many harmonics to leave in the null space.

## Value

A numeric vector of `2 * harmonics` entries.
