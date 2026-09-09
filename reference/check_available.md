# Check What a Smoother Asks For Against What Is Built

Signals an error for a smoother whose settings this version of the
package does not build. The arguments are on the constructor because
they are part of its interface; the values below reach arithmetic that
is not written yet, and an argument accepted and ignored would report a
fit of a model the caller did not ask for.

## Usage

``` r
check_available(sm)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## Value

`NULL`, invisibly. Called for the error it signals.
