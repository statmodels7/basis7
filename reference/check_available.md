# Check the Settings of a Smoother Against Each Other

Signals an error for a smoother that combines a `penalty` factory with
`null_space = "shrink"`, two settings that contradict each other.

## Usage

``` r
check_available(sm)
```

## Arguments

- sm:

  A
  [smoother](https://statmodels7.github.io/basis7/reference/smoother.md).

## Value

`NULL`, invisibly. The function is called for the error that it signals.
