# Check a Smoother's Penalty Factory

Validates the `penalty` argument of a smoother constructor: `NULL`, or a
function with at least one argument.

## Usage

``` r
check_penalty(penalty)
```

## Arguments

- penalty:

  The value given.

## Value

`penalty`, unchanged.

## Details

The check is deliberately weak because of the dependency graph: of the
toolkit basis7 imports numericals7 only, so it cannot name penalties7
and cannot test whether the function returns a penalty. It stores the
function and never calls it. The layer that builds the term calls it, at
the coefficient count that only the data settle, and checks the result
there.
