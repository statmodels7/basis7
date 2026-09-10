# Check a Smoother's Penalty Factory

Validates the `penalty` argument of a smoother constructor: `NULL`, or a
function of one argument.

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

The check is deliberately weak, and the reason is the dependency graph.
basis7 sits at the bottom of it and imports numericals7 alone, so it
cannot name penalties7 and cannot ask whether what the function returns
is a penalty. It stores the function and never calls it. Whichever layer
builds the term calls it, at the coefficient count only the data settle,
and checks the result there.
