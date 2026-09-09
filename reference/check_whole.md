# Check a Whole-Number Argument

Validates one of the whole-number arguments of a smoother constructor
and returns it as an integer.

## Usage

``` r
check_whole(v, nm, lo, why = "")
```

## Arguments

- v:

  The value.

- nm:

  The argument's name, for the message.

- lo:

  The smallest value admitted.

- why:

  An explanation appended to the message, or `""`.

## Value

`v` as a length-one integer.
