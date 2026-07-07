# Derive parameters for a function from the current environment

Extracts the required parameters for a given function from the current
environment.

## Usage

``` r
get_args_for(fun, env = parent.frame(), ..., dots = NULL)
```

## Source

<https://stackoverflow.com/a/51002887>

## Arguments

- fun:

  The function whose parameters are to be extracted.

- env:

  The environment from which to extract parameters (default is the
  parent environment).

- ...:

  Additional arguments passed to the function.

- dots:

  Optional argument for handling extra dots.

## Value

A list of arguments for the specified function.
