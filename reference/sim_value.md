# Random data simulation from predicted value.

Internal use only, predict response and simulate random data. The
simulated value will be restricted within the observed value range for
numeric values.

## Usage

``` r
sim_value(model, newdt)
```

## Arguments

- model:

  fitted objects defined in the \`spec_model\`.

- newdt:

  a data frame in which to look for variables with which to predict.

## Value

A simulated random vector using the predicted value from model and
newdt.
