# causalMed
Causal Mediation analysis for time fixed and time-varying mediator. This package is currently in development, please use with caution. 

## Installation

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("adayim/causalMed")
```

## Usage
For time fixed mediation analysis:
```r
dat <- readRDS("data/data.rds")
dt <- dat[dat$tij == 0, ]    # Select the first row
res <- iorw(data = dt, trt = "a", med = "mt",
            y = "status", time = "eventtime",
            family = "cox",
            cov = c("w1", "w2"))
```

For time-varying mediator:
```r
dat <- readRDS("data/data.rds")
res <- medlong(data = dat,
               trt = "a",
               med = "mt",
               y = "status",
               id = "id",
               time = "tij",
               cov = c("w1", "w2"),
               m.family  = "gaussian",
               y.family  = "binomial")
```
Data structure must be in longitudinal format, and only mediator is time-varying.

## TODO
- [ ] Time-varying treatment and covariates.
- [ ] More documentations.
- [ ] More tests is needed.


