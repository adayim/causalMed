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


## References
1. Tchetgen Tchetgen, E. J. (2013). Inverse odds ratio‐weighted estimation for causal mediation analysis. Statistics in medicine, 32(26), 4567-4580. [DOI:10.1002/sim.5864](https://doi.org/10.1002/sim.5864)
2. Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time‐varying exposures, mediators, and confounders. Statistics in medicine, 36(26), 4153-4166. [DOI:10.1002/sim.7426](https://doi.org/10.1002/sim.7426)
3. Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. Journal of causal inference, 5(2). [DOI:10.1515/jci-2016-0006](https://doi.org/10.1515/jci-2016-0006)

