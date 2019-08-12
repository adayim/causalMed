causalMed
=========

Causal Mediation analysis for time fixed and time-varying mediator. This
package is currently in development, please use with caution.

Installation
------------

    # Install the development version from GitHub:
    # install.packages("devtools")
    devtools::install_github("adayim/causalMed")

Usage
-----

For time fixed mediation analysis:

    library(causalMed)

    library(survival)
    data(lipdat)
    dtbase <- lipdat[lipdat$time == 0, ]   # Select the first row
    out <- iorw(coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase),
                exposure   = "smoke",
                mediator   = c("hdl", "ldl", "tg"),
                family     = "binomial")
    summary(out)

    ## Call:
    ## iorw(fitY = coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase), 
    ##     exposure = "smoke", mediator = c("hdl", "ldl", "tg"), family = "binomial")
    ## 
    ## Outcome Model Call:
    ## coxph(formula = Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase)
    ## 
    ## Exposure Model Call:
    ## glm(formula = smoke ~ bmi + age0 + smoke + hdl + ldl + tg, family = "binomial", 
    ##     data = dtbase)
    ## ------
    ## Natural effect model
    ## with standard errors based on the non-parametric bootstrap
    ## ---
    ## Exposure: smoke 
    ## Mediator(s): c, hdl, ldl, tg 
    ## ------
    ## Parameter estimates:
    ##                         Estimate    Bias Std.error conf.low conf.high
    ## Total effect              0.3136 -0.0183    0.4208  -0.4928    1.1566
    ## Natural Direct effect    -0.2574 -0.0011    0.5076  -1.2512    0.7386
    ## Natural Indirect effect   0.5710 -0.0172    0.2633   0.0720    1.1043
    ## ------
    ## Proportion Mediated: 182.0809%

For time-varying mediator:

    res <- medlong(data       = lipdat,
                   exposure   = "smoke",
                   mediator   = "hdl",
                   outcome    = "cvd",
                   id.var     = "id",
                   time.var   = "time",
                   covariates = c("bmi", "gender", "age0"),
                   m.family   = "gaussian",
                   y.family   = "binomial")
    summary(res)

    ## Call:
    ## medlong(data = lipdat, id.var = "id", exposure = "smoke", mediator = "hdl", 
    ##     outcome = "cvd", covariates = c("bmi", "gender", "age0"), 
    ##     time.var = "time", m.family = "gaussian", y.family = "binomial")
    ## ------
    ## Exposure: smoke 
    ## Mediator: hdlEstimation of standard errors based on the non-parametric bootstrap
    ## ---
    ## 
    ## ------
    ## Natural effect parameter estimates with G-formula estimation:
    ##                         Estimate    Bias Std.error conf.low conf.high
    ## Total effect              0.5828  0.0478    1.9370  -3.2614    4.3315
    ## Natural Indirect effect  -0.0033 -0.0433    0.8210  -1.5691    1.6490
    ## Natural Direct effect     0.5861  0.0910    2.1687  -3.7556    4.7457
    ## ------
    ## Proportion Mediated: -0.5701%
    ## ------
    ## 
    ## Natural effect parameter estimates with IPTW estimation:
    ##                         Estimate    Bias Std.error conf.low conf.high
    ## Total effect             -2.3198 -0.0575    0.4931  -3.2288   -1.2958
    ## Natural Indirect effect   0.2988  0.0199    0.2202  -0.1527    0.7104
    ## Natural Direct effect    -2.6185 -0.0774    0.4864  -3.4946   -1.5877
    ## ------
    ## Proportion Mediated: -12.8799%

Data structure must be in longitudinal format, and only mediator is
time-varying.

TODO
----

-   ☐ Time-varying treatment and covariates.
-   ☐ More documentations.
-   ☐ More tests is needed.

References
----------

1.  Tchetgen Tchetgen, E. J. (2013). Inverse odds ratio‐weighted
    estimation for causal mediation analysis. Statistics in medicine,
    32(26), 4567-4580.
    [DOI:10.1002/sim.5864](https://doi.org/10.1002/sim.5864)
2.  Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017).
    Mediation analysis for a survival outcome with time‐varying
    exposures, mediators, and confounders. Statistics in medicine,
    36(26), 4153-4166.
    [DOI:10.1002/sim.7426](https://doi.org/10.1002/sim.7426)
3.  Zheng, W., & van der Laan, M. (2017). Longitudinal mediation
    analysis with time-varying mediators and exposures, with application
    to survival outcomes. Journal of causal inference, 5(2).
    [DOI:10.1515/jci-2016-0006](https://doi.org/10.1515/jci-2016-0006)
