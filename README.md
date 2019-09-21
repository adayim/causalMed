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
    ## Total effect              0.3136 -0.0177    0.4038  -0.4602    1.1228
    ## Natural Direct effect    -0.2574  0.0026    0.4639  -1.1692    0.6492
    ## Natural Indirect effect   0.5710 -0.0203    0.2561   0.0894    1.0932
    ## ------
    ## Proportion Mediated: 182.0809%

For time-varying mediator:

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
