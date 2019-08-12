library(testthat)
library(causalMed)

test_check("causalMed")

data(lipdat)
dtbase <- lipdat[lipdat$time == 0, ]


out <- iorw(coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase),
            #data       = lipdat,
            exposure   = "smoke",
            mediator   = "hdl",
            family     = "binomial",
            stabilized = TRUE,
            R          = 1000)

res <- medlong(data       = lipdat,
               exposure   = "smoke",
               mediator   = "hdl",
               outcome    = "cvd",
               id.var     = "id",
               time.var   = "time",
               covariates = c("bmi", "gender", "age0"),
               m.family   = "gaussian",
               y.family   = "binomial")

