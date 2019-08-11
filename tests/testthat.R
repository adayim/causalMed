library(testthat)
library(causalMed)

test_check("causalMed")

data(lipdat)
dtbase <- lipdat[lipdat$time == 0, ]


out <- iorw(coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase),
            #data       = NULL,
            exposure   = "smoke",
            mediator   = "hdl",
            family     = "binomial",
            #ref        = NULL,
            stabilized = TRUE,
            R          = 1000)

res <- medlong(data = dat,
               trt = "a",
               med = "mt",
               y = "status",
               id = "id",
               time = "tij",
               cov = c("w1", "w2"),
               m.family  = "gaussian",
               y.family  = "binomial")



