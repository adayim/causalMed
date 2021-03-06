library(testthat)
library(causalMed)

set.seed(20190904)

test_check("causalMed")

# For time-fixed data
data(lipdat)
dtbase <- lipdat[lipdat$time == 0, ]
out <- iorw(coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase),
            #data       = lipdat,
            exposure   = "smoke",
            mediator   = "hdl",
            family     = "binomial",
            stabilized = TRUE,
            R          = 1000)

# For longitudinal data
df <- read.csv("data-raw/gvhd_data.csv")

mod_exp <- spec_model(gvhd ~ all + cmv + male + age + agecurs1 + agecurs2 +
                        platnormm1 + daysnoplatnorm + relapsem1 + daysnorelapse
                      + day + daysq + wait,
                      subset = gvhdm1 == 0,
                      family = "binomial",
                      type = "exposure",
                      recode = c('daysq = day^2', 'daycu = day^3'),
                      order = 1)

mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
                         agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
                       family = "binomial",
                       type = "covriate",
                       order = 2,
                       recode = c('platnorm = ifelse(platnormm1 == 1, 1, platnorm)',
                                  'daysnoplatnorm = ifelse(platnorm == 0, daysnoplatnorm + 1,
                                                          daysnoplatnorm)',
                                  'daysplatnorm = ifelse(platnorm == 0, daysplatnorm,
                                                        daysplatnorm + 1)'),
                       subset = platnormm1 == 0)

mod_cov2 <- spec_model(relapse ~ all + cmv + male + age + agecurs1 + agecurs2 +
                         gvhdm1 +  daysgvhd + platnormm1 + daysnoplatnorm + day +
                         daysq + wait,
                       family = "binomial",
                       type = "covriate",
                       order = 3,
                       recode = c('relapse = ifelse(relapsem1 == 1, 1, relapse)',
                                  'daysnorelapse = ifelse(relapse == 0, daysnorelapse + 1,
                                                            daysnorelapse)',
                                  'daysrelapse = ifelse(relapse == 0, daysrelapse,
                                                          daysrelapse + 1)'),
                       subset = relapsem1 == 0)

mod_cens <- spec_model(censlost ~ all + cmv + male + age + agesq + daysgvhd +
                         daysnoplatnorm + daysnorelapse + day + daysq +
                         daycu + wait,
                       family = "binomial",
                       type = "censor",
                       order = 4)

mod_out <- spec_model(d ~ gvhd + cmv + male + age  + agesq + day + daysq +
                        daycu + platnorm + daysnoplatnorm + relapse +
                        daysnorelapse + wait + all + day:gvhd +
                        daysq:gvhd + daycu:gvhd,
                      family = "binomial",
                      type = "outcome",
                      order = 5)

res <- Gformula(df,
                id.var = "id",
                base.vars = c("age", "agesq", "agecurs1", "agecurs2", "male",
                              "cmv", "all", "wait"),
                exposure = "gvhd",
                outcome = "d",
                time.var = "day",
                models = list(mod_exp, mod_cov1, mod_cov2, mod_cens, mod_out),
                init.recode = c("daysq = day^2", "daycu = day^3", "relapse=0", "gvhd=0",
                                "platnorm=0", "gvhdm1=0", "relapsem1=0", "platnormm1=0",
                                "daysnorelapse=0", "daysnoplatnorm=0", "daysnogvhd=0",
                                "daysrelapse=0", "daysplatnorm=0", "daysgvhd=0"),
                out.recode = c('platnormm1 = platnorm', 'relapsem1 = relapse', 'gvhdm1 = gvhd'),
                mc.sample = 13700)


