options(stringsAsFactors = FALSE)

df <- read.csv("data-raw/gvhd_data.csv")

devtools::load_all()

set.seed(12131231)


mod_cov1 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
                         agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
                    family = "binomial",
                    type = "covriate",
                    subset = platnormm1 == 0)

mod_cov2 <- spec_model(relapse ~ all + cmv + male + age + gvhdm1 + daysgvhd + platnormm1 +
                         daysnoplatnorm + agecurs1 + agecurs2 + day + daysq + wait,
                       family = "binomial",
                       type = "covriate",
                       subset = relapsem1 == 0)

mod_exp <- spec_model(gvhd ~ all + cmv + male + age + platnormm1 +
                        daysnoplatnorm + relapsem1 + daysnorelapse +
                        agecurs1 + agecurs2 + day + daysq + wait,
                      subset = gvhdm1 == 0,
                      family = "binomial",
                      type = "exposure")

mod_cens <- spec_model(censlost ~ all + cmv + male + age + daysgvhd +
                         daysnoplatnorm + daysnorelapse + agesq + day +
                         daycurs1 + daycurs2 + wait,
                       family = "binomial",
                       type = "censor")

mod_out <- spec_model(d ~ all + cmv + male + age + gvhd + platnorm +
                        daysnoplatnorm + relapse + daysnorelapse +
                        agesq + day + daycurs1 + daycurs2 + wait +
                        day*gvhd + daycurs1*gvhd + daycurs2*gvhd,
                      family = "binomial",
                      type = "outcome")

init_recode <- c("relapse=0", "gvhd=0", "platnorm=0", "gvhdm1=0",
                "relapsem1=0", "platnormm1=0", "daysnorelapse=0",
                "daysnoplatnorm=0", "daysnogvhd=0", "daysrelapse=0",
                "daysplatnorm=0", "daysgvhd=0")
in_recode = c('daysq = day^2', 'daycu = day^3',
              "daycurs1 = ((day>83.6)*((day-83.6)/83.6)^3)+((day>1862.2)*((day-1862.2)/83.6)^3)*(947.0-83.6) -((day>947.0)*((day-947.0)/83.6)^3)*(1862.2-83.6)/(1862.2-947.0)",
              "daycurs2 = ((day>401.4)*((day-401.4)/83.6)^3)+((day>1862.2)*((day-1862.2)/83.6)^3)*(947.0-401.4) -((day>947.0)*((day-947.0)/83.6)^3)*(1862.2-401.4)/(1862.2-947.0)")
out_recode = c('platnormm1 = platnorm',
               'relapsem1 = relapse',
               'gvhdm1 = gvhd',
               'daysnorelapse = ifelse(relapse == 0, daysnorelapse + 1,
                                                            daysnorelapse)',
               'daysrelapse = ifelse(relapse == 1, daysrelapse,
                                                          daysrelapse + 1)',

               'daysnoplatnorm = ifelse(platnorm == 0, daysnoplatnorm + 1,
                                                          daysnoplatnorm)',
               'daysplatnorm = ifelse(relapse == 1, daysplatnorm,
                                                          daysplatnorm + 1)',

               'daysnogvhd = ifelse(gvhd == 0, daysnogvhd + 1,
                                                          daysnogvhd)',
               'daysgvhd = ifelse(relapse == 1, daysgvhd,
                                                          daysgvhd + 1)',

               'platnorm = ifelse(platnormm1 == 1, 1, platnorm)',
               'relapse = ifelse(relapsem1 == 1, 1, relapse)',
               'gvhd = ifelse(gvhdm1 == 1, 1, gvhd)',
               'daysplatnorm = ifelse(platnorm == 0, daysplatnorm,
                                                        daysplatnorm + 1)')


debug(Gformula)
# debug(monte_g)
# debug(monte_sim)
natural <- Gformula(df,
                id.var = "id",
                base.vars = c("age", "agesq", "agecurs1", "agecurs2", "male",
                              "cmv", "all", "wait"),
                exposure = "gvhd",
                outcome = "d",
                time.var = "day",
                models = list(mod_cov1, mod_cov2, mod_exp, mod_cens, mod_out),
                init.recode = init_recode,
                in.recode = in_recode,
                out.recode = out_recode,
                mc.sample = 13700)

always <- Gformula(df,
                id.var = "id",
                base.vars = c("age", "agesq", "agecurs1", "agecurs2", "male",
                              "cmv", "all", "wait"),
                exposure = "gvhd",
                outcome = "d",
                time.var = "day",
                intervention = 1,
                models = list(mod_exp, mod_cov1, mod_cov2, mod_cens, mod_out),
                init.recode = init_recode,
                in.recode = in_recode,
                out.recode = out_recode,
                mc.sample = 13700)

never <- Gformula(df,
                id.var = "id",
                base.vars = c("age", "agesq", "agecurs1", "agecurs2", "male",
                              "cmv", "all", "wait"),
                exposure = "gvhd",
                outcome = "d",
                time.var = "day",
                intervention = 0,
                models = list(mod_exp, mod_cov1, mod_cov2, mod_cens, mod_out),
                init.recode = init_recode,
                in.recode = in_recode,
                out.recode = out_recode,
                mc.sample = 13700)

library(survminer)
library(survival)
library(tidyverse)
base_df <- df %>%
  group_by(id) %>%
  summarise(day = max(day),
            d = max(d)) %>%
  mutate(group = "Observed")

out_dat <- natural$out %>%
  select(id, day, d) %>%
  mutate(group = "Natural")

res_dat <- rbind.data.frame(base_df, out_dat)

fit <- survfit(Surv(day, d) ~ group,
               data = res_dat)

ggsurvplot(fit)


# Compare
out_dat <- always$out %>%
  select(id, day, d) %>%
  mutate(group = "Always") %>%
  bind_rows(., never$out %>%
              select(id, day, d) %>%
              mutate(group = "Never"))

res_cox <- coxph(Surv(day, d) ~ gvhd, data = out_dat, ties = "efron")
summary(res_cox)
fit2 <- survfit(Surv(day, d) ~ gvhd,
               data = out_dat)

ggsurvplot(fit2)



