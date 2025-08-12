
# Direct effect: 0.14
# Indirect effect: 0.274
# Total effect: 0.399
library(data.table)
devtools::load_all()

dat <- readRDS(file = "data/mediation-data.rds")


# Phi11
phi11 <- mean(dat$Y[dat$A == 1])
# Phi00
phi00 <- mean(dat$Y[dat$A == 0])
# Phi10
phi10 <- mean(dat$Y10)

# Direct effect
direct_effect <- phi11 - phi10
direct_effect
# Indirect effect
indirect_effect <- phi10 - phi00
indirect_effect
# Total effect
total_effect <- phi11 - phi00
total_effect


dat <- as.data.table(dat)
df <- melt(dat,
           measure.vars = patterns("^L", "^M"),
           variable.name = "time",
           value.name = c("L", "M"))
df <- df[order(id, time), ]
df[,`:=`(
  lag_L= shift(L),
  lag_M= shift(M)
), by = id]
df$lag_L[is.na(df$lag_L)] <- 0
df$lag_M[is.na(df$lag_M)] <- 0

df$time <- as.numeric(df$time)

mod1 <- spec_model(
  L ~ X + A  + lag_L,
  mod_type = "covariate"
)
mod2 <- spec_model(
  M ~ X + A + L + lag_M,
  mod_type = "mediator"
)
mod3 <- spec_model(
  Y ~ X + A + L + M,
  var_type = "binomial",
  mod_type = "outcome"
)

fit.y <- glm(Y ~ X + A + L + M,
             data = df[time == 4,],
             family = binomial())
summary(fit.y)
broom::tidy(fit.y, exponentiate = T)

res.med <- mediation(
  data = df,
  id_var = 'id',
  base_vars = "X",
  exposure = "A",
  time_var = "time",
  models = list(mod1, mod2, mod3),
  init_recode = c("lag_L = 0",
                  "lag_M = 0"),
  in_recode = NULL,
  out_recode = c("lag_L = L",
                 "lag_M = M"),
  return_fitted = TRUE,
  mc_sample = 250000,
  return_data = FALSE,
  R = 5,
  quiet = TRUE,
  seed = 1
)

res.med
res.med$fitted_models

