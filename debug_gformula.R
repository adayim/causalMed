
library(devtools)
load_all()
data(nonsurvivaldata)
  
mod1 <- spec_model(A ~ V, var_type = "binomial", mod_type = "exposure")
mod2 <- spec_model(L1 ~ A + V, var_type = "normal", mod_type = "covariate")
mod3 <- spec_model(Y_bin ~ A + L1 + V, var_type = "binomial", mod_type = "outcome")
models <- list(mod1, mod2, mod3)

ints <- list(natural = NULL, always = 1)

fit <- gformula(
  data = nonsurvivaldata,
  id_var = "id",
  base_vars = "V",
  exposure = "A",
  time_var = "time",
  models = models,
  intervention = ints,
  ref_int = 0, # natural
  mc_sample = 500, 
  R = 0, 
  quiet = TRUE
)

print("Estimate table:")
print(fit$estimate)
print("Unique Interventions in estimate:")
print(unique(fit$estimate$Intervention))
