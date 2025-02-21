
devtools::load_all()
data(gvhd)

gvhd <- within(gvhd, c(
  agesq <- age^2,
  agecurs1 <- (age > 17.0) * (age - 17.0)^3 - ((age > 30.0) * (age - 30.0)^3) * (41.4 - 17.0) / (41.4 - 30.0),
  agecurs2 <- (age > 25.4) * (age - 25.4)^3 - ((age > 41.4) * (age - 41.4)^3) * (41.4 - 25.4) / (41.4 - 30.0),
  daysq <- day^2,
  daycu <- day^3,
  daycurs1 <- ((day > 63) * ((day - 63) / 63)^3) + ((day > 716) * ((day - 716) / 63)^3) * (350.0 - 63) - ((day > 350) * ((day - 350) / 63)^3) * (716 - 63) / (716 - 350),
  daycurs2 <- ((day > 168) * ((day - 168) / 63)^3) + ((day > 716) * ((day - 716) / 63)^3) * (350 - 168) - ((day > 350) * ((day - 350) / 63)^3) * (716 - 168) / (716 - 350)
))

# Parametric G-formula coefficient estimation models
mod_cov1 <- spec_model(relapse ~ all + cmv + male + age + gvhdm1 + daysgvhd + platnormm1 +
  daysnoplatnorm + agecurs1 + agecurs2 + day + daysq + wait,
var_type = "binomial",
mod_type = "covriate",
subset = relapsem1 == 0
)

# Model for probability of platnorm=1 at day k
mod_cov2 <- spec_model(platnorm ~ all + cmv + male + age + agecurs1 +
  agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait,
var_type = "binomial",
mod_type = "covriate",
subset = platnormm1 == 0
)

# Model for probability of exposure=1 at day k
mod_exp <- spec_model(gvhd ~ all + cmv + male + age + platnormm1 +
  daysnoplatnorm + relapsem1 + daysnorelapse +
  agecurs1 + agecurs2 + day + daysq + wait,
subset = gvhdm1 == 0,
var_type = "binomial",
mod_type = "exposure"
)

# Model for probability of censoring=1 at day k
mod_cens <- spec_model(censlost ~ all + cmv + male + age + daysgvhd +
  daysnoplatnorm + daysnorelapse + agesq + day +
  daysq + daycu + wait,
var_type = "binomial",
mod_type = "censor"
)

# Model for probability of outcome=1 at day k
mod_out <- spec_model(d ~ all + cmv + male + age + gvhd + platnorm +
  daysnoplatnorm + relapse + daysnorelapse +
  agesq + wait + day * gvhd + daysq * gvhd +
  daycu * gvhd,
var_type = "binomial",
mod_type = "survival"
)

init_recode <- c(
  "relapse=0", "gvhd=0", "platnorm=0", "gvhdm1=0",
  "relapsem1=0", "platnormm1=0", "daysnorelapse=0",
  "daysnoplatnorm=0", "daysnogvhd=0", "daysrelapse=0",
  "daysplatnorm=0", "daysgvhd=0"
)

in_recode <- c(
  "daysq = day^2", "daycu = day^3",
  "platnormm1 = platnorm",
  "relapsem1 = relapse",
  "gvhdm1 = gvhd",
  "daycurs1 = ((day>83.6)*((day-83.6)/83.6)^3)+((day>1862.2)*((day-1862.2)/83.6)^3)*(947.0-83.6) -((day>947.0)*((day-947.0)/83.6)^3)*(1862.2-83.6)/(1862.2-947.0)",
  "daycurs2 = ((day>401.4)*((day-401.4)/83.6)^3)+((day>1862.2)*((day-1862.2)/83.6)^3)*(947.0-401.4) -((day>947.0)*((day-947.0)/83.6)^3)*(1862.2-401.4)/(1862.2-947.0)"
)

out_recode <- c(
  "daysnorelapse = ifelse(relapse == 0, daysnorelapse + 1, daysnorelapse)",
  "daysrelapse = ifelse(relapse == 1, daysrelapse, daysrelapse + 1)",
  "daysnoplatnorm = ifelse(platnorm == 0, daysnoplatnorm + 1, daysnoplatnorm)",
  "daysplatnorm = ifelse(relapse == 1, daysplatnorm, daysplatnorm + 1)",
  "daysnogvhd = ifelse(gvhd == 0, daysnogvhd + 1, daysnogvhd)",
  "daysgvhd = ifelse(relapse == 1, daysgvhd, daysgvhd + 1)",
  "platnorm = ifelse(platnormm1 == 1, 1, platnorm)",
  "relapse = ifelse(relapsem1 == 1, 1, relapse)",
  "gvhd = ifelse(gvhdm1 == 1, 1, gvhd)",
  "daysplatnorm = ifelse(platnorm == 0, daysplatnorm, daysplatnorm + 1)"
)

# devtools::load_all()
mods <- gformula(gvhd,
  id_var = "id",
  base_vars = c(
    "age", "agesq", "agecurs1", "agecurs2", "male",
    "cmv", "all", "wait"
  ),
  exposure = "gvhd",
  time_var = "day",
  models = list(mod_cov1, mod_cov2, mod_exp, mod_cens, mod_out),
  intervention = list(never = 0),
  init_recode = init_recode,
  in_recode = in_recode,
  out_recode = out_recode,
  mc_sample = 137000,
  R = 0
)

# save(mods, file = "data-raw/test.RData")
mods <- mods$gform.data
gc()
setkeyv(mods, c("Intervention", "new_ID", "day"))

dt <- mods[, .SD[.N], by = c("Intervention", "new_ID")]
dt$Intervention <- factor(dt$Intervention, levels = c("never", "natural"))

fit <- coxph(Surv(day, d) ~ Intervention, data = dt)
summary(fit)



