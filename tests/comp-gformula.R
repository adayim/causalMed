

library(gfoRmula)


id <- 'id'
time_points <- 7
time_name <- 't0'
covnames <- c('L1', 'L2', 'A')
outcome_name <- 'Y'
outcome_type <- 'survival'
covtypes <- c('binary', 'bounded normal', 'binary')
histories <- c(lagged, lagavg)
histvars <- list(c('A', 'L1', 'L2'), c('L1', 'L2'))
covparams <- list(covmodels = c(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
                                  L3 + t0,
                                L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
                                  lag_cumavg1_L2 + L3 + t0,
                                A ~ lag1_A + L1 + L2 + lag_cumavg1_L1 +
                                  lag_cumavg1_L2 + L3 + t0))
ymodel <- Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0
intvars <- list('A', 'A')
interventions <- list(list(c(static, rep(0, time_points))),
                      list(c(static, rep(1, time_points))))
int_descript <- c('Never treat', 'Always treat')
nsimul <- 10000

gform_basic <- gformula(obs_data = basicdata_nocomp,
                        id = id,
                        time_points = time_points,
                        time_name = time_name,
                        covnames = covnames,
                        outcome_name = outcome_name,
                        outcome_type = outcome_type,
                        covtypes = covtypes,
                        covparams = covparams,
                        ymodel = ymodel,
                        intvars = intvars,
                        interventions = interventions,
                        int_descript = int_descript,
                        histories = histories,
                        histvars = histvars,
                        basecovs = c('L3'),
                        nsimul = nsimul,
                        nsamples = 20,
                        seed = 1234)
gform_basic


devtools::load_all()

l1_mod <- spec_model(L1 ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
                       L3 + t0,
                       var_type = "binomial",
                       mod_type = "covriate"
)
l2_mod <- spec_model(L2 ~ lag1_A + L1 + lag_cumavg1_L1 +
                       lag_cumavg1_L2 + L3 + t0,
                     var_type = "normal",
                     mod_type = "covriate"
)
a_mod <- spec_model(A ~ lag1_A + lag_cumavg1_L1 + lag_cumavg1_L2 +
                       L3 + t0,
                     var_type = "binomial",
                     mod_type = "exposure"
)

y_mod <- spec_model(Y ~ A + L1 + L2 + L3 + lag1_A + lag1_L1 + lag1_L2 + t0,
                    var_type = "binomial",
                    mod_type = "survival"
)

dat <- copy(basicdata_nocomp)
dat[, c("lag1_A", "lag1_L1", "lag1_L2", "lag_cumavg1_L1", "lag_cumavg1_L2") :=  list(
  shift(A), shift(L1), shift(L2),
  cumsum(shift(L1, fill = 0)), cumsum(shift(L2, fill = 0))
), by = "id"]

init_recode <- c("lag1_A = 0",
                 "lag1_L1 = 0",
                 "lag1_L2 = 0",
                 "lag_cumavg1_L1 = 0",
                 "lag_cumavg1_L2 = 0")
out_recode <- c("lag1_A = A",
                "lag1_L1 = L1",
                "lag1_L2 = L2",
                "lag_cumavg1_L1 = L1 + lag_cumavg1_L1",
                "lag_cumavg1_L2 = L2 + lag_cumavg1_L2")

devtools::load_all()

debug(monte_g)

fit <- gformula(data = dat,
                id_var = "id",
                base_vars = c("A", "L3"),
                exposure = "A",
                time_var = "t0",
                models = list(l1_mod, l2_mod, a_mod, y_mod),
                intervention = list(never = 0, always = 1),
                ref_int = "natural",
                init_recode = init_recode,
                in_recode = NULL,
                out_recode = out_recode,
                return_fitted = FALSE,
                mc_sample = 10000,
                return_data = FALSE,
                R = 20,
                quiet = FALSE)

fit$estimate
fit$effect_size

gform_basic

