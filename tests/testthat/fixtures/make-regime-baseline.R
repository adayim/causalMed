# Generates tests/testthat/fixtures/regime-baseline.rds: the per-intervention
# estimates of mediation() under the DEFAULT regimes (always vs never
# exposed) for both mediation types.
out_file <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_file) || !nzchar(out_file))
  stop("Supply the output .rds path as the first argument.", call. = FALSE)

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
data("nonsurvivaldata", package = "causalMed")

m_L2  <- spec_model(L2    ~ A + V + time,      var_type = "binary", mod_type = "covariate")
m_med <- spec_model(L1    ~ A + V + L2 + time, var_type = "normal", mod_type = "mediator")
m_Y   <- spec_model(Y_bin ~ A + L1 + L2 + V,   var_type = "binary", mod_type = "outcome")
models <- list(m_L2, m_med, m_Y)

baseline <- lapply(c(I = "I", N = "N"), function(mtype) {
  fit <- suppressWarnings(mediation(
    data = nonsurvivaldata, id_var = "id", base_vars = "V", exposure = "A",
    outcome = "Y_bin", time_var = "time", models = models,
    mediation_type = mtype, mc_sample = 500L, R = 1L, quiet = TRUE, seed = 42L))
  list(effect_size = setNames(fit$effect_size$Est, fit$effect_size$Intervention),
       estimate    = setNames(fit$estimate$RD, fit$estimate$Effect))
})

saveRDS(baseline, out_file)
str(baseline)
