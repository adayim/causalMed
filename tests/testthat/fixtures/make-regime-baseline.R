# Generates tests/testthat/fixtures/regime-baseline.rds: the per-intervention
# estimates of mediation() under the DEFAULT regimes (always vs never
# exposed) for both mediation types, captured from the code as it stood
# before exposure_regime/reference_regime were added (2026-09-08).
# test-regressions.R asserts the defaults still reproduce these exactly.
# Regenerate ONLY when a deliberate engine change is meant to move them:
#   Rscript tests/testthat/fixtures/make-regime-baseline.R
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

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(baseline, "tests/testthat/fixtures/regime-baseline.rds")
str(baseline)
