
test_that("gformula runs successfully on nonsurvivaldata", {
  data(nonsurvivaldata)


  # Basic models
  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod2 <- spec_model(L1 ~ A + V, var_type = "normal", mod_type = "covariate")
  mod3 <- spec_model(Y_bin ~ A + L1 + V, var_type = "binary", mod_type = "outcome")
  models <- list(mod1, mod2, mod3)

  ints <- list(natural = NULL, always = 1)

  # Basic Run
  fit <- gformula(
    data = nonsurvivaldata,
    id_var = "id",
    base_vars = "V",
    exposure = "A",
    time_var = "time",
    models = models,
    intervention = ints,
    ref_int = 0, # natural
    mc_sample = 500, # small sample for speed
    R = 0, # no bootstrap for speed
    quiet = TRUE
  )

  expect_s3_class(fit, "gformula")
  expect_true(is.list(fit))
  expect_named(fit, c("call", "all.args", "estimate", "effect_size", "sim_data",
                      "fitted_models", "boot_estimates", "data_summary", "observed"))

  # No bootstrap (R = 0): per-replicate estimates absent, summaries present
  expect_null(fit$boot_estimates)
  expect_equal(fit$data_summary$n_id, length(unique(nonsurvivaldata$id)))
  expect_equal(fit$data_summary$n_obs, nrow(nonsurvivaldata))
  expect_true(is.finite(fit$observed$value))

  # Check effect_size
  expect_true(inherits(fit$effect_size, "data.table"))
  expect_true(all(c("Intervention", "Est") %in% names(fit$effect_size)))
  expect_equal(nrow(fit$effect_size), 2) # natural + always

  # Check estimate (contrast)
  expect_true(inherits(fit$estimate, "data.table"))
  expect_true(all(c("Intervention", "Risk_type", "Estimate") %in% names(fit$estimate)))
  # Should have Risk Difference and Risk Ratio for "always" vs "natural"
  expect_true(all(grepl("always", unique(fit$estimate$Intervention))))

  # Check result correctness roughly (not precise, just logic)
  expect_true(all(fit$effect_size$Est >= 0 & fit$effect_size$Est <= 1))

})

test_that("gformula works with return_data = TRUE and return_fitted = TRUE", {
  data(nonsurvivaldata)
  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V, var_type = "binary", mod_type = "outcome")
  models <- list(mod1, mod3)

  fit <- gformula(
    data = nonsurvivaldata,
    id_var = "id",
    base_vars = "V",
    exposure = "A",
    time_var = "time",
    models = models,
    intervention = list(natural = NULL),
    mc_sample = 200,
    R = 0,
    quiet = TRUE,
    return_data = TRUE,
    return_fitted = TRUE
  )

  expect_true(inherits(fit$sim_data, "data.table"))
  expect_true("Pred_Y" %in% names(fit$sim_data))
  expect_true("Intervention" %in% names(fit$sim_data))
  expect_equal(unique(fit$sim_data$Intervention), "natural")

  # Check fitted models returned as objects
  expect_true(inherits(fit$fitted_models[[1]], "glm")) # spec_model by default uses glm logic inside probably, or we check if it has recodes attr
  expect_true(!is.null(attr(fit$fitted_models[[1]], "mod_type")))

})


testthat::test_that("Fixed intervention", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("gfoRmula")

  suppressPackageStartupMessages(library(gfoRmula))

  #data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, time := as.integer(time)]
  data.table::setorder(dat, id, time)

  time_points <- max(dat$time, na.rm = TRUE) + 1L

  # ----- recodes for lags (causalMed side) ---------------------------------
  init_recodes <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
  in_recodes   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)

  # ----- models (temporal order) -------------------------------------------
  m_L1 <- spec_model(
    L1 ~ lag1_A + lag1_L1 + V + time,
    var_type = "normal",   mod_type = "covariate"
  )
  m_L2 <- spec_model(
    L2 ~ lag1_A + lag1_L2 + V + time,
    var_type = "binary", mod_type = "covariate"
  )
  m_A <- spec_model(
    A  ~ lag1_A + L1 + L2 + V + time,
    var_type = "binary", mod_type = "exposure"
  )
  m_Y <- spec_model(
    Y_bin  ~ A + L1 + L2 ,
    var_type = "binary", mod_type = "outcome"
  )
  models <- list(m_L1, m_L2, m_A, m_Y)

  # ----- causalMed::gformula (natural vs always=1) --------------------------
  mc_nsim <- 50000L  # Monte Carlo size; balance accuracy vs test time
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id",
      time_var      = "time",
      base_vars     = "V",
      exposure      = "A",
      models        = models,
      intervention  = list(natural = NULL, always = 1),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim,
      R             = 1,           # point estimates only
      return_data   = FALSE,
      return_fitted = FALSE,
      quiet         = TRUE,
      seed          = 20250915
    )
  )

  es_cm <- data.table::as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])

  # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
  dat2<-dat
  dat2$time<-as.integer(dat2$time)
  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat2,
      id            = "id",
      time_name     = "time",
      time_points   = 5,
      covnames      = c('L1', 'L2', 'A'),
      covtypes      = c("normal", "binary", "binary"),
      basecovs      = c("V"),
      covparams     = list(covmodels = c(
        L1 ~ lag1_A + lag1_L1 + V + time,
        L2 ~ lag1_A + lag1_L2 + V + time,
        A  ~ lag1_A + L1 + L2 + V + time
      )),
      ymodel        = Y_bin ~ A + L1 + L2 ,
      outcome_name  = "Y_bin",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c('A', 'L1', 'L2')),
      intervention1.A = list(gfoRmula::static, rep(1, 5)),
      int_descript  = "Always treat",
      nsimul        = 50000,
      seed          = 20250915
    )
  )

  gfm <- res_gf$result[["g-form mean"]]
  means_gf <- as.numeric(gfm[1:2])

  # ----- Compare with tolerance --------------------------------------------
  testthat::expect_equal(means_cm, means_gf, tolerance = 0.005)
})


testthat::test_that("Predefined intervention", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("gfoRmula")  

  suppressPackageStartupMessages(library(gfoRmula))

  data("nonsurvivaldata", package = "causalMed")
  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, time := as.integer(time)]
  data.table::setorder(dat, id, time)

  time_points <- max(dat$time, na.rm = TRUE) + 1L

  # ----- recodes for lags (causalMed side) ---------------------------------
  init_recodes <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
  in_recodes   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)

  # ----- models (temporal order) -------------------------------------------
  m_L1 <- spec_model(
    L1 ~ lag1_A + lag1_L1 + V + time,
    var_type = "normal",   mod_type = "covariate"
  )
  m_L2 <- spec_model(
    L2 ~ lag1_A + lag1_L2 + V + time,
    var_type = "binary", mod_type = "covariate"
  )
  m_A <- spec_model(
    A  ~ lag1_A + L1 + L2 + V + time,
    var_type = "binary", mod_type = "exposure"
  )
  m_Y <- spec_model(
    Y_bin  ~ A + L1 + L2 ,
    var_type = "binary", mod_type = "outcome"
  )
  models <- list(m_L1, m_L2, m_A, m_Y)

  # ----- causalMed::gformula (natural vs always=1) --------------------------
  mc_nsim <- 50000L  # Monte Carlo size; balance accuracy vs test time
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id",
      time_var      = "time",
      base_vars     = "V",
      exposure      = "A",
      models        = models,
      intervention  = list(Predefined=c(0,0,1,1,1)),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim,
      R             = 1,           # point estimates only
      return_data   = FALSE,
      return_fitted = FALSE,
      quiet         = TRUE,
      seed          = 20250915
    )
  )
  es_cm <- data.table::as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])

  # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
  dat2<-dat
  dat2$time<-as.integer(dat2$time)
  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat2,
      id            = "id",
      time_name     = "time",
      time_points   = 5,
      covnames      = c('L1', 'L2', 'A'),
      covtypes      = c("normal", "binary", "binary"),
      basecovs      = c("V"),
      covparams     = list(covmodels = c(
        L1 ~ lag1_A + lag1_L1 + V + time,
        L2 ~ lag1_A + lag1_L2 + V + time,
        A  ~ lag1_A + L1 + L2 + V + time
      )),
      ymodel        = Y_bin ~ A + L1 + L2 ,
      outcome_name  = "Y_bin",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c('A', 'L1', 'L2')),
      intervention1.A = list(gfoRmula::static, c(0,0,1,1,1)),
      int_descript  = "Always treat",
      nsimul        = mc_nsim,
      seed          = 20250915
    )
  )

  gfm <- res_gf$result[["g-form mean"]]
  means_gf <- as.numeric(gfm[1:2])

  # ----- Compare with tolerance --------------------------------------------
  testthat::expect_equal(means_cm, means_gf, tolerance = 0.005)
})



 testthat::test_that("Dynamic intervention: causalMed vs gfoRmula", {
   testthat::skip_on_cran()
   testthat::skip_if_not_installed("gfoRmula")

   suppressPackageStartupMessages(library(data.table))
   suppressPackageStartupMessages(library(gfoRmula))

   dat <- data.table::as.data.table(nonsurvivaldata)
   dat[, time := as.integer(time)]
   data.table::setorder(dat, id, time)
   
   time_points <- max(dat$time, na.rm = TRUE) + 1L
   
   # ----- recodes for lags (causalMed side) ---------------------------------
   init_recodes <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
   in_recodes   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)
   
   # ----- models (temporal order) -------------------------------------------
   m_L1 <- spec_model(
     L1 ~ lag1_A + lag1_L1 + V + time,
     var_type = "normal",   mod_type = "covariate"
   )
   m_L2 <- spec_model(
     L2 ~ lag1_A + lag1_L2 + V + time,
     var_type = "binary", mod_type = "covariate"
   )
   m_A <- spec_model(
     A  ~ lag1_A + L1 + L2 + V + time,
     var_type = "binary", mod_type = "exposure"
   )
   m_Y <- spec_model(
     Y_bin  ~ A + L1 + L2 ,
     var_type = "binary", mod_type = "outcome"
   )
   models <- list(m_L1, m_L2, m_A, m_Y)
   
   # ----- causalMed::gformula (natural vs always=1) --------------------------
   mc_nsim <- 50000L  # Monte Carlo size; balance accuracy vs test time
   fit_cm <- testthat::expect_no_error(
     causalMed::gformula(
       data          = dat,
       id_var        = "id",
       time_var      = "time",
       base_vars     = "V",
       exposure      = "A",
       models        = models,
       intervention  = list(dynamic = dyn_int(as.numeric(A > 0))),
       init_recode   = init_recodes,
       in_recode     = in_recodes,
       out_recode    = NULL,
       mc_sample     = mc_nsim,
       R             = 1,           # point estimates only
       return_data   = FALSE,
       return_fitted = FALSE,
       quiet         = TRUE,
       seed          = 20250915
     )
   )
   
   es_cm <- data.table::as.data.table(fit_cm$effect_size)
   means_cm <- as.numeric(es_cm$Est[1:2])
   
   # ----- gfoRmula::gformula (natural vs A=1 static) -------------------------
   dat2<-dat
   dat2$time<-as.integer(dat2$time)
   res_gf <- testthat::expect_no_error(
     gfoRmula::gformula(
       obs_data      = dat2,
       id            = "id",
       time_name     = "time",
       time_points   = 5,
       covnames      = c('L1', 'L2', 'A'),
       covtypes      = c("normal", "binary", "binary"),
       basecovs      = c("V"),
       covparams     = list(covmodels = c(
         L1 ~ lag1_A + lag1_L1 + V + time,
         L2 ~ lag1_A + lag1_L2 + V + time,
         A  ~ lag1_A + L1 + L2 + V + time
       )),
       ymodel        = Y_bin ~ A + L1 + L2 ,
       outcome_name  = "Y_bin",
       outcome_type  = "binary_eof",
       histories     = c(gfoRmula::lagged),
       histvars      = list(c('A', 'L1', 'L2')),
       intervention1.A =  list(function(newdf, pool, intvar, intvals, time_name, t){
         newdf[, (intvar) := as.numeric(newdf[[intvar]] > 0)]
       }),
       int_descript = "Treat=1 if natural treat>0; else 0",
       nsimul        = 50000,
       seed          = 20250915
     )
   )
   
   gfm <- res_gf$result[["g-form mean"]]
   means_gf <- as.numeric(gfm[1:2])
   
   # ----- Compare with tolerance --------------------------------------------
   testthat::expect_equal(means_cm, means_gf, tolerance = 0.005)
 })


testthat::test_that("Large-sample cross-validation: causalMed vs gfoRmula (tol=0.005)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("gfoRmula")  

  suppressPackageStartupMessages(library(gfoRmula))

  dat <- data.table::as.data.table(nonsurvivaldata)
  dat[, time := as.integer(time)]
  data.table::setorder(dat, id, time)

  # ----- recodes for lags ---------------------------------------------------
  init_recodes <- recodes(lag1_A = 0, lag1_L1 = 0, lag1_L2 = 0)
  in_recodes   <- recodes(lag1_A = A, lag1_L1 = L1, lag1_L2 = L2)

  # ----- models (temporal order) --------------------------------------------
  m_L1 <- spec_model(
    L1 ~ lag1_A + lag1_L1 + V + time,
    var_type = "normal",  mod_type = "covariate"
  )
  m_L2 <- spec_model(
    L2 ~ lag1_A + lag1_L2 + V + time,
    var_type = "binary", mod_type = "covariate"
  )
  m_A <- spec_model(
    A  ~ lag1_A + L1 + L2 + V + time,
    var_type = "binary", mod_type = "exposure"
  )
  m_Y <- spec_model(
    Y_bin ~ A + L1 + L2,
    var_type = "binary", mod_type = "outcome"
  )
  models <- list(m_L1, m_L2, m_A, m_Y)

  mc_nsim <- 50000L  # large sample to keep Monte Carlo noise below the tolerance

  # ----- causalMed::gformula ------------------------------------------------
  fit_cm <- testthat::expect_no_error(
    causalMed::gformula(
      data          = dat,
      id_var        = "id",
      time_var      = "time",
      base_vars     = "V",
      exposure      = "A",
      models        = models,
      intervention  = list(natural = NULL, always = 1),
      init_recode   = init_recodes,
      in_recode     = in_recodes,
      out_recode    = NULL,
      mc_sample     = mc_nsim,
      R             = 1,
      return_data   = FALSE,
      return_fitted = FALSE,
      quiet         = TRUE,
      seed          = 20250915
    )
  )

  es_cm    <- data.table::as.data.table(fit_cm$effect_size)
  means_cm <- as.numeric(es_cm$Est[1:2])

  # ----- gfoRmula::gformula -------------------------------------------------
  dat2      <- data.table::copy(dat)
  dat2$time <- as.integer(dat2$time)

  res_gf <- testthat::expect_no_error(
    gfoRmula::gformula(
      obs_data      = dat2,
      id            = "id",
      time_name     = "time",
      time_points   = 5,
      covnames      = c("L1", "L2", "A"),
      covtypes      = c("normal", "binary", "binary"),
      basecovs      = c("V"),
      covparams     = list(covmodels = c(
        L1 ~ lag1_A + lag1_L1 + V + time,
        L2 ~ lag1_A + lag1_L2 + V + time,
        A  ~ lag1_A + L1 + L2 + V + time
      )),
      ymodel        = Y_bin ~ A + L1 + L2,
      outcome_name  = "Y_bin",
      outcome_type  = "binary_eof",
      histories     = c(gfoRmula::lagged),
      histvars      = list(c("A", "L1", "L2")),
      intervention1.A = list(gfoRmula::static, rep(1, 5)),
      int_descript  = "Always treat",
      nsimul        = mc_nsim,
      seed          = 20250915
    )
  )

  gfm      <- res_gf$result[["g-form mean"]]
  means_gf <- as.numeric(gfm[1:2])

  # ----- Compare with tolerance --------------------------------------------
  testthat::expect_equal(means_cm, means_gf, tolerance = 0.005)
})


test_that("gformula retains per-replicate bootstrap estimates and prints new fields", {
  data(nonsurvivaldata)

  mod1 <- spec_model(A ~ V, var_type = "binary", mod_type = "exposure")
  mod3 <- spec_model(Y_bin ~ A + V, var_type = "binary", mod_type = "outcome")
  models <- list(mod1, mod3)

  fit <- gformula(
    data = nonsurvivaldata,
    id_var = "id",
    base_vars = "V",
    exposure = "A",
    time_var = "time",
    models = models,
    intervention = list(natural = NULL, always = 1),
    ref_int = 0,
    mc_sample = 200,
    R = 5,
    quiet = TRUE
  )

  # boot_estimates: 5 replicates x 2 interventions; 5 x 2 contrast rows (RD + RR)
  expect_true(is.list(fit$boot_estimates))
  bi <- fit$boot_estimates$interventions
  expect_true(all(c("replicate", "Intervention", "Est") %in% names(bi)))
  expect_equal(nrow(bi), 5L * 2L)
  expect_equal(sort(unique(bi$replicate)), 1:5)
  bc <- fit$boot_estimates$contrasts
  expect_true(all(c("replicate", "Intervention", "Risk_type", "Estimate") %in% names(bc)))
  expect_equal(nrow(bc), 5L * 2L)

  # print shows the new setup fields, observed benchmark, and CI footnote
  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "Outcome      : Y_bin")
  expect_match(out, "Data         :")
  expect_match(out, "Observed (nonparametric)", fixed = TRUE)
  expect_match(out, "95% CIs", fixed = TRUE)
  # per-intervention mean CIs must NOT be labelled as RD
  expect_false(grepl("RD 2.5%", out, fixed = TRUE))
})
