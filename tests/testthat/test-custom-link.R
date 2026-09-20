# Without custom_sim, a "normal" or "custom" node is drawn as lp + N(0, sigma):
# centred on the linear predictor, which is the fitted mean only under an
# identity link. A fit with any other link is refused rather than simulated on
# the link scale.

testthat::test_that("a non-identity-link fit without custom_sim is refused", {
  data("nonsurvivaldata", package = "causalMed")
  d <- data.table::as.data.table(nonsurvivaldata)
  set.seed(1)
  d[, cnt := stats::rpois(.N, exp(0.2 + 0.3 * V))]
  outcome <- spec_model(Y_bin ~ V + cnt, var_type = "binary", mod_type = "outcome")

  # custom_fit's default glm() with a Poisson family, and no custom_sim.
  testthat::expect_error(
    causalMed:::fit_spec_models(list(
      spec_model(cnt ~ V, var_type = "custom", mod_type = "covariate",
                 family = stats::poisson()),
      outcome), d),
    "fitted with a 'log' link", fixed = TRUE)

  # The same fit with a custom_sim that draws from the Poisson is accepted.
  testthat::expect_no_error(causalMed:::fit_spec_models(list(
    spec_model(cnt ~ V, var_type = "custom", mod_type = "covariate",
               family = stats::poisson(),
               custom_sim = function(fit, newdf)
                 stats::rpois(nrow(newdf),
                              stats::predict(fit, newdata = newdf, type = "response"))),
    outcome), d))

  # An identity link (the default gaussian glm) is unaffected.
  testthat::expect_no_error(causalMed:::fit_spec_models(list(
    spec_model(cnt ~ V, var_type = "custom", mod_type = "covariate"),
    outcome), d))
})
