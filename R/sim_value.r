
#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only, predict response and simulate random data. For numeric
#'  values the simulated value is restricted to the observed value range unless
#'  the model was created with \code{spec_model(truncate = FALSE)}.
#'
#' @param model fitted objects defined in the `spec_model`.
#' @param newdt a data frame in which to look for variables with which to predict.
#'
#' @return A simulated random vector using the predicted value from model and newdt.
#'
#' @keywords internal
#'
sim_value <- function(model, newdt) {
  var_type <- model$var_type

  # Clip simulated numeric values to the observed range of the response unless
  # the model was created with spec_model(truncate = FALSE). NULL means the
  # model object predates the `truncate` argument: keep the historical default.
  do_trunc <- is.null(model$truncate) || isTRUE(model$truncate)

  # Custom simulation function takes priority
  if (!is.null(model$custom_sim)) {
    out <- model$custom_sim(model$fitted, newdt)
    if (do_trunc && is.numeric(model$val_ran)) {
      out <- pmin(pmax(out, model$val_ran[1L]), model$val_ran[2L])
    }
    return(out)
  }

  # Categorical: multinom returns a probability matrix — predict() is required
  if (var_type == "categorical") {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      stop("Package 'Hmisc' is required for var_type = \"categorical\". ",
           "Please install it.", call. = FALSE)
    }
    pred <- predict(model$fitted, newdata = newdt, type = "probs")
    return(Hmisc::rMultinom(pred, 1))
  }

  # Fast linear predictor for binary and normal types.
  # model$Xterms and model$beta are pre-extracted once after fitting in .run_interventions,
  # avoiding the model.frame() + na.action overhead that predict() incurs.
  # model.matrix() + %*% is a direct BLAS call.
  mm <- model.matrix(model$Xterms, data = newdt)
  lp <- drop(mm %*% model$beta)

  if (var_type == "binary") {
    pred <- model$linkinv(lp)
    pred <- pmin(pmax(pred, 1e-5), 1 - 1e-5)
    return(rbinom(nrow(newdt), 1L, pred))
  } else {
    # Continuous/normal: lp + N(0, sigma) avoids a full predict() call.
    # model$sigma is the residual SD pre-extracted after fitting; it is NULL
    # when sigma() failed on the fitted object (e.g. a custom_fit class with
    # no sigma method and no custom_sim), where `NULL < eps` would give the
    # opaque "argument is of length zero" error.
    if (is.null(model$sigma)) {
      stop("No residual SD could be extracted from the fitted model for '",
           model$rsp_vars, "' (sigma() failed), so normal-type simulation ",
           "is not possible. Supply a custom_sim function in spec_model().",
           call. = FALSE)
    }
    if (model$sigma < .Machine$double.eps) {
      causalmed_env$warning <- c(
        causalmed_env$warning,
        sprintf("Zero residual variance in model for '%s'. Simulated values will equal the predicted mean.",
                model$rsp_vars)
      )
    }
    out <- lp + rnorm(nrow(newdt), 0, model$sigma)
    if (do_trunc) {
      out <- pmin(pmax(out, model$val_ran[1L]), model$val_ran[2L])
    }
    return(out)
  }
}
