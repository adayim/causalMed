
#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only, predict response and simulate random data. The simulated
#'  value will be restricted within the observed value range for numeric values.
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

  # Custom simulation function takes priority
  if (!is.null(model$custom_sim)) {
    out <- model$custom_sim(model$fitted, newdt)
    if (is.numeric(model$val_ran)) {
      out <- pmax(out, model$val_ran[1L])
      out <- pmin(out, model$val_ran[2L])
    }
    return(out)
  }

  # Categorical: multinom returns a probability matrix — predict() is required
  if (var_type == "categorical") {
    pred <- predict(model$fitted, newdata = newdt, type = "probs")
    return(Hmisc::rMultinom(pred, 1))
  }

  # Fast linear predictor for binary and normal types.
  # model$Xterms and model$beta are pre-extracted once after fitting in .gformula,
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
    # model$sigma is the residual SD pre-extracted after fitting.
    if (model$sigma < .Machine$double.eps) {
      warning("Zero residual variance in model for '",
              all.vars(formula(model$fitted)[[2]]),
              "'. Simulated values will equal the predicted mean.")
    }
    out <- lp + rnorm(nrow(newdt), 0, model$sigma)
    out <- pmax(out, model$val_ran[1L])
    out <- pmin(out, model$val_ran[2L])
    return(out)
  }
}
